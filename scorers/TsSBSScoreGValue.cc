// Scorer for TsSBSGValue
//
// ********************************************************************
// *                                                                  *
// * This file is part of the TOPAS-nBio extensions to the            *
// *   TOPAS Simulation Toolkit.                                      *
// * The TOPAS-nBio extensions are freely available under the license *
// *   agreement set forth at: https://topas-nbio.readthedocs.io/     *
// *                                                                  *
// ********************************************************************

#include "TsSBSScoreGValue.hh"

#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4DNAChemistryManager.hh"
#include "G4Scheduler.hh"
#include "G4MoleculeCounter.hh"
#include "G4MolecularConfiguration.hh"
#include "G4Molecule.hh"
#include "G4H2O.hh"

TsSBSScoreGValue::TsSBSScoreGValue(TsParameterManager* pM, TsMaterialManager* mM, TsGeometryManager* gM, TsScoringManager* scM, TsExtensionManager* eM,
                             G4String scorerName, G4String quantity, G4String outFileName, G4bool isSubScorer)
: TsVNtupleScorer(pM, mM, gM, scM, eM, scorerName, quantity, outFileName, isSubScorer),
fPm(pM), fEnergyDepositPerEvent(0), fSumEnergy(0), fSumEnergy2(0)
{
    if ( G4Threading::IsWorkerThread() )  {
        G4VMoleculeCounter::Instance()->Use(true);
        G4MoleculeCounter::Instance()->DontRegister(G4H2O::Definition());
    } else {
        G4DNAChemistryManager::Instance()->ResetCounterWhenRunEnds(false);
    }
    
    fNtuple->RegisterColumnD(&fGValue, "GValue: number of molecules per energy deposit per eV", "");
    fNtuple->RegisterColumnD(&fGValueError, "GValue statistical error", "");
    fNtuple->RegisterColumnD(&fTime,    "Time", "picosecond");
    fNtuple->RegisterColumnS(&fMoleculeName, "MoleculeName");
    // Appended after the original four columns so existing parsers keep working unchanged.
    fNtuple->RegisterColumnD(&fGValueRatio, "GValue ratio-of-sums: 100*sum(N)/sum(E in eV)", "");
    fNtuple->RegisterColumnD(&fGValueRatioError, "GValue ratio-of-sums statistical error", "");


    fTimeToRecord = fPm->GetDoubleVector(GetFullParmName("TimeToRecord"),"Time");
    fNbTimeToRecord = fPm->GetVectorLength(GetFullParmName("TimeToRecord"));
    
    fEnergyLossKill = fPm->GetDoubleParameter(GetFullParmName("KillPrimaryIfEnergyLossExceeds"),"Energy");
    fEnergyLossAbort = fPm->GetDoubleParameter(GetFullParmName("AbortEventIfPrimaryEnergyLossExceeds"), "Energy");
    fMaximumTrackLength = fPm->GetDoubleParameter(GetFullParmName("KillPrimaryBasedOnTrackLength"), "Length");
    
    fNbOfMoleculesToScavenge = -1;
    if (fPm->ParameterExists(GetFullParmName("ScavengeTheseMolecules")) )  {
        fNbOfMoleculesToScavenge = fPm->GetVectorLength(GetFullParmName("ScavengeTheseMolecules"));
        fMoleculeIDToScavenge    = fPm->GetIntegerVector(GetFullParmName("ScavengeTheseMolecules"));
        fScavengingCapacity      = fPm->GetDoubleVector(GetFullParmName("ScavengingCapacities"), "perTime");
    }
    
    fNbOfScoredEvents = 0;
    fEnergyLoss = 0.0;
    fTotalTrackLength = 0.0;
    
}

TsSBSScoreGValue::~TsSBSScoreGValue()
{;}


G4bool TsSBSScoreGValue::ProcessHits(G4Step* aStep, G4TouchableHistory*)
{
    if (!fIsActive) {
        fSkippedWhileInactive++;
        return false;
    }
    
    G4Track* aTrack = aStep->GetTrack();
    if ( -1 < aTrack->GetTrackID()) {
        if ( aTrack->GetParentID() == 0 ) {
            G4double kineticE = aStep->GetPostStepPoint()->GetKineticEnergy();
            G4double eLoss = aStep->GetPreStepPoint()->GetKineticEnergy() - kineticE;
            fEnergyLoss += eLoss;
            fTotalTrackLength += aStep->GetStepLength();
            
            if ( fEnergyLoss > fEnergyLossAbort ) {
                std::cout << " -- Aborting event " << GetEventID() << std::endl;
                G4RunManager::GetRunManager()->AbortEvent();
                return true;
            }
            if ( fEnergyLoss >= fEnergyLossKill ) {
                std::cout << " -- Killing primary track of event " << GetEventID() << std::endl;
                aStep->GetTrack()->SetTrackStatus(fStopAndKill);
                return true;
            }
            if ( fTotalTrackLength >= fMaximumTrackLength ) {
                aStep->GetTrack()->SetTrackStatus(fStopAndKill);
                std::cout << " Track killed with track: " << fTotalTrackLength/nm << " nm with energy lost " << fEnergyLoss/keV << " keV "  << std::endl;
                return true;
            }
        }
        
        G4double edep = aStep->GetTotalEnergyDeposit();
        if ( edep > 0 ) {
            fEnergyDepositPerEvent += edep * aStep->GetPreStepPoint()->GetWeight();
            return true;
        }
    } else {
        if ( 0 < fNbOfMoleculesToScavenge ) {
            G4Molecule* molecule = GetMolecule(aTrack);
            G4int moleculeID = molecule->GetMoleculeID();
            // The survival probability against a first-order sink of capacity S must be
            // evaluated over THIS STEP's duration. Using the track's absolute global time
            // here gives 1 - exp(-S t), the cumulative probability of having been scavenged
            // by time t, and draws against it afresh at every step, so a molecule is offered
            // its whole accumulated chance of dying over and over. The effective removal rate
            // is then far above S, and increasingly so the longer the molecule lives, which
            // bends exactly the late-time part of a G value that a scavenger measurement is
            // meant to probe.
            G4double dt = aStep->GetDeltaTime();
            for ( int i = 0; i < fNbOfMoleculesToScavenge; i++ ) {
                if ( moleculeID == fMoleculeIDToScavenge[i] ) {
                    G4double probability = 1. - std::exp( -fScavengingCapacity[i] * dt );
                    if ( G4UniformRand() < probability ) {
                        aStep->GetTrack()->SetTrackStatus(fStopAndKill);
                        return true;
                    }
                }
            }
        }
    }
    
    return false;
}


void TsSBSScoreGValue::UserHookForEndOfEvent() {
    
    if (!G4EventManager::GetEventManager()->GetConstCurrentEvent()->IsAborted() && fEnergyDepositPerEvent > 0 ) {
        G4MoleculeCounter::RecordedMolecules species = G4MoleculeCounter::Instance()->GetRecordedMolecules();
        if ( species.get() == 0 ) return;
        
        G4double energyInEV = fEnergyDepositPerEvent / eV;

        auto molecule_it = species->begin();
        while ( molecule_it != species->end() ) {
            for ( int iTime = 0; iTime < fNbTimeToRecord; iTime++ ) {
                G4int moleculeAtSpecificTime = G4MoleculeCounter::Instance()->GetNMoleculesAtTime((*molecule_it), fTimeToRecord[iTime]);
                G4String moleculeID = (*molecule_it)->GetName();
                G4double gvalue = 100.0 * moleculeAtSpecificTime/energyInEV;
                fGValuePerSpeciePerTime[moleculeID][fTimeToRecord[iTime]] += gvalue;
                fGValuePerSpeciePerTime2[moleculeID][fTimeToRecord[iTime]] += gvalue * gvalue;

                G4double n = moleculeAtSpecificTime;
                fMoleculesPerSpeciePerTime[moleculeID][fTimeToRecord[iTime]] += n;
                fMoleculesPerSpeciePerTime2[moleculeID][fTimeToRecord[iTime]] += n * n;
                fMoleculesTimesEnergy[moleculeID][fTimeToRecord[iTime]] += n * energyInEV;
            }
            molecule_it++;
        }

        fSumEnergy += energyInEV;
        fSumEnergy2 += energyInEV * energyInEV;
        fNbOfScoredEvents++;
    }
    fTotalTrackLength = 0.0;
    fEnergyDepositPerEvent = 0.0;
    fEnergyLoss = 0.0;
    G4MoleculeCounter::Instance()->ResetCounter();
    
}


// Merge a worker's species/time sums into the master's. Unlike the two hand-written loops below
// this does not assume the master already holds every key the worker holds.
static void MergeSpeciesTimeMap(std::map<G4String, std::map<G4double, G4double> >& master,
                                std::map<G4String, std::map<G4double, G4double> >& worker) {
    for ( auto& speciesAndTimes : worker )
        for ( auto& timeAndValue : speciesAndTimes.second )
            master[speciesAndTimes.first][timeAndValue.first] += timeAndValue.second;
    worker.clear();
}


void TsSBSScoreGValue::AbsorbResultsFromWorkerScorer(TsVScorer* workerScorer) {
    TsVNtupleScorer::AbsorbResultsFromWorkerScorer(workerScorer);

    TsSBSScoreGValue* workerGvalueScorer = dynamic_cast<TsSBSScoreGValue*>(workerScorer);

    fNbOfScoredEvents += workerGvalueScorer->fNbOfScoredEvents;

    MergeSpeciesTimeMap(fMoleculesPerSpeciePerTime,  workerGvalueScorer->fMoleculesPerSpeciePerTime);
    MergeSpeciesTimeMap(fMoleculesPerSpeciePerTime2, workerGvalueScorer->fMoleculesPerSpeciePerTime2);
    MergeSpeciesTimeMap(fMoleculesTimesEnergy,       workerGvalueScorer->fMoleculesTimesEnergy);
    fSumEnergy  += workerGvalueScorer->fSumEnergy;
    fSumEnergy2 += workerGvalueScorer->fSumEnergy2;
    workerGvalueScorer->fSumEnergy = 0.0;
    workerGvalueScorer->fSumEnergy2 = 0.0;


    std::map<G4String, std::map<G4double, G4double> >::iterator wIter;
    std::map<G4String, std::map<G4double, G4double> >::iterator mIter;
    
    for ( wIter = workerGvalueScorer->fGValuePerSpeciePerTime.begin(); wIter != workerGvalueScorer->fGValuePerSpeciePerTime.end(); wIter++) {
        mIter = fGValuePerSpeciePerTime.find( wIter->first );
        if ( mIter == fGValuePerSpeciePerTime.end() ) {
            fGValuePerSpeciePerTime.insert(std::pair<G4String, std::map<G4double, G4double> > ( wIter->first, wIter->second));
        } else {
            std::map<G4double, G4double>::iterator witer;
            std::map<G4double, G4double>::iterator miter;
            for ( witer = (wIter->second).begin(); witer != (wIter->second).end(); witer++ ) {
                miter = (mIter->second).find( witer->first );
                miter->second += witer->second;
            }
        }
    }
    
    std::map<G4String, std::map<G4double, G4double> >::iterator wIter2;
    std::map<G4String, std::map<G4double, G4double> >::iterator mIter2;
    
    for ( wIter2 = workerGvalueScorer->fGValuePerSpeciePerTime2.begin(); wIter2 != workerGvalueScorer->fGValuePerSpeciePerTime2.end(); wIter2++) {
        mIter2 = fGValuePerSpeciePerTime2.find( wIter2->first );
        if ( mIter2 == fGValuePerSpeciePerTime2.end() ) {
            fGValuePerSpeciePerTime2.insert(std::pair<G4String, std::map<G4double, G4double> > ( wIter2->first, wIter2->second));
        } else {
            std::map<G4double, G4double>::iterator witer2;
            std::map<G4double, G4double>::iterator miter2;
            for ( witer2 = (wIter2->second).begin(); witer2 != (wIter2->second).end(); witer2++ ) {
                miter2 = (mIter2->second).find( witer2->first );
                miter2->second += witer2->second;
            }
        }
    }
    
    workerGvalueScorer->fGValuePerSpeciePerTime.clear();
    workerGvalueScorer->fGValuePerSpeciePerTime2.clear();
    workerGvalueScorer->fNbOfScoredEvents = 0;
    workerGvalueScorer->fEnergyDepositPerEvent = 0.0;
    
}


void TsSBSScoreGValue::Output() {
    
    if ( fNbOfScoredEvents == 0 )
        return;
    
    std::map<G4String, std::map<G4double, G4double> >::iterator wIter;
    std::map<G4String, std::map<G4double, G4double> >::iterator wIter2;
    std::map<G4double, G4double>::iterator iter;
    std::map<G4double, G4double>::iterator iter2;
    
    for ( wIter = fGValuePerSpeciePerTime.begin(); wIter != fGValuePerSpeciePerTime.end(); wIter++ ) {
        wIter2 = fGValuePerSpeciePerTime2.find(wIter->first);
        
        for ( iter = (wIter->second).begin(); iter != (wIter->second).end(); iter++ ) {
            iter2 = (wIter2->second).find( iter->first );
            
            fGValue = iter->second/fNbOfScoredEvents;
            if ( fNbOfScoredEvents > 1 ) {
                fGValueError = sqrt( (1.0/(fNbOfScoredEvents-1)) * ( (iter2->second)/fNbOfScoredEvents - fGValue*fGValue));
            } else {
                fGValueError = 1.0;
            }
            fTime = iter->first;
            fMoleculeName = wIter->first;

            // Ratio-of-sums estimator and its uncertainty. The variance is the standard delta
            // method for a ratio of two correlated means, which is what makes this estimator so
            // much steadier than the per-event one: n and E are strongly positively correlated,
            // so the covariance term cancels most of the straggling that the per-event ratio
            // leaves in.
            fGValueRatio = 0.0;
            fGValueRatioError = 0.0;
            if ( fSumEnergy > 0 ) {
                G4double n = fNbOfScoredEvents;
                G4double sumN = fMoleculesPerSpeciePerTime[fMoleculeName][fTime];
                G4double sumN2 = fMoleculesPerSpeciePerTime2[fMoleculeName][fTime];
                G4double sumNE = fMoleculesTimesEnergy[fMoleculeName][fTime];
                G4double meanN = sumN / n;
                G4double meanE = fSumEnergy / n;
                fGValueRatio = 100.0 * sumN / fSumEnergy;
                if ( n > 1 && meanN > 0 ) {
                    G4double varN = (sumN2 - n*meanN*meanN) / (n - 1);
                    G4double varE = (fSumEnergy2 - n*meanE*meanE) / (n - 1);
                    G4double covNE = (sumNE - n*meanN*meanE) / (n - 1);
                    G4double relVar = varN/(meanN*meanN) + varE/(meanE*meanE)
                                      - 2.0*covNE/(meanN*meanE);
                    if ( relVar > 0 )
                        fGValueRatioError = fGValueRatio * sqrt(relVar / n);
                }
            }
            fNtuple->Fill();
        }
    }
    
    fNtuple->Write();
}


void TsSBSScoreGValue::Clear() {
    fGValuePerSpeciePerTime.clear();
    fGValuePerSpeciePerTime2.clear();
    fMoleculesPerSpeciePerTime.clear();
    fMoleculesPerSpeciePerTime2.clear();
    fMoleculesTimesEnergy.clear();
    fSumEnergy = 0.0;
    fSumEnergy2 = 0.0;
    fNbOfScoredEvents = 0;
    UpdateFileNameForUpcomingRun();
}

