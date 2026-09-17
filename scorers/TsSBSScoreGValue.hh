//
// ********************************************************************
// *                                                                  *
// * This file is part of the TOPAS-nBio extensions to the            *
// *   TOPAS Simulation Toolkit.                                      *
// * The TOPAS-nBio extensions are freely available under the license *
// *   agreement set forth at: https://topas-nbio.readthedocs.io/     *
// *                                                                  *
// ********************************************************************

#ifndef TsSBSScoreGValue_hh
#define TsSBSScoreGValue_hh

#include "TsVNtupleScorer.hh"
#include "Randomize.hh"

#include <stdint.h>

class G4MolecularConfiguration;

class TsSBSScoreGValue : public TsVNtupleScorer
{
public:
    TsSBSScoreGValue(TsParameterManager* pM, TsMaterialManager* mM, TsGeometryManager* gM, TsScoringManager* scM, TsExtensionManager* eM,
                  G4String scorerName, G4String quantity, G4String outFileName, G4bool isSubScorer);
    ~TsSBSScoreGValue();
    
    virtual G4bool ProcessHits(G4Step*, G4TouchableHistory*);
    void AbsorbResultsFromWorkerScorer(TsVScorer* workerScorer);
    void UserHookForEndOfEvent();
   
protected:
//    void AccumulateEvent();
    
    void Output();
    void Clear();
    
    // Output variables
    G4double fGValue;
    G4double fGValueError;
    G4double fTime;
    G4String fMoleculeName;
    G4double fGValueRatio;
    G4double fGValueRatioError;

    std::map<G4String, std::map<G4double, G4double> > fGValuePerSpeciePerTime;
    std::map<G4String, std::map<G4double, G4double> > fGValuePerSpeciePerTime2;

    // Sums for the ratio-of-sums estimator, G = 100 * sum(N) / sum(E). The per-event ratio the
    // columns above report is a different estimator of the same quantity and a worse one: a G
    // value is molecules per 100 eV ABSORBED, so the energy belongs in the denominator of one
    // ratio taken over the whole run, not inside a per-event ratio that is then averaged. Because
    // E varies event to event through straggling, averaging n/E over-weights the low-energy
    // events, and a single event that deposits almost nothing while making one molecule sends the
    // estimate to infinity. Seen in practice: 8 seeds of the shipped 300-history example gave
    // G(H2O2) at 1 us of 0.40 to 0.48 on six of them and 8.96 on one.
    std::map<G4String, std::map<G4double, G4double> > fMoleculesPerSpeciePerTime;   // sum n
    std::map<G4String, std::map<G4double, G4double> > fMoleculesPerSpeciePerTime2;  // sum n^2
    std::map<G4String, std::map<G4double, G4double> > fMoleculesTimesEnergy;        // sum n*E

private:
    TsParameterManager* fPm;
    
    G4double fEnergyDepositPerEvent;
    G4double fSumEnergy;    // sum of per-event deposited energy, eV
    G4double fSumEnergy2;   // sum of its square
    G4double* fTimeToRecord;
    G4int fNbTimeToRecord;
    G4int fNbOfScoredEvents;
    G4double fEnergyLossKill;
    G4double fEnergyLossAbort;
    G4double fEnergyLoss;
    G4double fMaximumTrackLength;
    G4double fTotalTrackLength;
    G4int fNbOfMoleculesToScavenge;
    G4int* fMoleculeIDToScavenge;
    G4double* fScavengingCapacity;

    std::map<G4double, G4double> fScavengerProducts;
};

#endif

