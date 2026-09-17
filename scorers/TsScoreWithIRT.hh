//
// ********************************************************************
// *                                                                  *
// * This file is part of the TOPAS-nBio extensions to the            *
// *   TOPAS Simulation Toolkit.                                      *
// * The TOPAS-nBio extensions are freely available under the license *
// *   agreement set forth at: https://topas-nbio.readthedocs.io/     *
// *                                                                  *
// ********************************************************************
//

#ifndef TsScoreWithIRT_hh
#define TsScoreWithIRT_hh

#include "TsVNtupleScorer.hh"
#include "TsIRTConfiguration.hh"

#include "G4Timer.hh"
#include "Randomize.hh"

#include <stdint.h>

class TsIRTManager;

class TsScoreWithIRT : public TsVNtupleScorer
{
public:
    TsScoreWithIRT(TsParameterManager* pM, TsMaterialManager* mM, TsGeometryManager* gM, TsScoringManager* scM, TsExtensionManager* eM,
                      G4String scorerName, G4String quantity, G4String outFileName, G4bool isSubScorer);
    ~TsScoreWithIRT();
    
    virtual G4bool ProcessHits(G4Step*, G4TouchableHistory*);
    void AbsorbResultsFromWorkerScorer(TsVScorer* workerScorer);
 
    void UserHookForEndOfEvent();
	
	virtual void UserHookForPreTimeStepAction();
	//virtual void UserHookForBeginOfChemicalTrack();

	//virtual void FillIRT(const G4Track*);
	
protected:
    
    void Output();
    void Clear();
    
    // Output variables
    G4double fGValue;
    G4double fGValueError;
    G4double fGValueRatio;
    G4double fGValueRatioError;
    G4double fTime;
    G4String fMoleculeName;
    G4double fMolecules;
    G4double fMoleculesError;
	
    std::map<G4String, std::map<G4double, G4double> > fGValuePerSpeciePerTime;
    std::map<G4String, std::map<G4double, G4double> > fGValuePerSpeciePerTime2;

    std::map<G4String, std::map<G4double, G4double> > fMoleculesPerSpeciePerTime;
    std::map<G4String, std::map<G4double, G4double> > fMoleculesPerSpeciePerTime2;

    // Sums for the ratio-of-sums estimator, G = 100 * sum(N) / sum(E). See the same block in
    // TsSBSScoreGValue.hh: sum(N) is already kept above, so only the energy sums and the N*E
    // cross term are new. Kept identical in both scorers so an IRT-vs-step-by-step comparison
    // is never a comparison of two different estimators.
    std::map<G4String, std::map<G4double, G4double> > fMoleculesTimesEnergy;

    std::map<G4int, std::map<G4double, G4double> > fDeltaGPerReactionPerTime;
    std::map<G4int, std::map<G4double, G4double> > fDeltaGPerReactionPerTime2;
   
private:
    TsParameterManager* fPm;
	TsIRTManager* fIRT;

	std::vector<TsIRTConfiguration::TsMolecule> fSpecies;

    G4Timer fIRTTimer[1];
    G4Timer fGilTimer[1];
    G4double fIRTExecutionTime;
    G4double fIRTExecutionTimeStdv;
	
    G4double fEnergyDepositPerEvent;
    G4double fSumEnergy;    // sum of per-event deposited energy, eV
    G4double fSumEnergy2;   // sum of its square
    G4int fNbOfScoredEvents;
    G4double fEnergyLossKill;
    G4double fEnergyLossAbort;
    G4double fEnergyLoss;
	G4String fName;
	
    G4double fTotalTrackLength;
    G4double fMaximumTrackLength;
    
    G4int  fMaximumNumberOfSteps;
    G4bool fUseMaximumNumberOfSteps;
    G4int  fNumberOfSteps;
    
    G4double fEnergyDepositPlusEnergyKinetic;
    G4double fLET;
    G4double fVolume;
		
	std::vector<G4double> fVEnergyDeposit;

    G4bool   fReportDelta;
    G4String fOutputFile;
    
    G4bool fTestIsInside;
    G4String fSensitiveVolume;
};

#endif

