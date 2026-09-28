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

#ifndef TsScoreMoleculePositionsAtTime_hh
#define TsScoreMoleculePositionsAtTime_hh

#include "TsVNtupleScorer.hh"
#include "G4ThreeVector.hh"

#include <map>

// One ntuple row per live molecule, taken once per event at the first chemical time step
// that reaches Sc/<name>/RecordTime (default 1 ps). This is the spatial pattern the physical
// stage hands to the chemistry, which the G-value scorer only ever sees as a count.
//
// The parent linkage: a solvated electron's G4Track parent is the sub-excitation electron it
// thermalised from, and that electron's vertex is the ionisation it was ejected at, which is
// where the parent cation sits. Physical tracks are gone by 1 ps, so their vertices are kept
// in a per-event map filled from ProcessHits during the physical stage.
class TsScoreMoleculePositionsAtTime : public TsVNtupleScorer
{
public:
    TsScoreMoleculePositionsAtTime(TsParameterManager* pM, TsMaterialManager* mM, TsGeometryManager* gM, TsScoringManager* scM, TsExtensionManager* eM,
                                   G4String scorerName, G4String quantity, G4String outFileName, G4bool isSubScorer);
    virtual ~TsScoreMoleculePositionsAtTime();

    virtual G4bool ProcessHits(G4Step*, G4TouchableHistory*);
    virtual void UserHookForPostTimeStepAction();
    virtual void UserHookForEndOfEvent();

private:
    G4double fRecordTime;
    G4bool   fRecorded;
    std::map<G4int, G4ThreeVector> fVertexOfPhysicalTrack;

    // ntuple columns
    G4int    fEvt;
    G4String fMoleculeName;
    G4float  fPosX, fPosY, fPosZ;
    G4int    fTrackID;
    G4int    fParentID;
    G4bool   fParentVertexKnown;
    G4float  fParentVertexX, fParentVertexY, fParentVertexZ;
    G4float  fTime;
};

#endif
