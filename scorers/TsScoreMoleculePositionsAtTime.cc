// Scorer for MoleculePositionsAtTime
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

#include "TsScoreMoleculePositionsAtTime.hh"

#include "G4SystemOfUnits.hh"
#include "G4Scheduler.hh"
#include "G4ITTrackHolder.hh"
#include "G4Molecule.hh"
#include "G4H2O.hh"
#include "G4EventManager.hh"
#include "G4Event.hh"

TsScoreMoleculePositionsAtTime::TsScoreMoleculePositionsAtTime(TsParameterManager* pM, TsMaterialManager* mM, TsGeometryManager* gM, TsScoringManager* scM, TsExtensionManager* eM,
                                                               G4String scorerName, G4String quantity, G4String outFileName, G4bool isSubScorer)
: TsVNtupleScorer(pM, mM, gM, scM, eM, scorerName, quantity, outFileName, isSubScorer), fRecorded(false)
{
    SetUnit("");

    fRecordTime = 1.0 * ps;
    if ( fPm->ParameterExists(GetFullParmName("RecordTime")) )
        fRecordTime = fPm->GetDoubleParameter(GetFullParmName("RecordTime"), "Time");

    fNtuple->RegisterColumnI(&fEvt, "EventID");
    fNtuple->RegisterColumnS(&fMoleculeName, "MoleculeName");
    fNtuple->RegisterColumnF(&fPosX, "Position X", "nm");
    fNtuple->RegisterColumnF(&fPosY, "Position Y", "nm");
    fNtuple->RegisterColumnF(&fPosZ, "Position Z", "nm");
    fNtuple->RegisterColumnI(&fTrackID, "Track ID");
    fNtuple->RegisterColumnI(&fParentID, "Parent track ID");
    fNtuple->RegisterColumnB(&fParentVertexKnown, "Parent vertex known");
    fNtuple->RegisterColumnF(&fParentVertexX, "Parent vertex X", "nm");
    fNtuple->RegisterColumnF(&fParentVertexY, "Parent vertex Y", "nm");
    fNtuple->RegisterColumnF(&fParentVertexZ, "Parent vertex Z", "nm");
    fNtuple->RegisterColumnF(&fTime, "Time recorded", "ps");
}


TsScoreMoleculePositionsAtTime::~TsScoreMoleculePositionsAtTime() {;}


G4bool TsScoreMoleculePositionsAtTime::ProcessHits(G4Step* aStep, G4TouchableHistory*)
{
    if (!fIsActive) {
        fSkippedWhileInactive++;
        return false;
    }
    // Physical tracks carry positive IDs; molecules carry negative ones. Keep every physical
    // track's vertex so a solvated electron can be tied back to where its parent was ejected.
    G4Track* aTrack = aStep->GetTrack();
    if ( aTrack->GetTrackID() > 0 )
        fVertexOfPhysicalTrack[aTrack->GetTrackID()] = aTrack->GetVertexPosition();
    return false;
}


void TsScoreMoleculePositionsAtTime::UserHookForPostTimeStepAction()
{
    if ( fRecorded ) return;
    G4double now = G4Scheduler::Instance()->GetGlobalTime();
    if ( now < fRecordTime ) return;
    // An aborted event (primary energy loss overshot the abort threshold) still runs its
    // chemistry, but on a truncated physical stage; the G-value scorer drops it, so drop it here.
    if ( G4EventManager::GetEventManager()->GetConstCurrentEvent()->IsAborted() ) { fRecorded = true; return; }

    G4TrackManyList* mainList = G4ITTrackHolder::Instance()->GetMainList();
    if ( !mainList ) return;

    // Water that has not yet dissociated is not a radical; wait one more step for it.
    for ( auto it = mainList->begin(); it != mainList->end(); ++it )
        if ( GetMolecule(*it)->GetDefinition() == G4H2O::Definition() ) return;

    // Values are stored in Geant4 internal units; TsVNtuple applies the registered unit itself.
    fEvt = GetEventID();
    fTime = now;
    for ( auto it = mainList->begin(); it != mainList->end(); ++it ) {
        G4Track* track = *it;
        G4ThreeVector pos = track->GetPosition();
        fMoleculeName = GetMolecule(track)->GetName();
        fPosX = pos.x();
        fPosY = pos.y();
        fPosZ = pos.z();
        fTrackID = track->GetTrackID();
        fParentID = track->GetParentID();
        auto v = fVertexOfPhysicalTrack.find(fParentID);
        fParentVertexKnown = ( v != fVertexOfPhysicalTrack.end() );
        fParentVertexX = fParentVertexKnown ? v->second.x() : 0.0;
        fParentVertexY = fParentVertexKnown ? v->second.y() : 0.0;
        fParentVertexZ = fParentVertexKnown ? v->second.z() : 0.0;
        fNtuple->Fill();
    }
    fRecorded = true;
}


void TsScoreMoleculePositionsAtTime::UserHookForEndOfEvent()
{
    fRecorded = false;
    fVertexOfPhysicalTrack.clear();
}
