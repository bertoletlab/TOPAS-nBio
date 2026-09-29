// Extra Class for TsEmDNAChemistry

#include "TsDNAMaxTimeStep.hh"

#include "G4Molecule.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4VSolid.hh"
#include "G4Box.hh"
#include "G4VisExtent.hh"
#include "G4RotationMatrix.hh"

#include <algorithm>
#include <cmath>

#ifndef State
#define State(theXInfo) (GetState<MaxTimeStepState>()->theXInfo)
#endif

void TsDNAMaxTimeStep::Create()
{
    pParticleChange = &fParticleChange;
    enableAtRestDoIt    = false;
    enableAlongStepDoIt = false;
    enablePostStepDoIt  = true;

    SetProcessSubType(61);

    G4VITProcess::SetInstantiateProcessState(false);

    fProposesTimeStep = true;
    fMaxTimeStep = DBL_MAX;
    fVolumePrefix = "";
    fMargin = 0.;
    fRegionsCollected = false;

    verboseLevel = 0;
}

TsDNAMaxTimeStep::TsDNAMaxTimeStep(const G4String& aName, G4ProcessType type)
: G4VITProcess(aName, type)
{
    Create();
}

TsDNAMaxTimeStep::TsDNAMaxTimeStep(const TsDNAMaxTimeStep& rhs)
: G4VITProcess(rhs)
{
    Create();
    fMaxTimeStep = rhs.fMaxTimeStep;
    fVolumePrefix = rhs.fVolumePrefix;
    fMargin = rhs.fMargin;
}

TsDNAMaxTimeStep::~TsDNAMaxTimeStep()
{;}

TsDNAMaxTimeStep& TsDNAMaxTimeStep::operator=(const TsDNAMaxTimeStep& rhs)
{
    if (this == &rhs) return *this;
    fMaxTimeStep = rhs.fMaxTimeStep;
    fVolumePrefix = rhs.fVolumePrefix;
    fMargin = rhs.fMargin;
    fRegionsCollected = false;
    fRegions.clear();
    return *this;
}

void TsDNAMaxTimeStep::BuildPhysicsTable(const G4ParticleDefinition&)
{;}

void TsDNAMaxTimeStep::StartTracking(G4Track* track)
{
    G4VProcess::StartTracking(track);
    G4VITProcess::fpState.reset(new MaxTimeStepState());
    G4VITProcess::StartTracking(track);
}

// Collected on first use, once the geometry is closed. Positions are the placement
// translations, so the mother of the matching volumes is assumed to sit at the world origin
// without rotation, which is how the plasmid decks place their Solution box.
void TsDNAMaxTimeStep::CollectRegions()
{
    fRegionsCollected = true;
    fRegions.clear();
    if (fVolumePrefix.empty()) return;

    G4PhysicalVolumeStore* store = G4PhysicalVolumeStore::GetInstance();
    for (auto it = store->begin(); it != store->end(); ++it) {
        G4VPhysicalVolume* pv = *it;
        // TOPAS names a component's envelope after the component ("Plasmid0") and every
        // sub-volume with a slash after it ("Plasmid0/Base1_"), so the prefix alone would match
        // each DNA sphere: 208136 of them on an 8-plasmid lattice, placed in their envelope's
        // frame, not the world's. Envelopes only.
        const G4String& name = pv->GetName();
        if (name.compare(0, fVolumePrefix.size(), fVolumePrefix) != 0) continue;
        if (name.find('/') != std::string::npos) continue;

        G4VSolid* solid = pv->GetLogicalVolume()->GetSolid();
        Region r;
        G4ThreeVector offset(0., 0., 0.);
        if (G4Box* box = dynamic_cast<G4Box*>(solid)) {
            r.halfLengths = G4ThreeVector(box->GetXHalfLength(), box->GetYHalfLength(), box->GetZHalfLength());
        } else {
            G4VisExtent e = solid->GetExtent();
            r.halfLengths = G4ThreeVector(0.5 * (e.GetXmax() - e.GetXmin()),
                                          0.5 * (e.GetYmax() - e.GetYmin()),
                                          0.5 * (e.GetZmax() - e.GetZmin()));
            offset = G4ThreeVector(0.5 * (e.GetXmax() + e.GetXmin()),
                                   0.5 * (e.GetYmax() + e.GetYmin()),
                                   0.5 * (e.GetZmax() + e.GetZmin()));
        }
        r.centre = pv->GetObjectTranslation() + offset;
        r.radius = r.halfLengths.mag();
        r.isSphere = !pv->GetObjectRotationValue().isIdentity();
        fRegions.push_back(r);
    }

    G4cout << "-- TsDNAMaxTimeStep: cap " << G4BestUnit(fMaxTimeStep, "Time")
           << " inside " << fRegions.size() << " volume(s) named " << fVolumePrefix
           << "*, grown by " << G4BestUnit(fMargin, "Length")
           << "; outside, d^2 / (8 D) to the nearest of them." << G4endl;
}

G4double TsDNAMaxTimeStep::DistanceToNearestRegion(const G4ThreeVector& p) const
{
    G4double best = DBL_MAX;
    for (const Region& r : fRegions) {
        G4double d;
        if (r.isSphere) {
            d = (p - r.centre).mag() - r.radius;
        } else {
            G4ThreeVector q = p - r.centre;
            G4double qx = std::max(std::fabs(q.x()) - r.halfLengths.x(), 0.);
            G4double qy = std::max(std::fabs(q.y()) - r.halfLengths.y(), 0.);
            G4double qz = std::max(std::fabs(q.z()) - r.halfLengths.z(), 0.);
            d = std::sqrt(qx * qx + qy * qy + qz * qz);
        }
        if (d < best) best = d;
    }
    return best;
}

G4double TsDNAMaxTimeStep::PostStepGetPhysicalInteractionLength(const G4Track& track,
                                                                G4double,
                                                                G4ForceCondition* pForceCond)
{
    *pForceCond = NotForced;

    G4Molecule* mol = GetMolecule(track);
    if (!mol) return DBL_MAX;

    if (!fRegionsCollected) CollectRegions();

    G4double dt = fMaxTimeStep;
    if (!fRegions.empty()) {
        G4double d = DistanceToNearestRegion(track.GetPosition()) - fMargin;
        if (d > 0.) {
            G4double D = mol->GetDiffusionCoefficient();
            if (D > 0.) dt = std::max(fMaxTimeStep, d * d / (8. * D));
        }
    }

    // A time proposal, in the convention TsDNAFirstOrderReaction uses: the value is returned
    // negated, and the step processor reads it as a time rather than a length. It is proposed
    // afresh at every step, so a shorter step taken for any other reason costs nothing here.
    fpState->currentInteractionLength = dt;
    fpState->theNumberOfInteractionLengthLeft = 1.0;
    return -dt;
}

G4VParticleChange* TsDNAMaxTimeStep::PostStepDoIt(const G4Track& track, const G4Step&)
{
    fParticleChange.Initialize(track);
    return &fParticleChange;
}
