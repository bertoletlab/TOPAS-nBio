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
    fNearDNAMinStep = 0.;
    fDNARadius = 0.;
    fCell = 2.0 * CLHEP::nanometer;

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
    fNearDNAMinStep = rhs.fNearDNAMinStep;
    fDNARadius = rhs.fDNARadius;
}

TsDNAMaxTimeStep::~TsDNAMaxTimeStep()
{;}

TsDNAMaxTimeStep& TsDNAMaxTimeStep::operator=(const TsDNAMaxTimeStep& rhs)
{
    if (this == &rhs) return *this;
    fMaxTimeStep = rhs.fMaxTimeStep;
    fVolumePrefix = rhs.fVolumePrefix;
    fMargin = rhs.fMargin;
    fNearDNAMinStep = rhs.fNearDNAMinStep;
    fDNARadius = rhs.fDNARadius;
    fRegionsCollected = false;
    fRegions.clear();
    fBasePairs.clear();
    fGrid.clear();
    return *this;
}

long long TsDNAMaxTimeStep::CellKey(const G4ThreeVector& p, int dx, int dy, int dz) const
{
    const long long ix = (long long)std::floor(p.x() / fCell) + dx;
    const long long iy = (long long)std::floor(p.y() / fCell) + dy;
    const long long iz = (long long)std::floor(p.z() / fCell) + dz;
    // 21 bits per axis, offset to keep them positive: +-1e6 cells of 2 nm is +-2 mm
    return ((ix + 1048576LL) << 42) | ((iy + 1048576LL) << 21) | (iz + 1048576LL);
}

G4double TsDNAMaxTimeStep::DistanceSquaredToNearestBasePair(const G4ThreeVector& p) const
{
    G4double best = DBL_MAX;
    for (int dx = -1; dx <= 1; dx++)
        for (int dy = -1; dy <= 1; dy++)
            for (int dz = -1; dz <= 1; dz++) {
                auto it = fGrid.find(CellKey(p, dx, dy, dz));
                if (it == fGrid.end()) continue;
                for (int idx : it->second) {
                    const G4double d2 = (fBasePairs[idx] - p).mag2();
                    if (d2 < best) best = d2;
                }
            }
    return best;
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

        // Base-pair positions for the near-DNA refinement: one daughter per base pair is
        // enough, and Base1_ is placed by every DNA model at the base-pair's own frame.
        if (fNearDNAMinStep > 0.) {
            G4LogicalVolume* lv = pv->GetLogicalVolume();
            const G4RotationMatrix rot = pv->GetObjectRotationValue();
            const G4ThreeVector trans = pv->GetObjectTranslation();
            for (size_t i = 0; i < lv->GetNoDaughters(); i++) {
                G4VPhysicalVolume* d = lv->GetDaughter(i);
                if (d->GetName().find("Base1_") == std::string::npos) continue;
                fBasePairs.push_back(rot * d->GetObjectTranslation() + trans);
            }
        }
    }
    for (size_t i = 0; i < fBasePairs.size(); i++)
        fGrid[CellKey(fBasePairs[i])].push_back((int)i);

    G4cout << "-- TsDNAMaxTimeStep: cap " << G4BestUnit(fMaxTimeStep, "Time")
           << " inside " << fRegions.size() << " volume(s) named " << fVolumePrefix
           << "*, grown by " << G4BestUnit(fMargin, "Length")
           << "; outside, d^2 / (8 D) to the nearest of them." << G4endl;
    if (fNearDNAMinStep > 0.)
        G4cout << "-- TsDNAMaxTimeStep: near-DNA refinement, " << fBasePairs.size()
               << " base-pair positions in a " << G4BestUnit(fCell, "Length") << " grid; a molecule "
               << "within reach of one steps by (d - " << G4BestUnit(fDNARadius, "Length")
               << ")^2 / (8 D), not below " << G4BestUnit(fNearDNAMinStep, "Time") << G4endl;
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
        const G4ThreeVector& p = track.GetPosition();
        G4double D = mol->GetDiffusionCoefficient();
        G4double d = DistanceToNearestRegion(p) - fMargin;
        if (d > 0.) {
            if (D > 0.) dt = std::max(fMaxTimeStep, d * d / (8. * D));
        } else if (fNearDNAMinStep > 0. && D > 0.) {
            // inside an envelope: refine where a coarse step could reach a solid
            const G4double d2 = DistanceSquaredToNearestBasePair(p);
            if (d2 < DBL_MAX) {
                const G4double gap = std::sqrt(d2) - fDNARadius;
                if (gap <= 0.) dt = fNearDNAMinStep;
                else dt = std::min(fMaxTimeStep, std::max(fNearDNAMinStep, gap * gap / (8. * D)));
            }
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
