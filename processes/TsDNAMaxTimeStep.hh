// Extra Class for TsEmDNAChemistry
//
// A cap on the chemistry time step, attached to every diffusing molecule.
//
// Why it exists. The step-by-step scheduler takes, for each molecule, the time to its nearest
// reactive partner, and nothing bounds that time from above: G4Scheduler::SetMaxTimeStep is
// assigned at the start of every step and overwritten before use, and the user time-step table
// is a floor, not a ceiling. Once a track has thinned out, a radical with no partner nearby is
// moved by one Gaussian jump covering the whole remaining interval, and Brownian transport
// checks for volume boundaries only along the straight segment of that jump. A hydroxyl radical
// jumping 100 nm in one step therefore crosses a 0.3 nm DNA sphere only if the segment happens to
// pass through it, whereas the real path over the same interval sweeps a volume of order
// 4 pi D R t. Encounters with DNA volumes are undercounted by roughly R / sqrt(D dt). Reactions
// between molecules are unaffected, since the encounter stepper brings partners into contact by
// construction; only geometry targets are lost. Measured 2026-09-29 on an 8-plasmid lattice at
// 30 us: hydroxyl attacks on DNA 48, 86, 118, 178 for no cap, 10 ns, 1 ns, 0.1 ns.
//
// What it does. Proposes a time at every step and does nothing when selected, so it never
// lengthens a step another process or a reaction would have made shorter. With only
// Ch/<Name>/MaximumTimeStep given, the cap is the same everywhere. With
// Ch/<Name>/MaximumTimeStepVolumePrefix also given, the cap applies inside every physical volume
// whose name starts with that prefix, grown by Ch/<Name>/MaximumTimeStepMargin; outside, the
// step is the time for the molecule to diffuse to the nearest such volume with 95% confidence,
// d^2 / (8 D), the rule G4DNABrownianTransportation uses for its own boundary estimate, but
// never less than the cap. DNA lives inside those volumes, so the walk is fine where it has to
// be and coarse where nothing can be hit. Off unless MaximumTimeStep is given.

#ifndef TsDNAMaxTimeStep_HH
#define TsDNAMaxTimeStep_HH

#include "G4VITProcess.hh"
#include "G4ThreeVector.hh"
#include <vector>
#include <unordered_map>

class TsDNAMaxTimeStep : public G4VITProcess
{
public:
    TsDNAMaxTimeStep(const G4String& aName = "DNAMaxTimeStep", G4ProcessType type = fDecay);
    virtual ~TsDNAMaxTimeStep();

    G4IT_ADD_CLONE(G4VITProcess, TsDNAMaxTimeStep)

    TsDNAMaxTimeStep(const TsDNAMaxTimeStep&);
    TsDNAMaxTimeStep& operator=(const TsDNAMaxTimeStep&);

    void StartTracking(G4Track*);

    void SetMaximumTimeStep(G4double maxTimeStep) { fMaxTimeStep = maxTimeStep; }
    // Empty prefix (the default) means the cap applies everywhere.
    void SetVolumePrefix(const G4String& prefix, G4double margin) { fVolumePrefix = prefix; fMargin = margin; }
    // Refinement near the DNA solids. Within the named envelopes the base-pair positions are
    // collected (the daughter named Base1_ of each base pair, in world coordinates) into a
    // grid, and a molecule whose coarse step could reach a base pair steps by
    // max(minStep, (d - dnaRadius)^2 / (8 D)) instead of the cap, so the walk is fine only
    // where a coarse step could have crossed a solid unseen. Measured on the sphere model at a
    // 0.1 ns cap, hydroxyl attacks at 1e9 /s rose 1.66-fold when the cap was lowered to 0.01 ns
    // everywhere, at seven times the wall time; this spends the fine steps where they count.
    // Off unless minStep > 0.
    void SetNearDNA(G4double minStep, G4double dnaRadius) { fNearDNAMinStep = minStep; fDNARadius = dnaRadius; }

    virtual void BuildPhysicsTable(const G4ParticleDefinition&);

    virtual G4double PostStepGetPhysicalInteractionLength(const G4Track& track,
                                                          G4double previousStepSize,
                                                          G4ForceCondition* condition);

    virtual G4VParticleChange* PostStepDoIt(const G4Track&, const G4Step&);

    virtual G4double AtRestGetPhysicalInteractionLength(const G4Track&, G4ForceCondition*) { return -1.0; }
    virtual G4VParticleChange* AtRestDoIt(const G4Track&, const G4Step&) { return 0; }
    virtual G4double AlongStepGetPhysicalInteractionLength(const G4Track&, G4double, G4double,
                                                           G4double&, G4GPILSelection*) { return -1.0; }
    virtual G4VParticleChange* AlongStepDoIt(const G4Track&, const G4Step&) { return 0; }

protected:
    struct MaxTimeStepState : public G4ProcessState
    {
        MaxTimeStepState() : G4ProcessState() {}
        virtual ~MaxTimeStepState() {}
    };

    // One axis-aligned box per matching physical volume, in world coordinates. A rotated
    // volume is covered by its bounding sphere instead (fIsSphere), which only ever enlarges
    // the fine-step region.
    struct Region { G4ThreeVector centre; G4ThreeVector halfLengths; G4double radius; G4bool isSphere; };

    void CollectRegions();
    G4double DistanceToNearestRegion(const G4ThreeVector& p) const;
    // squared distance to the nearest collected base-pair position, or DBL_MAX when none lies
    // within the grid's reach (one cell in every direction)
    G4double DistanceSquaredToNearestBasePair(const G4ThreeVector& p) const;

    G4ParticleChange fParticleChange;

private:
    void Create();
    G4double fMaxTimeStep;
    G4String fVolumePrefix;
    G4double fMargin;
    G4bool fRegionsCollected;
    std::vector<Region> fRegions;

    G4double fNearDNAMinStep;     // 0 = off
    G4double fDNARadius;
    G4double fCell;               // grid cell size; the search reaches one cell each way
    std::vector<G4ThreeVector> fBasePairs;
    std::unordered_map<long long, std::vector<int>> fGrid;
    long long CellKey(const G4ThreeVector& p, int dx = 0, int dy = 0, int dz = 0) const;
};

#endif
