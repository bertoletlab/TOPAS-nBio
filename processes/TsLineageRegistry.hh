// Radial classes: molecular tracks do not carry the ID of the physical track that made them (their ParentID is 0), so the origin of a molecule
// is classed by the distance of its birth point from the axis of the primary ion (set from the first step of the primary of the event): bins
// below 2, 2 to 5, 5 to 10, 10 to 30 nm and above 30 nm; class 5 holds molecules whose birth was not seen.
//
// Lineage classes of physical and molecular tracks, for the optional core/delta split of the G-value scorer (TsSBSScoreGValue).
//
// A physical track belongs to class 0 when it is the primary ion (or a charge-change continuation of it) and, when it descends from an
// electron made DIRECTLY by that ion, to the class of that electron's initial kinetic energy: 1 below T1, 2 from T1 to T2, 3 from T2 to T3,
// 4 above T3 (descendants inherit the class of the root electron). A molecular track inherits the class of its parent track. Class 5 holds
// tracks whose parent was not classified. Thread-local, since TOPAS workers run their own scorer and processes.
#ifndef TsLineageRegistry_hh
#define TsLineageRegistry_hh

#include "G4Track.hh"
#include "G4ParticleDefinition.hh"
#include "globals.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include <map>
#include <vector>

class TsLineageRegistry {
public:
    static const G4int kNClasses = 6;
    static TsLineageRegistry& Get() { static thread_local TsLineageRegistry r; return r; }

    G4bool active = false;
    G4double t1 = 100. * eV, t2 = 500. * eV, t3 = 1000. * eV;
    std::map<G4int, G4int> cls;                                   // trackID -> class
    G4bool axisSet = false; G4ThreeVector p0, dir;
    G4double edges[4] = {2. * nm, 5. * nm, 10. * nm, 30. * nm};
    std::map<G4int, G4int> radCls;                                // molecule trackID -> radial class of its birth point
    std::vector<std::pair<G4double, G4int> > captured;            // (time, class) of every product made by a first-order reaction

    void Clear() { cls.clear(); radCls.clear(); captured.clear(); axisSet = false; }
    void SetAxis(const G4ThreeVector& p, const G4ThreeVector& d) { if (!axisSet) { p0 = p; dir = d.unit(); axisSet = true; } }
    G4int RadialBin(const G4ThreeVector& p) const {
        if (!axisSet) return kNClasses - 1;
        G4ThreeVector v = p - p0; G4double r = (v - v.dot(dir) * dir).mag();
        for (int i = 0; i < 4; i++) if (r < edges[i]) return i;
        return 4;
    }
    G4int RadialClassOfMolecule(G4int id) const { auto it = radCls.find(id); return it == radCls.end() ? kNClasses - 1 : it->second; }
    void SetIon(G4int id) { cls[id] = 0; }

    // class of a track; creates and memoizes it from the parent on first sight (isNew tells whether it was created now)
    G4int ClassOf(const G4Track* track, G4bool* isNew = nullptr) {
        G4int id = track->GetTrackID();
        auto it = cls.find(id);
        if (it != cls.end()) { if (isNew) *isNew = false; return it->second; }
        if (isNew) *isNew = true;
        G4int parent = track->GetParentID(), c = kNClasses - 1;
        if (parent == 0) c = 0;
        else {
            auto pit = cls.find(parent);
            if (pit != cls.end()) {
                c = pit->second;
                if (id > 0 && c == 0 && track->GetDefinition()->GetParticleName() == "e-") {
                    G4double e = track->GetVertexKineticEnergy();
                    c = (e < t1) ? 1 : (e < t2) ? 2 : (e < t3) ? 3 : 4;
                }
            }
        }
        cls[id] = c;
        return c;
    }
};

#endif
