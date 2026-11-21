// Multiple ionisation of water by light ions: the alpha(x) weights, the thresholds, the multiply ionised state names and the creation of the
// H2O(n+) molecule (fork branch feat/multiple-ionisation, default off).
//
// Scheme of Meesungnoen and Jay-Gerin 2005 as implemented by Baba 2022 (thesis, Hokkaido) and Okada et al. 2025 (MPEXS2.1-DNA), read by Puck
// (research/knowledge/water-radiolysis-yields/extractions/2026-11-21_multiple-ionization-implementation-inputs.md): per ionisation of a He2+ ion, the
// relative weights single : double : triple are 1 : alpha : alpha^2, with x = (E / u) / Z. Still unsourced and therefore PLACEHOLDERS: the energy of
// the second and third electron (here the single-ionisation sampler at the projectile energy), the shells of the extra vacancies (here drawn from the
// single-ionisation shell cross sections, K shell excluded), whether the single-ionisation cross section is reduced (here not), the displacement
// scales, and every alpha outside the plotted range of the He2+ curve (0.27 to 5.9 MeV/u, off outside) and for He+ and He0 (off). With the model
// inactive no random number is drawn, so a deck that does not set the parameters is byte-identical.
#ifndef TsDNAMultipleIonisation_hh
#define TsDNAMultipleIonisation_hh 1

#include "globals.hh"
#include <atomic>
#include <vector>

class G4Track;

class TsDNAMultipleIonisation
{
public:
    void Activate(G4bool on) { fActive = on; }
    G4bool IsActive() const { return fActive; }
    // Multiplies alpha(x) (sensitivity knob, default 1).
    void SetScale(G4double scale) { fScale = scale; }
    // Highest order of ionisation per event: 2 (double) or 3 (triple). Quadruple ionisation (Okada only; transcription of its products doubtful) is not coded.
    void SetMaxOrder(G4int n) { fMaxOrder = (n < 2) ? 2 : (n > 3 ? 3 : n); }

    // alpha of a He2+ ion ("alpha") with kinetic energy kineticEnergy, from the curve of Okada et al. 2025 Fig. S2 as recovered by Puck (log-log
    // interpolation in x = E/u / 2); 0 for any other projectile and outside x = 0.135 to 2.97 (0.27 to 5.9 MeV/u).
    G4double Alpha(const G4String& particle, G4double kineticEnergy) const;
    // Order of this ionisation event: 1 (single), 2 or 3, drawn with the weights 1 : a : a^2 (a = Scale alpha, orders above MaxOrder dropped); no random
    // number is drawn when the model is inactive or alpha is 0.
    G4int SampleOrder(const G4String& particle, G4double kineticEnergy) const;
    // Smallest total energy transfer (binding plus ejected energies) of an event of this order: 40 eV double, 65 eV triple (Okada et al. 2025 Table S4).
    static G4double Threshold(G4int order);

    // Name of the multiply ionised water state with holes in the given shells (shell 0 = 1b1 ... 4 = K; orbit = 4 - shell), independent of their order.
    static G4String StateName(std::vector<G4int> shells);
    // Builds the H2O(n+) track at the place of the incoming track, as G4DNAChemistryManager::CreateWaterMolecule does for one hole. Does nothing when the
    // chemistry is inactive; raises a fatal exception when the state was not registered (b:Ch/<name>/IncludeMultipleIonisation).
    static void CreateMultiplyIonisedWater(const G4Track* incomingTrack, const std::vector<G4int>& shells);

    // Event counters (this process), printed by the model at the end of the run: ionisation events of order 1, 2, 3 and the ion energy lost in them.
    static std::atomic<long> nEvents[4];
    static std::atomic<long> ionEnergyLost_meV;   // in units of 1e-3 eV
    // Registers (once) a summary printed at process exit: events by order, ion energy lost, and the multiple ionisations per 100 eV of ion energy lost.
    static void RegisterSummary();

private:
    G4double fScale{1.0};
    G4int fMaxOrder{3};
    G4bool fActive{false};
};

#endif
