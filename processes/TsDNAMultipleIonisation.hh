// Multiple ionisation of water by light ions: the alpha(x) weights, the thresholds, the multiply ionised state and the creation of the H2O(n+) molecule
// (fork branch feat/multiple-ionisation, default off). Two switchable, internally consistent schemes:
//
//  MJG05  (Meesungnoen and Jay-Gerin 2005, J Phys Chem A 109:6406): sigma_si is kept, sigma_di, sigma_ti, sigma_qi = R, R^2, R^3 times sigma_si are ADDED on
//         top (the model cross section is multiplied by 1 + R + R^2 + R^3, and the LET then rises; MJG05 rescaled all ion cross sections to Watt's LET, which
//         is NOT done here: the ratio is reported and a rescaling needs PI approval). The energy transfer W of the event is drawn from the single-ionisation
//         spectrum under the condition W >= 40, 65, 88 eV (double, triple, quadruple); W minus that threshold is shared equally among the ejected electrons,
//         emitted 180, 120, 90 degrees apart; the threshold is deposited locally. Helium is He2+ at every energy (Z = 2 for He+ and He0 too).
//  GER06  (Gervais, Beuve, Olivera and Galassi 2006, Radiat Phys Chem 75:493): the total ionisation cross section is split among the degrees, so a q-fold event
//         replaces about q single events (event rate unchanged, weights 1 : R : R^2 : R^3 renormalised; GER06 tabulates no He, so R of He2+ stands in, which
//         agrees with their proton ratio sigma2/sigma1 within 15 percent at 0.25 to 1 MeV); each electron is drawn from the single-ionisation spectrum of its own
//         orbital (the K shell is excluded: Auger events stay single, as in Geant4).
//
// R = alpha(x) of He2+, x = (E / u) / 2, from the curve of Meesungnoen and Jay-Gerin 2005 Fig. 1 (values read by Puck from Okada et al. 2025 Fig. S2 and
// checked against the raster of the original within 3 percent); zero outside 0.3 to 6 MeV/u. Inputs: research/knowledge/water-radiolysis-yields/extractions/
// 2026-11-21_multiple-ionization-implementation-inputs.md. With the model inactive no random number is drawn and the cross section is unchanged.
#ifndef TsDNAMultipleIonisation_hh
#define TsDNAMultipleIonisation_hh 1

#include "globals.hh"
#include <atomic>
#include <vector>

class G4Track;

class TsDNAMultipleIonisation
{
public:
    enum Scheme { MJG05 = 0, Gervais = 1 };

    void Activate(G4bool on) { fActive = on; }
    G4bool IsActive() const { return fActive; }
    void SetScheme(Scheme s) { fScheme = s; }
    Scheme GetScheme() const { return fScheme; }
    // Multiplies alpha(x) (sensitivity knob, default 1).
    void SetScale(G4double scale) { fScale = scale; }
    // Highest order of ionisation per event, 2 to 4 (default 4: quadruple ionisation, MJG05 channel R4 and the O2+ route of GER06).
    void SetMaxOrder(G4int n) { fMaxOrder = (n < 2) ? 2 : (n > 4 ? 4 : n); }
    G4int GetMaxOrder() const { return fMaxOrder; }

    // alpha of He2+ (applied to alpha, alpha+ and helium alike, Z = 2) at this kinetic energy: log-log interpolation of the curve of MJG05 Fig. 1; 0 outside
    // x = 0.15 to 2.97 (0.3 to 5.9 MeV/u) and for any other projectile.
    G4double Alpha(const G4String& particle, G4double kineticEnergy) const;
    // Factor of the ionisation cross section (event rate): MJG05 1 + a + a^2 + a^3 (orders up to MaxOrder, a = Scale alpha, multiples added on top); Gervais the same sum
    // divided by 1 + 2a + 3a^2 + 4a^3 (a q-fold event replaces q singles); 1 where alpha is 0.
    G4double CrossSectionFactor(const G4String& particle, G4double kineticEnergy) const;
    // Order 1 to MaxOrder of an event, drawn with the weights 1 : a : a^2 : a^3; no random number when inactive or alpha is 0.
    G4int SampleOrder(const G4String& particle, G4double kineticEnergy) const;
    // Smallest energy transfer of an event of this order, 40, 65, 88 eV (MJG05 p. 6408).
    static G4double Threshold(G4int order);

    // Name of the multiply ionised water state of this order (IonisationMulti_<order>; holes in orbits 4, 4, 3, 3).
    static G4String StateName(G4int order);
    // Builds the H2O(n+) track at the place of the incoming track (as G4DNAChemistryManager::CreateWaterMolecule does for one hole). Does nothing when the
    // chemistry is inactive; raises a fatal exception when the state was not registered (b:Ch/<name>/IncludeMultipleIonisation).
    static void CreateMultiplyIonisedWater(const G4Track* incomingTrack, G4int order);

    // Event counters (this process): ionisation events of order 1 to 4, the ion energy lost in them (units of 1e-3 eV) and the extra events per 100 eV.
    static std::atomic<long> nEvents[5];
    static std::atomic<long> ionEnergyLost_meV;
    // Registers (once) a summary printed at process exit.
    static void RegisterSummary();

private:
    G4double fScale{1.0};
    G4int fMaxOrder{4};
    G4bool fActive{false};
    Scheme fScheme{MJG05};
};

#endif
