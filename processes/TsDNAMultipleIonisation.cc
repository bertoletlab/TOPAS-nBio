// Extra Class for TsEmDNAChemistry
// Multiple ionisation of water by light ions (fork branch feat/multiple-ionisation). See the header.
#include "TsDNAMultipleIonisation.hh"

#include "G4DNAChemistryManager.hh"
#include "G4H2O.hh"
#include "G4Molecule.hh"
#include "G4MolecularConfiguration.hh"
#include "G4MoleculeTable.hh"
#include "G4ITTrackHolder.hh"
#include "G4Track.hh"
#include "G4VUserPulseInfo.hh"
#include "G4SystemOfUnits.hh"
#include "G4Exception.hh"
#include "Randomize.hh"
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>

std::atomic<long> TsDNAMultipleIonisation::nEvents[5] = {0, 0, 0, 0, 0};
std::atomic<long> TsDNAMultipleIonisation::ionEnergyLost_meV{0};

namespace {
    // He2+ alpha(x), x = (E/u)/2, recovered from the vector geometry of Fig. S2 of Okada et al. 2025 (Sci Rep, DOI 10.1038/s41598-025-00875-x), the curve of
    // Meesungnoen and Jay-Gerin 2005 (adjusted to reproduce measured G(HO2 + O2-)); no measurement of sigma_di / sigma_si in liquid water stands behind it.
    const G4double kX[10]     = {0.15, 0.2, 0.3, 0.4, 0.5, 0.7, 1.0, 1.5, 2.0, 2.5};
    const G4double kAlpha[10] = {0.1415, 0.1307, 0.1037, 0.0801, 0.0631, 0.0414, 0.0255, 0.0141, 0.0089, 0.0059};
    const G4double kXmin = 0.15, kXmax = 2.97;   // MJG05 plotted He2+ from x = 0.15 (0.3 MeV/u); the curve ends at 2.97
    const G4double kAlphaAtXmin = 0.1415 * 1.0;   // see Alpha(): clamped to the first and last tabulated point inside the span
    const G4double kAlphaAtXmax = 0.0059 * 1.0;
}

G4double TsDNAMultipleIonisation::Alpha(const G4String& particle, G4double kineticEnergy) const
{
    if ( particle != "alpha" && particle != "alpha+" && particle != "helium" ) return 0.;
    const G4double massU = 4.0015;                       // He-4 in u
    const G4double x = (kineticEnergy / MeV / massU) / 2.0;   // Z = 2 at every charge state (MJG05 neglected charge exchange)
    if ( x < kXmin || x > kXmax ) return 0.;
    G4double a;
    if ( x <= kX[0] ) a = kAlphaAtXmin;
    else if ( x >= kX[9] ) a = kAlphaAtXmax;
    else {
        int i = int(std::upper_bound(kX, kX + 10, x) - kX);   // kX[i-1] <= x < kX[i]
        const G4double t = (std::log(x) - std::log(kX[i-1])) / (std::log(kX[i]) - std::log(kX[i-1]));
        a = std::exp(std::log(kAlpha[i-1]) + t * (std::log(kAlpha[i]) - std::log(kAlpha[i-1])));
    }
    return fScale * a;
}

G4double TsDNAMultipleIonisation::CrossSectionFactor(const G4String& particle, G4double kineticEnergy) const
{
    if ( !fActive ) return 1.;
    const G4double a = Alpha(particle, kineticEnergy);
    if ( a <= 0. ) return 1.;
    G4double f = 1., w = 1.;
    for ( G4int n = 2; n <= fMaxOrder; ++n ) { w *= a; f += w; }
    return f;
}

G4int TsDNAMultipleIonisation::SampleOrder(const G4String& particle, G4double kineticEnergy) const
{
    if ( !fActive ) return 1;
    const G4double a = Alpha(particle, kineticEnergy);
    if ( a <= 0. ) return 1;
    G4double wsum = 1., w = 1.;
    for ( G4int n = 2; n <= fMaxOrder; ++n ) { w *= a; wsum += w; }
    const G4double u = G4UniformRand() * wsum;
    G4double acc = 1.; w = 1.;
    if ( u < acc ) return 1;
    for ( G4int n = 2; n <= fMaxOrder; ++n ) { w *= a; acc += w; if ( u < acc ) return n; }
    return fMaxOrder;
}

G4double TsDNAMultipleIonisation::Threshold(G4int order)
{
    return order >= 4 ? 88. * eV : (order == 3 ? 65. * eV : (order == 2 ? 40. * eV : 0.));
}

G4String TsDNAMultipleIonisation::StateName(G4int order)
{
    return "IonisationMulti_" + std::to_string(order);
}

void TsDNAMultipleIonisation::CreateMultiplyIonisedWater(const G4Track* pIncomingTrack, G4int order)
{
    if ( !G4DNAChemistryManager::Instance()->IsChemistryActivated() ) return;
    auto pH2OMolecule = new G4Molecule(G4H2O::Definition());
    const G4int orbits[4] = {4, 4, 3, 3};   // holes of the order-n state: two in the outer orbital, then two in the next
    for ( G4int k = 0; k < order && k < 4; ++k ) pH2OMolecule->IonizeMolecule(orbits[k]);
    if ( pH2OMolecule->GetMolecularConfiguration()->GetLabel() != StateName(order) ) {
        // the registered, named state carries the decay channels; an unnamed configuration would never dissociate
        delete pH2OMolecule;
        G4Exception("TsDNAMultipleIonisation::CreateMultiplyIonisedWater", "TsMI002", FatalException,
                    "A multiply ionised water state was requested but not registered: set b:Ch/<name>/IncludeMultipleIonisation = \"True\".");
        return;
    }

    G4double delayedTime = 0.;
    if ( pIncomingTrack->GetUserInformation() != nullptr ) {
        auto pPulseInfo = dynamic_cast<G4VUserPulseInfo*>(pIncomingTrack->GetUserInformation());
        if ( pPulseInfo != nullptr ) delayedTime = pPulseInfo->GetDelayedTime();
    }
    G4Track* pH2OTrack = pH2OMolecule->BuildTrack(picosecond + delayedTime, pIncomingTrack->GetPosition());
    pH2OTrack->SetParentID(pIncomingTrack->GetTrackID());
    pH2OTrack->SetTrackStatus(fStopButAlive);
    pH2OTrack->SetKineticEnergy(0.);
    G4ITTrackHolder::Instance()->Push(pH2OTrack);
}

static void PrintMultipleIonisationSummary()
{
    const long n1 = TsDNAMultipleIonisation::nEvents[1], n2 = TsDNAMultipleIonisation::nEvents[2], n3 = TsDNAMultipleIonisation::nEvents[3], n4 = TsDNAMultipleIonisation::nEvents[4];
    const double lostEv = TsDNAMultipleIonisation::ionEnergyLost_meV.load() / 1000.0;
    std::cout << "TsDNAMultipleIonisation summary: ionisation events single " << n1 << ", double " << n2 << ", triple " << n3 << ", quadruple " << n4
              << "; ion energy lost in them " << lostEv << " eV";
    if ( lostEv > 0. ) std::cout << "; per 100 eV of ion energy lost: single " << 100. * n1 / lostEv << ", double " << 100. * n2 / lostEv
                                 << ", triple " << 100. * n3 / lostEv << ", quadruple " << 100. * n4 / lostEv << "; double fraction of events " << (n1 + n2 + n3 + n4 > 0 ? double(n2) / (n1 + n2 + n3 + n4) : 0.)
                                 << ", triple fraction " << (n1 + n2 + n3 + n4 > 0 ? double(n3) / (n1 + n2 + n3 + n4) : 0.) << ", quadruple fraction " << (n1 + n2 + n3 + n4 > 0 ? double(n4) / (n1 + n2 + n3 + n4) : 0.);
    std::cout << std::endl;
}

void TsDNAMultipleIonisation::RegisterSummary()
{
    static bool registered = false;
    if ( registered ) return;
    registered = true;
    std::atexit(PrintMultipleIonisationSummary);
}
