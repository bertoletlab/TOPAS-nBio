// A plasmid built with the nucleus's own DNA geometry
// ********************************************************************
// *																  *
// * This file is part of the TOPAS-nBio extensions to the			  *
// *   TOPAS Simulation Toolkit.									  *
// * The TOPAS-nBio extensions are freely available under the license *
// *   agreement set forth at: https://topas-nbio.readthedocs.io/	  *
// *																  *
// ********************************************************************

#ifndef TsPlasmidSphereDNA_hh
#define TsPlasmidSphereDNA_hh

#include "TsVGeometryComponent.hh"
#include "G4ThreeVector.hh"

#include <vector>

class TsPlasmidSphereDNA : public TsVGeometryComponent
{
public:
	TsPlasmidSphereDNA(TsParameterManager* pM, TsExtensionManager* eM, TsMaterialManager* mM,
					   TsGeometryManager* gM, TsVGeometryComponent* parentComponent,
					   G4VPhysicalVolume* parentVolume, G4String& name);
	~TsPlasmidSphereDNA();

	G4VPhysicalVolume* Construct();

private:
	// The physical DNA model (TsPhysicalDNASolids): mass-faithful slabs, no overlaps,
	// continuous strands. Selected by Ge/<name>/DNAModel = "Physical". The path arrives
	// centred, closed and in Geant4 units.
	G4VPhysicalVolume* ConstructPhysical(const std::vector<G4ThreeVector>& path);

	G4int fNumberOfBasePairs;
};

#endif
