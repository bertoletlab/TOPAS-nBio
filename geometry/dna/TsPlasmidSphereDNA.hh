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

class TsPlasmidSphereDNA : public TsVGeometryComponent
{
public:
	TsPlasmidSphereDNA(TsParameterManager* pM, TsExtensionManager* eM, TsMaterialManager* mM,
					   TsGeometryManager* gM, TsVGeometryComponent* parentComponent,
					   G4VPhysicalVolume* parentVolume, G4String& name);
	~TsPlasmidSphereDNA();

	G4VPhysicalVolume* Construct();

private:
	G4int fNumberOfBasePairs;
};

#endif
