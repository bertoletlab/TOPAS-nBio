// Component for TsPlasmidSphereDNA
// ********************************************************************
// *																  *
// * This file is part of the TOPAS-nBio extensions to the			  *
// *   TOPAS Simulation Toolkit.									  *
// * The TOPAS-nBio extensions are freely available under the license *
// *   agreement set forth at: https://topas-nbio.readthedocs.io/	  *
// *																  *
// ********************************************************************
//
// A supercoiled plasmid whose DNA volumes are placed by TsBuildSphereDNAPlacement, the same
// routine TsNucleus uses. It exists so that a parameter calibrated against a naked-DNA
// experiment can be used in a cell nucleus.
//
// TsPlasmidSupercoiled also offers a DNA model called "Sphere", but it is a different geometry:
// helix radius 1.15 nm against 1.2, backbone radial position 0.860 nm against 0.9285, base
// radial position 0.270 nm against 0.3290, the opposite helix handedness, the opposite plane
// normal, and no rotation orienting the base within its own base-pair plane, which a
// non-spherical base needs. Measured at a backbone scavenging probability of 0.25 and a base
// probability of 1.0, the backbone share of radical attacks came out 0.1328 +/- 0.0058 there
// against 0.2950 +/- 0.0096 in the nucleus. Turning histones and the hydration shell off in the
// nucleus changed nothing, which is what identified the placement as the cause rather than the
// chromatin environment.
//
// The experimental anchors for the indirect chain are naked-plasmid measurements, while the
// model is applied to a nucleus, so that discrepancy made the calibration untransferable. This
// component removes it by construction rather than by keeping two sets of constants in step.

#include "TsPlasmidSphereDNA.hh"
#include "TsSphereDNAPlacement.hh"

#include "TsParameterManager.hh"

#include "G4Box.hh"
#include "G4Orb.hh"
#include "G4Ellipsoid.hh"
#include "G4LogicalVolume.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

#include <fstream>
#include <vector>

TsPlasmidSphereDNA::TsPlasmidSphereDNA(TsParameterManager* pM, TsExtensionManager* eM,
									   TsMaterialManager* mM, TsGeometryManager* gM,
									   TsVGeometryComponent* parentComponent,
									   G4VPhysicalVolume* parentVolume, G4String& name)
	: TsVGeometryComponent(pM, eM, mM, gM, parentComponent, parentVolume, name),
	  fNumberOfBasePairs(0)
{
}

TsPlasmidSphereDNA::~TsPlasmidSphereDNA()
{
}

G4VPhysicalVolume* TsPlasmidSphereDNA::Construct()
{
	BeginConstruction();

	// ---- the plasmid path ------------------------------------------------------------
	G4String fileName = fPm->GetStringParameter(GetFullParmName("FileName"));
	std::ifstream infile(fileName);
	if (!infile)
	{
		G4cerr << "Cannot open plasmid path file: " << fileName << G4endl;
		fPm->AbortSession(true);
	}

	std::vector<G4ThreeVector> path;
	G4double x, y, z;
	G4double xmin = 1.0*mm, xmax = -1.0*mm;
	G4double ymin = 1.0*mm, ymax = -1.0*mm;
	G4double zmin = 1.0*mm, zmax = -1.0*mm;
	while (true)
	{
		infile >> x >> y >> z;
		if (!infile.good())
			break;
		x *= nm; y *= nm; z *= nm;
		path.push_back(G4ThreeVector(x, y, z));
		if (x < xmin) xmin = x;
		if (y < ymin) ymin = y;
		if (z < zmin) zmin = z;
		if (x > xmax) xmax = x;
		if (y > ymax) ymax = y;
		if (z > zmax) zmax = z;
	}
	if (path.size() < 2)
	{
		G4cerr << "Plasmid path file holds fewer than two points: " << fileName << G4endl;
		fPm->AbortSession(true);
	}
	path.push_back(path[0]);   // close the loop

	G4ThreeVector offset(0.5*(xmax+xmin), 0.5*(ymax+ymin), 0.5*(zmax+zmin));
	for (size_t u = 0; u < path.size(); u++)
		path[u] -= offset;

	G4double HLX = 0.5*(xmax-xmin), HLY = 0.5*(ymax-ymin), HLZ = 0.5*(zmax-zmin);
	G4Box* envelope = new G4Box(fName, HLX + 0.5*3.4*nm, HLY + 0.5*3.4*nm, HLZ + 0.5*3.4*nm);
	fEnvelopeLog  = CreateLogicalVolume(envelope);
	fEnvelopePhys = CreatePhysicalVolume(fEnvelopeLog);

	// ---- the DNA solids, as TsNucleus builds them ------------------------------------
	G4double backboneRadius = TsSphereDNADefaults::kBackboneRadius;
	G4double baseSemiX = 0.328*nm, baseSemiY = 0.328*nm, baseSemiZ = 0.185*nm;
	if (fPm->ParameterExists(GetFullParmName("SphereBackboneRadius")))
		backboneRadius = fPm->GetDoubleParameter(GetFullParmName("SphereBackboneRadius"),"Length");
	if (fPm->ParameterExists(GetFullParmName("SphereBaseSemiAxisX")))
		baseSemiX = fPm->GetDoubleParameter(GetFullParmName("SphereBaseSemiAxisX"),"Length");
	if (fPm->ParameterExists(GetFullParmName("SphereBaseSemiAxisY")))
		baseSemiY = fPm->GetDoubleParameter(GetFullParmName("SphereBaseSemiAxisY"),"Length");
	if (fPm->ParameterExists(GetFullParmName("SphereBaseSemiAxisZ")))
		baseSemiZ = fPm->GetDoubleParameter(GetFullParmName("SphereBaseSemiAxisZ"),"Length");
	G4double helixRadius = TsSphereDNADefaults::kHelixRadius;
	if (fPm->ParameterExists(GetFullParmName("HelixRadius")))
		helixRadius = fPm->GetDoubleParameter(GetFullParmName("HelixRadius"),"Length");

	G4Ellipsoid* gBase     = new G4Ellipsoid("DNA_base", baseSemiX, baseSemiY, baseSemiZ);
	G4Orb*       gBackbone = new G4Orb("DNA_backbone", backboneRadius);

	G4LogicalVolume* lBase1 = CreateLogicalVolume("Base1", gBase);
	G4LogicalVolume* lBase2 = CreateLogicalVolume("Base2", gBase);
	G4LogicalVolume* lBack1 = CreateLogicalVolume("Backbone1", gBackbone);
	G4LogicalVolume* lBack2 = CreateLogicalVolume("Backbone2", gBackbone);

	// ---- place them, using the nucleus's placement -----------------------------------
	std::vector<G4ThreeVector> steps = TsSegmentPathToBasePairs(path);
	std::vector<TsSphereDNABasePair> bps =
		TsBuildSphereDNAPlacement(steps, helixRadius, backboneRadius, baseSemiX);

	for (size_t i = 0; i < bps.size(); i++)
	{
		TsSphereDNABasePair& e = bps[i];
		fNumberOfBasePairs++;

		G4ThreeVector* pBase1 = new G4ThreeVector(e.base1);
		G4ThreeVector* pBase2 = new G4ThreeVector(e.base2);
		G4ThreeVector* pBack1 = new G4ThreeVector(e.back1);
		G4ThreeVector* pBack2 = new G4ThreeVector(e.back2);
		G4RotationMatrix* rBase = new G4RotationMatrix(e.rotBase);
		G4RotationMatrix* rBack = new G4RotationMatrix(e.rotBackbone);

		CreatePhysicalVolume("Base1_", e.bpID, true, lBase1, rBase, pBase1, fEnvelopePhys);
		CreatePhysicalVolume("Base2_", e.bpID, true, lBase2, rBase, pBase2, fEnvelopePhys);
		CreatePhysicalVolume("Backbone1_", e.bpID, true, lBack1, rBack, pBack1, fEnvelopePhys);
		CreatePhysicalVolume("Backbone2_", e.bpID, true, lBack2, rBack, pBack2, fEnvelopePhys);
	}

	G4cout << G4endl << "Plasmid built with the nucleus sphere DNA model: "
		   << fNumberOfBasePairs << " bp" << G4endl << G4endl;

	InstantiateChildren(fEnvelopePhys);
	return fEnvelopePhys;
}
