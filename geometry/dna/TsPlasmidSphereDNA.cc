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
#include "TsPhysicalDNASolids.hh"

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

	// ---- which DNA model: the sphere solids (default) or the physical slabs --------------
	G4String dnaModel = "Sphere";
	if (fPm->ParameterExists(GetFullParmName("DNAModel")))
		dnaModel = fPm->GetStringParameter(GetFullParmName("DNAModel"));
	if (dnaModel == "Physical")
		return ConstructPhysical(path);
	if (dnaModel != "Sphere")
	{
		G4cerr << GetFullParmName("DNAModel") << " must be \"Sphere\" or \"Physical\", not \""
			   << dnaModel << "\"" << G4endl;
		fPm->AbortSession(true);
	}

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

	// ---- the hydration shell, as TsNucleus builds it, off unless asked for -----------
	// Quasi-direct damage is charge transfer from an ionised shell to the backbone, and the
	// scorer can only see it where a volume named HydrationShell1 or HydrationShell2 exists;
	// without the shell that channel is silently zero on a plasmid. The same parameter names as
	// the nucleus, and the same wedges, from the shared builder. The nucleus defaults the shell
	// on; this defaults it off so that every plasmid deck written before it existed builds the
	// geometry it always built.
	G4bool addHydrationShell = false;
	if (fPm->ParameterExists(GetFullParmName("AddHydrationShell")))
		addHydrationShell = fPm->GetBooleanParameter(GetFullParmName("AddHydrationShell"));
	G4LogicalVolume* lShell1 = NULL;
	G4LogicalVolume* lShell2 = NULL;
	if (addHydrationShell)
	{
		G4double shellThickness = TsSphereDNADefaults::kHydrationShellThickness;
		if (fPm->ParameterExists(GetFullParmName("HydrationShellThickness")))
			shellThickness = fPm->GetDoubleParameter(GetFullParmName("HydrationShellThickness"),"Length");
		G4double shellPhiSpan = TsSphereDNADefaults::kHydrationShellPhiSpan;
		if (fPm->ParameterExists(GetFullParmName("HydrationShellPhiSpan")))
			shellPhiSpan = fPm->GetDoubleParameter(GetFullParmName("HydrationShellPhiSpan"),"Angle");
		if (shellPhiSpan <= 0. || shellPhiSpan > 180*deg)
		{
			G4cerr << "TOPAS is exiting due to a serious error in Geometry setup." << G4endl;
			G4cerr << GetFullParmName("HydrationShellPhiSpan")
				   << " must lie in (0, 180] deg, or the two wedges overlap each other."
				   << G4endl;
			fPm->AbortSession(true);
		}
		TsSphereDNAHydrationShellSolids gShell =
			TsBuildSphereDNAHydrationShellSolids(backboneRadius, shellThickness, shellPhiSpan);
		lShell1 = CreateLogicalVolume("HydrationShell1", gShell.shell1);
		lShell2 = CreateLogicalVolume("HydrationShell2", gShell.shell2);
	}

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

		if (addHydrationShell)
		{
			// on the backbone positions with the base's orientation, as the nucleus places it
			G4ThreeVector* pShell1 = new G4ThreeVector(e.back1);
			G4ThreeVector* pShell2 = new G4ThreeVector(e.back2);
			G4RotationMatrix* rShell = new G4RotationMatrix(e.rotBase);
			CreatePhysicalVolume("HydrationShell1_", e.bpID, true, lShell1, rShell, pShell1, fEnvelopePhys);
			CreatePhysicalVolume("HydrationShell2_", e.bpID, true, lShell2, rShell, pShell2, fEnvelopePhys);
		}
	}

	G4cout << G4endl << "Plasmid built with the nucleus sphere DNA model: "
		   << fNumberOfBasePairs << " bp"
		   << (addHydrationShell ? ", with hydration shell" : ", no hydration shell")
		   << G4endl << G4endl;

	InstantiateChildren(fEnvelopePhys);
	return fEnvelopePhys;
}

G4VPhysicalVolume* TsPlasmidSphereDNA::ConstructPhysical(const std::vector<G4ThreeVector>& path)
{
	// ---- parameters, each defaulting to the design values of TsPhysicalDNAParams ----------
	TsPhysicalDNAParams p;
	auto length = [&](const char* name, G4double& v) {
		if (fPm->ParameterExists(GetFullParmName(name))) v = fPm->GetDoubleParameter(GetFullParmName(name), "Length"); };
	auto angle = [&](const char* name, G4double& v) {
		if (fPm->ParameterExists(GetFullParmName(name))) v = fPm->GetDoubleParameter(GetFullParmName(name), "Angle"); };
	length("PhysicalBaseRadius", p.baseRadius);
	length("PhysicalBackboneOuterRadius", p.backboneOuter);
	length("PhysicalShellOuterRadius", p.shellOuter);
	angle("PhysicalBackboneSpan", p.backboneSpan);
	angle("PhysicalStrandSeparation", p.strandSeparation);
	if (fPm->ParameterExists(GetFullParmName("PhysicalBasePairsPerTurn")))
		p.basePairsPerTurn = fPm->GetUnitlessParameter(GetFullParmName("PhysicalBasePairsPerTurn"));
	if (p.baseRadius <= 0 || p.backboneOuter <= p.baseRadius || p.shellOuter <= p.backboneOuter
		|| p.backboneSpan <= 0 || p.backboneSpan > p.strandSeparation
		|| p.strandSeparation + p.backboneSpan > CLHEP::twopi)
	{
		G4cerr << GetFullParmName("Physical*") << ": radii must increase base < backbone < shell and the "
			   << "two backbone sectors must not meet (span <= separation, separation + span <= 360 deg)" << G4endl;
		fPm->AbortSession(true);
	}

	G4bool addHydrationShell = false;
	if (fPm->ParameterExists(GetFullParmName("AddHydrationShell")))
		addHydrationShell = fPm->GetBooleanParameter(GetFullParmName("AddHydrationShell"));
	// The envelope margin (1.7 nm a side) must hold the outermost solid.
	const G4double outermost = addHydrationShell ? p.shellOuter : p.backboneOuter;
	if (outermost > 1.7 * nm)
	{
		G4cerr << GetFullParmName("PhysicalShellOuterRadius") << " exceeds the 1.7 nm envelope margin" << G4endl;
		fPm->AbortSession(true);
	}

	// ---- placement: smooth the bead path, one point per rise, then the slabs ------------
	G4double pathScale = 1.0;
	G4double maxBend = 8.0 * deg;
	angle("PhysicalMaxBendPerBasePair", maxBend);
	std::vector<G4ThreeVector> steps = TsSmoothResampleClosedPath(path, p.rise, pathScale, maxBend);
	if (steps.size() < 3)
	{
		G4cerr << "TsPlasmidSphereDNA (Physical): the path needs at least three distinct points" << G4endl;
		fPm->AbortSession(true);
	}
	std::vector<TsPhysicalDNASolidSet> sets;
	G4String error;
	std::vector<TsPhysicalDNABasePair> bps = TsBuildPhysicalDNAPlacement(steps, p, true, sets, error);
	if (!error.empty())
	{
		G4cerr << "TsPlasmidSphereDNA (Physical): " << error << G4endl;
		fPm->AbortSession(true);
	}

	// ---- logical volumes: one per solid set and moiety, materials by the usual names -----
	G4String matBase1 = fPm->GetStringParameter(GetFullParmName("Base1/Material"));
	G4String matBase2 = fPm->GetStringParameter(GetFullParmName("Base2/Material"));
	G4String matBack1 = fPm->GetStringParameter(GetFullParmName("Backbone1/Material"));
	G4String matBack2 = fPm->GetStringParameter(GetFullParmName("Backbone2/Material"));
	G4String matShell1, matShell2;
	if (addHydrationShell)
	{
		matShell1 = fPm->GetStringParameter(GetFullParmName("HydrationShell1/Material"));
		matShell2 = fPm->GetStringParameter(GetFullParmName("HydrationShell2/Material"));
	}
	std::vector<G4LogicalVolume*> lBase1, lBase2, lBack1, lBack2, lShell1, lShell2;
	for (size_t k = 0; k < sets.size(); k++)
	{
		lBase1.push_back(CreateLogicalVolume("Base1", matBase1, sets[k].base1));
		lBase2.push_back(CreateLogicalVolume("Base2", matBase2, sets[k].base2));
		lBack1.push_back(CreateLogicalVolume("Backbone1", matBack1, sets[k].back1));
		lBack2.push_back(CreateLogicalVolume("Backbone2", matBack2, sets[k].back2));
		if (addHydrationShell)
		{
			lShell1.push_back(CreateLogicalVolume("HydrationShell1", matShell1, sets[k].shell1));
			lShell2.push_back(CreateLogicalVolume("HydrationShell2", matShell2, sets[k].shell2));
		}
	}

	// ---- placement. G4PVPlacement takes the rotation of the mother frame seen from the
	// daughter, the inverse of the active rotation that carries the solid into the world. ----
	for (size_t i = 0; i < bps.size(); i++)
	{
		const TsPhysicalDNABasePair& e = bps[i];
		fNumberOfBasePairs++;
		G4RotationMatrix* rot = new G4RotationMatrix(e.active.inverse());
		G4ThreeVector* pos = new G4ThreeVector(e.center);
		CreatePhysicalVolume("Base1_", e.bpID, true, lBase1[e.key], rot, pos, fEnvelopePhys);
		CreatePhysicalVolume("Base2_", e.bpID, true, lBase2[e.key], rot, pos, fEnvelopePhys);
		CreatePhysicalVolume("Backbone1_", e.bpID, true, lBack1[e.key], rot, pos, fEnvelopePhys);
		CreatePhysicalVolume("Backbone2_", e.bpID, true, lBack2[e.key], rot, pos, fEnvelopePhys);
		if (addHydrationShell)
		{
			CreatePhysicalVolume("HydrationShell1_", e.bpID, true, lShell1[e.key], rot, pos, fEnvelopePhys);
			CreatePhysicalVolume("HydrationShell2_", e.bpID, true, lShell2[e.key], rot, pos, fEnvelopePhys);
		}
	}

	G4double maxTheta = 0;
	for (size_t i = 0; i < bps.size(); i++) maxTheta = std::max(maxTheta, bps[i].theta);
	G4cout << G4endl << "Plasmid built with the physical DNA solids: " << fNumberOfBasePairs << " bp, "
		   << "bead path smoothed and scaled by " << pathScale << " to its contour length, "
		   << sets.size() << " solid sets (one per base pair), largest bend " << maxTheta / deg << " deg, "
		   << "base radius " << p.baseRadius / nm << " nm, backbone to " << p.backboneOuter / nm
		   << " nm over " << p.backboneSpan / deg << " deg per strand, "
		   << p.basePairsPerTurn << " bp per turn"
		   << (addHydrationShell ? ", hydration shell to " : ", no hydration shell")
		   << (addHydrationShell ? std::to_string(p.shellOuter / nm) + " nm" : "")
		   << G4endl << G4endl;

	InstantiateChildren(fEnvelopePhys);
	return fEnvelopePhys;
}
