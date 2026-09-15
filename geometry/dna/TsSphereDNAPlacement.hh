// Shared placement for the sphere DNA volume model
// ********************************************************************
// *																  *
// * This file is part of the TOPAS-nBio extensions to the			  *
// *   TOPAS Simulation Toolkit.									  *
// * The TOPAS-nBio extensions are freely available under the license *
// *   agreement set forth at: https://topas-nbio.readthedocs.io/	  *
// *																  *
// ********************************************************************
//
// One implementation of where the sphere model's four DNA volumes sit around a path.
//
// It exists because there were two. TsNucleus and TsPlasmidSupercoiled both offered a DNA model
// called "Sphere" and built different geometries under that name, differing in the helix radius,
// both radial positions, the helix handedness, the plane normal, and whether the base was
// oriented within its base-pair plane at all. Measured at a backbone scavenging probability of
// 0.25 and a base probability of 1.0, the backbone share of radical attacks came out
// 0.1328 +/- 0.0058 in the plasmid against 0.2950 +/- 0.0096 in the nucleus. Switching histones
// and the hydration shell off in the nucleus moved nothing, which is what identified the
// placement rather than the environment as the cause.
//
// The consequence is that a scavenging probability calibrated in one geometry was not valid in
// the other, which matters because the experimental anchors for that parameter are naked-plasmid
// measurements while the model is applied to a cell nucleus. Sharing the placement makes the two
// identical by construction instead of by keeping constants in step by hand.

#ifndef TsSphereDNAPlacement_hh
#define TsSphereDNAPlacement_hh

#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include "G4SystemOfUnits.hh"

#include <vector>

// Where one base pair's volumes go. The two rotations differ: the base carries an extra
// rotateZ so a non-spherical base is oriented within its own base-pair plane, which matters
// because the sphere model's base is an ellipsoid, not a sphere, despite the model's name.
struct TsSphereDNABasePair
{
	G4int            bpID;
	G4ThreeVector    base1, base2, back1, back2;
	G4RotationMatrix rotBase;       // bases and, in the nucleus, the hydration shell
	G4RotationMatrix rotBackbone;   // backbones
};

// Default geometry of the sphere model, as built by TsNucleus.
namespace TsSphereDNADefaults
{
	const G4double kHelixRadius    = 1.2   * CLHEP::nanometer;
	const G4double kBackboneRadius = 0.271 * CLHEP::nanometer;
	// radial extent of the base, i.e. the ellipsoid semi-axis that points along the helix
	// radius. TsNucleus's base ellipsoid is 0.328 / 0.328 / 0.185 nm.
	const G4double kBaseRadialExtent = 0.328 * CLHEP::nanometer;
	const G4int    kBasePairsPerTurn = 10;
}

// Build the placement for every base pair along a path. The path is consumed pairwise, so the
// result has one entry fewer than the path has points, exactly as both original loops did.
// Resample an arbitrary path onto one point per base pair, 0.34 nm apart. The nucleus's
// Hilbert path already arrives at that pitch; a plasmid path read from file does not.
std::vector<G4ThreeVector> TsSegmentPathToBasePairs(
	const std::vector<G4ThreeVector>& path,
	G4double rise = 0.34 * CLHEP::nanometer);

std::vector<TsSphereDNABasePair> TsBuildSphereDNAPlacement(
	const std::vector<G4ThreeVector>& path,
	G4double helixRadius      = TsSphereDNADefaults::kHelixRadius,
	G4double backboneRadius   = TsSphereDNADefaults::kBackboneRadius,
	G4double baseRadialExtent = TsSphereDNADefaults::kBaseRadialExtent,
	G4int    basePairsPerTurn = TsSphereDNADefaults::kBasePairsPerTurn);

#endif
