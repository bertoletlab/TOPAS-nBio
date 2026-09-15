// Extra Class for use by TsNucleus and TsPlasmidSphereDNA
// ********************************************************************
// *																  *
// * This file is part of the TOPAS-nBio extensions to the			  *
// *   TOPAS Simulation Toolkit.									  *
// * The TOPAS-nBio extensions are freely available under the license *
// *   agreement set forth at: https://topas-nbio.readthedocs.io/	  *
// *																  *
// ********************************************************************
//
// The body below is TsNucleus::PlaceDNASphere's placement arithmetic, moved out unchanged so
// that any component wanting the sphere model gets that one geometry. The radial positions are
// derived from the solid sizes rather than written as literals, which is what keeps a resized
// sphere in the place its size implies; TsNucleus previously wrote 0.2715 and 0.3285 nm as
// constants beside a backbone radius of 0.271 and a base semi-axis of 0.328, so changing either
// solid left its placement behind.

#include "TsSphereDNAPlacement.hh"

#include <cmath>
#include <cstdlib>
#include <fstream>

namespace
{
	// Rotate a 3-vector by a rotation matrix. Kept identical to TsNucleus::ApplyRotation,
	// including its use of the transpose convention, so the placement is bit-for-bit the same.
	void ApplyRotation(G4ThreeVector& rotated, const G4ThreeVector& v,
					   const G4RotationMatrix& rot)
	{
		rotated[0] = v[0]*rot.xx() + v[1]*rot.yx() + v[2]*rot.zx();
		rotated[1] = v[0]*rot.xy() + v[1]*rot.yy() + v[2]*rot.zy();
		rotated[2] = v[0]*rot.xz() + v[1]*rot.yz() + v[2]*rot.zz();
	}
}

std::vector<G4ThreeVector> TsSegmentPathToBasePairs(
	const std::vector<G4ThreeVector>& path, G4double rise)
{
	// Same walk as DNACoordinates::Segment: step along each leg in whole rises, carrying the
	// remainder into the next leg rather than restarting at each control point.
	std::vector<G4ThreeVector> out;
	if (path.empty())
		return out;
	out.push_back(path[0]);
	G4int counter = 0;
	for (size_t i = 0; i + 1 < path.size(); i++)
	{
		G4ThreeVector v = path[i+1] - out[counter];
		G4int nDiv = (G4int)(v.mag() / rise);
		G4ThreeVector unit = v.unit() * rise;
		for (G4int j = 0; j < nDiv; j++)
		{
			out.push_back(out[counter] + unit);
			counter++;
		}
	}
	return out;
}

std::vector<TsSphereDNABasePair> TsBuildSphereDNAPlacement(
	const std::vector<G4ThreeVector>& path,
	G4double helixRadius,
	G4double backboneRadius,
	G4double baseRadialExtent,
	G4int    basePairsPerTurn)
{
	std::vector<TsSphereDNABasePair> out;
	if (path.size() < 2)
		return out;

	const G4double rotPair = (2.0 * CLHEP::pi) / (G4double)basePairsPerTurn;
	const G4double rBack   = helixRadius - backboneRadius;
	const G4double rBase   = rBack - backboneRadius - baseRadialExtent;

	const G4int nBP = (G4int)path.size();
	out.reserve(nBP - 1);

	for (G4int bp = 0; bp < nBP - 1; bp++)
	{
		const G4double angle1 = (G4double)bp * rotPair;
		const G4double angle2 = angle1 + CLHEP::pi;   // strand 2, across the helix

		G4ThreeVector back1temp(rBack*std::cos(angle1), rBack*std::sin(angle1), 0.0);
		G4ThreeVector back2temp(rBack*std::cos(angle2), rBack*std::sin(angle2), 0.0);
		G4ThreeVector base1temp(rBase*std::cos(angle1), rBase*std::sin(angle1), 0.0);
		G4ThreeVector base2temp(rBase*std::cos(angle2), rBase*std::sin(angle2), 0.0);

		// Rotation carrying this base pair's plane onto the next step of the path
		G4ThreeVector vecNext = (path[bp] - path[bp+1]).unit();
		G4ThreeVector norm(0., 0., -1.);            // Geant4 builds planes facing -z
		G4double AngBetween = std::acos(norm.dot(vecNext));
		G4ThreeVector cross = (vecNext.cross(norm)).unit();

		TsSphereDNABasePair e;
		e.bpID = bp + 1;

		const G4bool turning = (cross.x() != 0 || cross.y() != 0 || cross.z() != 0);
		if (turning)
		{
			G4RotationMatrix rot;
			rot.rotate(AngBetween, cross);
			ApplyRotation(e.back1, back1temp, rot);
			ApplyRotation(e.back2, back2temp, rot);
			ApplyRotation(e.base1, base1temp, rot);
			ApplyRotation(e.base2, base2temp, rot);
			e.rotBackbone.rotate(AngBetween, cross);
			e.rotBase.rotate(AngBetween, cross);
		}
		else
		{
			e.back1 = back1temp;
			e.back2 = back2temp;
			e.base1 = base1temp;
			e.base2 = base2temp;
		}

		// The base carries one further turn about its own axis so that a non-spherical base
		// lies correctly within the base-pair plane. The sphere model's base is an ellipsoid,
		// so omitting this silently misorients every base, which is what the plasmid did.
		e.rotBase.rotateZ(-angle1 + CLHEP::pi/2.);

		e.back1 += path[bp];
		e.back2 += path[bp];
		e.base1 += path[bp];
		e.base2 += path[bp];

		out.push_back(e);
	}

	// Debug facility: set TSNBIO_DUMP_DNA to a filename to append the placed positions, so the
	// geometry two components actually build can be compared directly rather than inferred
	// from yields. Capped, because a nucleus places far too many base pairs to write them all.
	if (const char* dumpPath = std::getenv("TSNBIO_DUMP_DNA"))
	{
		static G4int written = 0;
		const G4int cap = 200000;
		if (written < cap)
		{
			std::ofstream f(dumpPath, std::ios::app);
			for (size_t i = 0; i < out.size() && written < cap; i++, written++)
			{
				const TsSphereDNABasePair& e = out[i];
				f << e.base1.x()/CLHEP::nanometer << " " << e.base1.y()/CLHEP::nanometer << " "
				  << e.base1.z()/CLHEP::nanometer << " "
				  << e.back1.x()/CLHEP::nanometer << " " << e.back1.y()/CLHEP::nanometer << " "
				  << e.back1.z()/CLHEP::nanometer << "\n";
			}
		}
	}
	return out;
}
