// Physical DNA solids: base pairs as mass-faithful, non-overlapping, continuous slabs
// ********************************************************************
// *																  *
// * This file is part of the TOPAS-nBio extensions to the			  *
// *   TOPAS Simulation Toolkit.									  *
// * The TOPAS-nBio extensions are freely available under the license *
// *   agreement set forth at: https://topas-nbio.readthedocs.io/	  *
// *																  *
// ********************************************************************
//
// The sphere DNA model (TsSphereDNAPlacement) carries 0.31 of the base-pair mass in four small
// solids with gaps between them, so direct damage sees a third of the DNA and a diffusing
// radical can pass between the solids where a duplex would have intercepted it. Measured on
// the plasmid against Milligan 1993: a direct plateau at 0.35 of the measurement, an indirect
// yield that rose 1.7-fold when the solids were scaled up, and a backbone attack count that
// was not converged in the chemistry time step because a 0.27 nm orb is stepped over.
//
// This model builds each base pair as a slab of one rise along the path, partitioned by radius
// and angle into two half-cylinder bases, two annular backbone sectors on the strands and two
// hydration-shell sectors outside them, leaving the grooves open. Within a slab nothing can
// overlap; between slabs, the end faces are the planes perpendicular to the path legs at the
// leg midpoints, shared by consecutive slabs (G4CutTubs), so neighbours meet on one plane. The
// frame is parallel-transported along the path and, for a closed path, the closure mismatch
// is spread over the base pairs as a small extra twist, so the strands close without a seam.
//
// At the density of DNA (1.7 g/cm3) the base-pair mass fits inside a duplex of 1.0 nm radius
// with the grooves open; the density is a material property set in the deck, not here.
//
// Volume names, copy numbers and strand materials are those of the sphere model, so
// TsScoreDNADamageSBS, the SDD output and the attack counters need no change.

#ifndef TsPhysicalDNASolids_hh
#define TsPhysicalDNASolids_hh

#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include "G4SystemOfUnits.hh"

#include <vector>
#include <map>

class G4VSolid;

struct TsPhysicalDNAParams
{
	G4double rise             = 0.34  * CLHEP::nanometer;   // slab height per base pair
	G4double baseRadius       = 0.486 * CLHEP::nanometer;   // bases: 0 .. baseRadius
	G4double backboneOuter    = 1.00  * CLHEP::nanometer;   // backbone: baseRadius .. backboneOuter
	G4double shellOuter       = 1.47  * CLHEP::nanometer;   // hydration shell: backboneOuter .. shellOuter
	G4double backboneSpan     = 86.7  * CLHEP::degree;      // angular width of each backbone sector
	G4double strandSeparation = 180.0 * CLHEP::degree;      // strand 2 sector centre relative to strand 1
	G4double basePairsPerTurn = 10.5;
};

// One base pair's slab: where it sits, how it is oriented, and which solid set it uses.
struct TsPhysicalDNABasePair
{
	G4int            bpID;
	G4ThreeVector    center;      // slab centre (the path point, shifted along the axis when the two legs differ in length)
	G4RotationMatrix active;      // local -> world; columns are the world directions of local x, y, z
	G4double         twist;       // angle of strand 1 about the slab axis, already folded into `active`
	G4double         theta;       // bend angle between the incoming and outgoing legs
	G4double         phiLocal;    // azimuth of the bending direction in the slab frame
	G4double         dz;          // half-height of the cut tube
	G4double         shift;       // centre shift along the axis from the path point
	G4int            key;         // index into the solid-set table (one exact set per base pair)
};

struct TsPhysicalDNASolidSet
{
	G4VSolid* base1 = nullptr;
	G4VSolid* base2 = nullptr;
	G4VSolid* back1 = nullptr;
	G4VSolid* back2 = nullptr;
	G4VSolid* shell1 = nullptr;
	G4VSolid* shell2 = nullptr;
	G4double theta = 0, phiLocal = 0, dz = 0;   // the bin values the set was built for
};

// Smooth a closed polyline of control points into a curve and resample it at one point per
// rise, uniformly by arc length. A path file lists beads several nanometres apart, and the
// corner between two consecutive beads can be tens of degrees where the chain loops, which a
// slab of one rise cannot follow; a Catmull-Rom spline through the beads spreads that bend over
// the rise steps in between (a loop of 5 nm radius bends 3.9 deg per base pair). The curve is
// scaled about the origin so that its arc length equals the polyline's, which is the contour
// length the path file encodes, and the returned points number round(length / rise), with the
// last leg closing on the first. `scale` reports the factor applied. Where the smoothed curve
// bends more than `maxBend` per rise it is relaxed locally until it does not (0 disables).
std::vector<G4ThreeVector> TsSmoothResampleClosedPath(
	const std::vector<G4ThreeVector>& controlPoints, G4double rise, G4double& scale,
	G4double maxBend = 8.0 * CLHEP::degree);

// Build the placement for every base pair along a path already resampled to one point per
// rise (TsSmoothResampleClosedPath for a closed path, TsSegmentPathToBasePairs for an open
// one). For a closed path the last leg runs back to the first point; a closing gap shorter
// than half a rise is absorbed into the previous leg. One solid set is built per base pair,
// its end faces computed from that base pair's own legs, so that two slabs sharing a face
// describe it identically; a shared table of binned bend angles was tried first and left
// picometre crossings at the shared faces. Aborts (returns an empty vector and sets `error`)
// if any bend is too sharp for a cut tube of the shell radius, since then consecutive slabs
// would overlap on the inside of the bend.
std::vector<TsPhysicalDNABasePair> TsBuildPhysicalDNAPlacement(
	const std::vector<G4ThreeVector>& path,
	const TsPhysicalDNAParams& p,
	G4bool closed,
	std::vector<TsPhysicalDNASolidSet>& solidSets,   // filled: one entry per base pair
	G4String& error);

// Build the six solids for one base pair (theta, phiLocal, dz), named with the key.
TsPhysicalDNASolidSet TsBuildPhysicalDNASolidSet(
	const TsPhysicalDNAParams& p, G4double theta, G4double phiLocal, G4double dz, G4int key);

#endif
