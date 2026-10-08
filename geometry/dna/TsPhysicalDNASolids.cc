// Extra Class for use by TsPlasmidSphereDNA and TsNucleus
// ********************************************************************
// *																  *
// * This file is part of the TOPAS-nBio extensions to the			  *
// *   TOPAS Simulation Toolkit.									  *
// * The TOPAS-nBio extensions are freely available under the license *
// *   agreement set forth at: https://topas-nbio.readthedocs.io/	  *
// *																  *
// ********************************************************************

#include "TsPhysicalDNASolids.hh"

#include "G4CutTubs.hh"
#include "G4Tubs.hh"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <tuple>

namespace
{
	// Rotation taking unit vector a onto unit vector b by the smallest angle; identity when
	// they coincide. Used to parallel-transport the reference vector of the slab frame.
	G4RotationMatrix MinimalRotation(const G4ThreeVector& a, const G4ThreeVector& b)
	{
		G4RotationMatrix r;
		G4ThreeVector axis = a.cross(b);
		G4double s = axis.mag();
		G4double c = a.dot(b);
		if (s < 1e-12)
			return r;              // parallel (anti-parallel cannot occur between consecutive legs of a sampled path)
		r.rotate(std::atan2(s, c), axis.unit());
		return r;
	}

}

namespace
{
	// n points uniformly spaced by arc length along a closed polyline
	std::vector<G4ThreeVector> ResampleClosedUniform(const std::vector<G4ThreeVector>& curve, G4int n)
	{
		std::vector<G4ThreeVector> out;
		const G4int nd = (G4int)curve.size();
		G4double length = 0;
		for (G4int i = 0; i < nd; i++)
			length += (curve[(i + 1) % nd] - curve[i]).mag();
		const G4double step = length / n;
		out.reserve(n);
		G4double walked = 0;
		G4int seg = 0;
		G4ThreeVector a = curve[0], b = curve[1 % nd];
		G4double segLen = (b - a).mag();
		for (G4int i = 0; i < n; i++)
		{
			const G4double target = i * step;
			while (walked + segLen < target && seg < nd)
			{
				walked += segLen;
				seg++;
				a = curve[seg % nd]; b = curve[(seg + 1) % nd];
				segLen = (b - a).mag();
			}
			const G4double f = segLen > 0 ? (target - walked) / segLen : 0;
			out.push_back(a + f * (b - a));
		}
		return out;
	}
}

std::vector<G4ThreeVector> TsSmoothResampleClosedPath(
	const std::vector<G4ThreeVector>& controlIn, G4double rise, G4double& scale, G4double maxBend)
{
	std::vector<G4ThreeVector> ctrl = controlIn;
	scale = 1.0;
	// a path file may repeat its first point at the end; the spline closes itself
	while (ctrl.size() > 2 && (ctrl.back() - ctrl.front()).mag() < 1e-6 * CLHEP::nanometer)
		ctrl.pop_back();
	const G4int m = (G4int)ctrl.size();
	std::vector<G4ThreeVector> out;
	if (m < 3)
		return out;

	// contour length of the polyline, including the closing leg
	G4double polyLength = 0;
	for (G4int i = 0; i < m; i++)
		polyLength += (ctrl[(i + 1) % m] - ctrl[i]).mag();

	// dense sampling of the uniform Catmull-Rom spline through the beads
	const G4double fine = 0.02 * CLHEP::nanometer;
	std::vector<G4ThreeVector> dense;
	for (G4int i = 0; i < m; i++)
	{
		const G4ThreeVector& p0 = ctrl[(i - 1 + m) % m];
		const G4ThreeVector& p1 = ctrl[i];
		const G4ThreeVector& p2 = ctrl[(i + 1) % m];
		const G4ThreeVector& p3 = ctrl[(i + 2) % m];
		const G4int k = std::max(2, (G4int)std::ceil((p2 - p1).mag() / fine));
		for (G4int j = 0; j < k; j++)
		{
			const G4double t = (G4double)j / k, t2 = t * t, t3 = t2 * t;
			dense.push_back(0.5 * ((2. * p1) + (-p0 + p2) * t
								   + (2. * p0 - 5. * p1 + 4. * p2 - p3) * t2
								   + (-p0 + 3. * p1 - 3. * p2 + p3) * t3));
		}
	}
	const G4int nd = (G4int)dense.size();
	G4double splineLength = 0;
	for (G4int i = 0; i < nd; i++)
		splineLength += (dense[(i + 1) % nd] - dense[i]).mag();

	// scale the curve so its length is the polyline's contour length
	scale = polyLength / splineLength;
	for (G4int i = 0; i < nd; i++)
		dense[i] *= scale;

	// resample uniformly by arc length, n = round(length / rise), closing on the first point
	const G4int n = std::max(3, (G4int)std::floor(polyLength / rise + 0.5));
	out = ResampleClosedUniform(dense, n);

	// Curvature cap. A spline through beads a few nanometres apart can turn more sharply at
	// a corner than DNA does (a 2 nm radius is 10 deg per base pair, below the persistence
	// length by far), and sharper than a slab of the shell radius can follow. Where the bend
	// between consecutive legs exceeds maxBend, the point is relaxed toward its neighbours and
	// the curve resampled, repeated until the cap holds; the contour length is restored each
	// round. The relaxation moves a point by a fraction of a nanometre and only where needed.
	if (maxBend > 0)
	{
		for (G4int round = 0; round < 200; round++)
		{
			G4double worst = 0;
			std::vector<G4ThreeVector> relaxed = out;
			for (G4int i = 0; i < n; i++)
			{
				const G4ThreeVector lin = (out[i] - out[(i - 1 + n) % n]).unit();
				const G4ThreeVector lout = (out[(i + 1) % n] - out[i]).unit();
				const G4double bend = std::acos(std::max(-1., std::min(1., lin.dot(lout))));
				worst = std::max(worst, bend);
				if (bend > maxBend)
					relaxed[i] = 0.5 * out[i] + 0.25 * (out[(i - 1 + n) % n] + out[(i + 1) % n]);
			}
			if (worst <= maxBend)
				break;
			// restore the contour length, then re-space the points uniformly
			G4double len = 0;
			for (G4int i = 0; i < n; i++)
				len += (relaxed[(i + 1) % n] - relaxed[i]).mag();
			const G4double s = polyLength / len;
			for (G4int i = 0; i < n; i++)
				relaxed[i] *= s;
			scale *= s;
			out = ResampleClosedUniform(relaxed, n);
		}
	}
	return out;
}

namespace
{
	// n points uniformly spaced by arc length along an open polyline, end points kept
	std::vector<G4ThreeVector> ResampleOpenUniform(const std::vector<G4ThreeVector>& curve, G4int n)
	{
		std::vector<G4ThreeVector> out;
		const G4int nd = (G4int)curve.size();
		G4double length = 0;
		for (G4int i = 0; i + 1 < nd; i++)
			length += (curve[i + 1] - curve[i]).mag();
		const G4double step = length / (n - 1);
		out.reserve(n);
		G4double walked = 0;
		G4int seg = 0;
		G4double segLen = (curve[1] - curve[0]).mag();
		for (G4int i = 0; i < n - 1; i++)
		{
			const G4double target = i * step;
			while (walked + segLen < target && seg + 2 < nd)
			{
				walked += segLen;
				seg++;
				segLen = (curve[seg + 1] - curve[seg]).mag();
			}
			const G4double f = segLen > 0 ? (target - walked) / segLen : 0;
			out.push_back(curve[seg] + f * (curve[seg + 1] - curve[seg]));
		}
		out.push_back(curve[nd - 1]);
		return out;
	}

	G4double LargestOpenBend(const std::vector<G4ThreeVector>& pts)
	{
		G4double worst = 0;
		for (size_t i = 1; i + 1 < pts.size(); i++)
		{
			const G4ThreeVector lin = (pts[i] - pts[i - 1]).unit();
			const G4ThreeVector lout = (pts[i + 1] - pts[i]).unit();
			worst = std::max(worst, std::acos(std::max(-1., std::min(1., lin.dot(lout)))));
		}
		return worst;
	}
}

std::vector<G4ThreeVector> TsCapOpenPathBends(
	const std::vector<G4ThreeVector>& path, G4double maxBend, G4double& bendBefore, G4double& bendAfter)
{
	std::vector<G4ThreeVector> out = path;
	const G4int n = (G4int)out.size();
	bendBefore = bendAfter = LargestOpenBend(out);
	if (n < 3 || maxBend <= 0)
		return out;
	for (G4int round = 0; round < 500 && bendAfter > maxBend; round++)
	{
		std::vector<G4ThreeVector> relaxed = out;
		for (G4int i = 1; i + 1 < n; i++)
		{
			const G4ThreeVector lin = (out[i] - out[i - 1]).unit();
			const G4ThreeVector lout = (out[i + 1] - out[i]).unit();
			const G4double bend = std::acos(std::max(-1., std::min(1., lin.dot(lout))));
			if (bend > maxBend)
				relaxed[i] = 0.5 * out[i] + 0.25 * (out[i - 1] + out[i + 1]);
		}
		out = ResampleOpenUniform(relaxed, n);
		bendAfter = LargestOpenBend(out);
	}
	return out;
}

TsPhysicalDNASolidSet TsBuildPhysicalDNASolidSet(
	const TsPhysicalDNAParams& p, G4double theta, G4double phiLocal, G4double dz, G4int key)
{
	TsPhysicalDNASolidSet s;
	s.theta = theta; s.phiLocal = phiLocal; s.dz = dz;

	// End-face normals in the slab frame. The slab axis is the bisector of the two legs, the
	// bending direction u lies at azimuth phiLocal in the slab's x-y plane, and the legs are
	// tilted by theta/2 either side of the axis within the plane spanned by z and u. The low
	// face is perpendicular to the incoming leg, the high face to the outgoing leg.
	const G4double h = theta / 2.;
	G4ThreeVector u(std::cos(phiLocal), std::sin(phiLocal), 0.);
	G4ThreeVector lowNorm  = -std::cos(h) * G4ThreeVector(0, 0, 1) + std::sin(h) * u;
	G4ThreeVector highNorm =  std::cos(h) * G4ThreeVector(0, 0, 1) + std::sin(h) * u;

	std::ostringstream suffix; suffix << "_k" << key;
	const G4double span = p.backboneSpan;
	const G4double a1 = -span / 2.;                                // strand 1 centred on local +x
	const G4double a2 = p.strandSeparation - span / 2.;            // strand 2

	// An unbent slab (the ends of an open path, or a straight run) has end faces perpendicular
	// to its axis, which is a plain G4Tubs; G4CutTubs warns on such normals and asks for one.
	auto tube = [&](const char* name, G4double rmin, G4double rmax, G4double sphi, G4double dphi) -> G4VSolid* {
		if (theta < 1e-9)
			return new G4Tubs(name + suffix.str(), rmin, rmax, dz, sphi, dphi);
		return new G4CutTubs(name + suffix.str(), rmin, rmax, dz, sphi, dphi, lowNorm, highNorm);
	};
	s.base1  = tube("DNA_base1", 0., p.baseRadius, -90. * CLHEP::degree, 180. * CLHEP::degree);
	s.base2  = tube("DNA_base2", 0., p.baseRadius, 90. * CLHEP::degree, 180. * CLHEP::degree);
	s.back1  = tube("DNA_backbone1", p.baseRadius, p.backboneOuter, a1, span);
	s.back2  = tube("DNA_backbone2", p.baseRadius, p.backboneOuter, a2, span);
	s.shell1 = tube("DNA_shell1", p.backboneOuter, p.shellOuter, a1, span);
	s.shell2 = tube("DNA_shell2", p.backboneOuter, p.shellOuter, a2, span);
	return s;
}

std::vector<TsPhysicalDNABasePair> TsBuildPhysicalDNAPlacement(
	const std::vector<G4ThreeVector>& pathIn,
	const TsPhysicalDNAParams& p,
	G4bool closed,
	std::vector<TsPhysicalDNASolidSet>& solidSets,
	G4String& error)
{
	std::vector<TsPhysicalDNABasePair> out;
	error = "";
	std::vector<G4ThreeVector> path = pathIn;
	if (path.size() < 2) { error = "path holds fewer than two points"; return out; }

	// A closed path's last point can sit within a fraction of a rise of the first; a slab
	// there would be a sliver, so it is absorbed into the previous leg.
	if (closed && (path.back() - path.front()).mag() < 0.5 * p.rise && path.size() > 2)
		path.pop_back();

	const G4int n = (G4int)path.size();
	const G4int nLegs = closed ? n : n - 1;

	// Legs. For an open path the two ends reuse their single leg on both sides.
	std::vector<G4ThreeVector> leg(nLegs);
	std::vector<G4double> legLen(nLegs);
	for (G4int i = 0; i < nLegs; i++)
	{
		G4ThreeVector v = path[(i + 1) % n] - path[i];
		legLen[i] = v.mag();
		if (legLen[i] < 1e-9 * CLHEP::nanometer) { error = "zero-length leg in the path"; return out; }
		leg[i] = v.unit();
	}
	auto legIn  = [&](G4int i) { return closed ? leg[(i - 1 + nLegs) % nLegs] : leg[std::max(i - 1, 0)]; };
	auto legOut = [&](G4int i) { return closed ? leg[i % nLegs] : leg[std::min(i, nLegs - 1)]; };
	auto lenIn  = [&](G4int i) { return closed ? legLen[(i - 1 + nLegs) % nLegs] : legLen[std::max(i - 1, 0)]; };
	auto lenOut = [&](G4int i) { return closed ? legLen[i % nLegs] : legLen[std::min(i, nLegs - 1)]; };

	// Slab axes (bisectors) and a parallel-transported reference vector e, perpendicular to
	// the axis, so the twist is measured against a frame that does not jump at bends.
	std::vector<G4ThreeVector> axis(n), eRef(n);
	for (G4int i = 0; i < n; i++)
		axis[i] = (legIn(i) + legOut(i)).unit();
	{
		G4ThreeVector seed = std::fabs(axis[0].x()) < 0.9 ? G4ThreeVector(1, 0, 0) : G4ThreeVector(0, 1, 0);
		eRef[0] = (seed - seed.dot(axis[0]) * axis[0]).unit();
		for (G4int i = 1; i < n; i++)
		{
			G4RotationMatrix r = MinimalRotation(axis[i - 1], axis[i]);
			G4ThreeVector e = r * eRef[i - 1];
			eRef[i] = (e - e.dot(axis[i]) * axis[i]).unit();
		}
	}

	// Twist per base pair. For a closed path the transported frame returns rotated by some
	// angle delta relative to where it started, so the nominal twist is corrected by
	// (delta + 2 pi m) / n, with m chosen to make the correction smallest, and the strands
	// close on themselves.
	G4double twistPerBp = CLHEP::twopi / p.basePairsPerTurn;
	if (closed)
	{
		G4RotationMatrix r = MinimalRotation(axis[n - 1], axis[0]);
		G4ThreeVector eBack = r * eRef[n - 1];
		eBack = (eBack - eBack.dot(axis[0]) * axis[0]).unit();
		G4ThreeVector f0 = axis[0].cross(eRef[0]);
		G4double delta = std::atan2(eBack.dot(f0), eBack.dot(eRef[0]));   // frame mismatch at closure
		// strand angle after n steps: n * twist + delta (relative to eRef[0]) must be 0 mod 2 pi
		G4double total = n * twistPerBp + delta;
		G4double m = std::floor(total / CLHEP::twopi + 0.5);
		twistPerBp -= (total - m * CLHEP::twopi) / n;
	}

	// Bend limit: on the inside of a bend the slab thickness at the shell radius is
	// rise - shellOuter * theta; it must stay positive with margin or neighbours overlap.
	const G4double thetaMax = 0.8 * p.rise / p.shellOuter;

	out.reserve(n);
	solidSets.reserve(n);
	for (G4int i = 0; i < n; i++)
	{
		const G4ThreeVector lin = legIn(i), lout = legOut(i);
		const G4double cosT = std::max(-1., std::min(1., lin.dot(lout)));
		const G4double theta = std::acos(cosT);
		if (theta > thetaMax)
		{
			std::ostringstream ss;
			ss << "bend of " << theta / CLHEP::degree << " deg at base pair " << i + 1
			   << " exceeds the " << thetaMax / CLHEP::degree
			   << " deg a slab of shell radius " << p.shellOuter / CLHEP::nanometer
			   << " nm can follow without overlapping its neighbour";
			error = ss.str();
			out.clear();
			return out;
		}
		const G4double h = theta / 2.;
		const G4double cosH = std::cos(h);

		// end planes along the axis: low at -lenIn/2/cosH, high at +lenOut/2/cosH
		const G4double zLow  = -0.5 * lenIn(i) / cosH;
		const G4double zHigh =  0.5 * lenOut(i) / cosH;
		const G4double shift = 0.5 * (zLow + zHigh);
		const G4double dz    = 0.5 * (zHigh - zLow);

		TsPhysicalDNABasePair e;
		e.bpID  = i + 1;
		e.twist = i * twistPerBp;
		e.theta = theta;
		e.dz    = dz;
		e.shift = shift;

		// frame: x along strand 1 (reference rotated by the twist about the axis), z the axis
		G4ThreeVector f = axis[i].cross(eRef[i]);
		G4ThreeVector X = std::cos(e.twist) * eRef[i] + std::sin(e.twist) * f;
		G4ThreeVector Z = axis[i];
		G4ThreeVector Y = Z.cross(X);
		e.active = G4RotationMatrix(X, Y, Z);
		e.center = path[i] + shift * axis[i];

		// bending direction in the slab frame
		G4ThreeVector u = lout - lin;
		if (theta < 1e-9) { e.phiLocal = 0.; }
		else
		{
			u = (u - u.dot(Z) * Z);
			e.phiLocal = std::atan2(u.dot(Y), u.dot(X));
		}

		// One exact solid set per base pair. Sharing sets across a bend-angle grid was tried
		// first; the two solids meeting on a face then computed that face from slightly
		// different binned angles, and Geant4's overlap check found them crossing by up to
		// a picometre. Exact faces cost one G4CutTubs per solid and base pair, which a plasmid
		// of a few thousand base pairs affords.
		e.key = (G4int)solidSets.size();
		solidSets.push_back(TsBuildPhysicalDNASolidSet(p, theta, e.phiLocal, dz, e.key));
		out.push_back(e);
	}

	// Debug facility: set TSNBIO_DUMP_PHYSICAL_DNA to a filename to write every slab's centre,
	// frame, twist, bend and half-height, so a script can rebuild the solids and test overlap
	// and connectivity independently of Geant4. Capped at 200000 base pairs.
	if (const char* dumpPath = std::getenv("TSNBIO_DUMP_PHYSICAL_DNA"))
	{
		static G4int written = 0;
		const G4int cap = 200000;
		std::ofstream fdump(dumpPath, std::ios::app);
		fdump.precision(12);   // the checker tests shared planes to 1e-6 nm
		for (size_t i = 0; i < out.size() && written < cap; i++, written++)
		{
			const TsPhysicalDNABasePair& e = out[i];
			const G4RotationMatrix& A = e.active;
			fdump << e.bpID << " "
				  << e.center.x()/CLHEP::nanometer << " " << e.center.y()/CLHEP::nanometer << " " << e.center.z()/CLHEP::nanometer << " "
				  << A.xx() << " " << A.yx() << " " << A.zx() << " "     // world direction of local x
				  << A.xy() << " " << A.yy() << " " << A.zy() << " "     // local y
				  << A.xz() << " " << A.yz() << " " << A.zz() << " "     // local z
				  << e.twist << " " << e.theta << " " << e.phiLocal << " "
				  << e.dz/CLHEP::nanometer << " " << e.key << "\n";
		}
	}
	return out;
}
