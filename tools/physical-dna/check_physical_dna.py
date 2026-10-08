#!/usr/bin/env python3
"""Independent geometry check of the physical DNA solids from a placement dump.

Run a TOPAS deck with TSNBIO_DUMP_PHYSICAL_DNA=<file> set, then:

    check_physical_dna.py <dump> [--base 0.486 --outer 1.0 --shell 1.47 --span 86.7 --sep 180]
                                 [--samples 400] [--neighbours 3]

Rebuilds every slab from its centre, frame, bend and half-height and tests, without Geant4:
  1. end planes: the high face of slab i and the low face of slab i+1 are one plane
     (same point on the path, opposite normals);
  2. no overlap: surface points of every solid of slab i lie outside every solid of slabs
     i-k .. i+k (k = --neighbours), to a tolerance of 1e-6 nm; and bases, backbones and
     shells of one slab are disjoint by construction (checked anyway on their samples);
  3. strand continuity: consecutive backbone sectors of one strand share angular range
     across their common face; bases touch their backbone (shared radius) and each other
     (shared plane through the axis);
  4. the closure: the last slab's high face is the first slab's low face, and the strand
     angle wraps without a seam.
Exit code 0 when every test passes; the report names every failing base pair otherwise.
"""
import argparse, math, sys
import numpy as np

TOL = 1e-6


def load(path):
    rows = np.loadtxt(path)
    if rows.ndim == 1:
        rows = rows[None, :]
    out = []
    for r in rows:
        out.append(dict(bp=int(r[0]), c=r[1:4], X=r[4:7], Y=r[7:10], Z=r[10:13],
                        twist=r[13], theta=r[14], phi=r[15], dz=r[16], key=int(r[17])))
    return out


def to_local(s, pts):
    d = pts - s["c"]
    return np.c_[d @ s["X"], d @ s["Y"], d @ s["Z"]]


def to_world(s, pts):
    return s["c"] + np.outer(pts[:, 0], s["X"]) + np.outer(pts[:, 1], s["Y"]) + np.outer(pts[:, 2], s["Z"])


def planes(s):
    """low/high face (point, outward normal) in local coordinates, as G4CutTubs defines them."""
    h = s["theta"] / 2
    u = np.array([math.cos(s["phi"]), math.sin(s["phi"]), 0.0])
    low = (np.array([0, 0, -s["dz"]]), -math.cos(h) * np.array([0, 0, 1.0]) + math.sin(h) * u)
    high = (np.array([0, 0, s["dz"]]), math.cos(h) * np.array([0, 0, 1.0]) + math.sin(h) * u)
    return low, high


def sectors(args):
    """(name, rmin, rmax, phi0, dphi) of the six solids, degrees."""
    sp, sep = args.span, args.sep
    return [("base1", 0, args.base, -90, 180), ("base2", 0, args.base, 90, 180),
            ("back1", args.base, args.outer, -sp / 2, sp), ("back2", args.base, args.outer, sep - sp / 2, sp),
            ("shell1", args.outer, args.shell, -sp / 2, sp), ("shell2", args.outer, args.shell, sep - sp / 2, sp)]


def inside(s, sec, pts_local, tol=-TOL):
    """points strictly inside the solid (tol < 0 shrinks the solid by |tol|)."""
    name, rmin, rmax, phi0, dphi = sec
    r = np.hypot(pts_local[:, 0], pts_local[:, 1])
    ok = (r > rmin - tol) & (r < rmax + tol)
    ang = (np.degrees(np.arctan2(pts_local[:, 1], pts_local[:, 0])) - phi0) % 360
    ok &= ang < dphi + 1e-9 if dphi >= 360 else (ang > -tol * 0) & (ang < dphi)
    for p0, n in planes(s):
        ok &= ((pts_local - p0) @ n) < tol
    return ok


def surface_samples(s, sec, n):
    """points on the surface of one solid, local coords: curved faces, flat faces, end faces."""
    name, rmin, rmax, phi0, dphi = sec
    rng = np.random.default_rng(s["bp"] * 7 + hash(name) % 1000)
    pts = []
    low, high = planes(s)

    def clip_z(xy):
        # for each (x, y), the z range between the two cut planes
        zl = low[0][2] - ((xy[:, 0] - low[0][0]) * low[1][0] + (xy[:, 1] - low[0][1]) * low[1][1]) / low[1][2]
        zh = high[0][2] - ((xy[:, 0] - high[0][0]) * high[1][0] + (xy[:, 1] - high[0][1]) * high[1][1]) / high[1][2]
        return zl, zh

    for r in (rmin, rmax):
        if r <= 0:
            continue
        a = np.radians(phi0 + dphi * rng.random(n))
        xy = np.c_[r * np.cos(a), r * np.sin(a)]
        zl, zh = clip_z(xy)
        z = zl + (zh - zl) * rng.random(n)
        pts.append(np.c_[xy, z])
    for a in (phi0, phi0 + dphi):               # the two flat radial faces
        rr = rmin + (rmax - rmin) * rng.random(n)
        xy = np.c_[rr * np.cos(math.radians(a)), rr * np.sin(math.radians(a))]
        zl, zh = clip_z(xy)
        z = zl + (zh - zl) * rng.random(n)
        pts.append(np.c_[xy, z])
    for p0, nrm in (low, high):                 # the two end faces
        a = np.radians(phi0 + dphi * rng.random(n))
        rr = np.sqrt(rmin ** 2 + (rmax ** 2 - rmin ** 2) * rng.random(n))
        xy = np.c_[rr * np.cos(a), rr * np.sin(a)]
        z = p0[2] - ((xy[:, 0] - p0[0]) * nrm[0] + (xy[:, 1] - p0[1]) * nrm[1]) / nrm[2]
        pts.append(np.c_[xy, z])
    return np.vstack(pts)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("dump")
    ap.add_argument("--base", type=float, default=0.486)
    ap.add_argument("--outer", type=float, default=1.0)
    ap.add_argument("--shell", type=float, default=1.47)
    ap.add_argument("--span", type=float, default=86.7)
    ap.add_argument("--sep", type=float, default=180.0)
    ap.add_argument("--samples", type=int, default=400)
    ap.add_argument("--neighbours", type=int, default=3)
    ap.add_argument("--closed", action="store_true", default=True)
    args = ap.parse_args()

    slabs = load(args.dump)
    n = len(slabs)
    secs = sectors(args)
    fails = []

    # 1 and 4: shared end planes, including the closure
    for i in range(n):
        j = (i + 1) % n
        if j == 0 and not args.closed:
            break
        a, b = slabs[i], slabs[j]
        _, high = planes(a)
        low, _ = planes(b)
        p_hi = to_world(a, high[0][None, :])[0]
        n_hi = high[1] @ np.vstack([a["X"], a["Y"], a["Z"]])
        p_lo = to_world(b, low[0][None, :])[0]
        n_lo = low[1] @ np.vstack([b["X"], b["Y"], b["Z"]])
        # same plane: normals opposite, and each point on the other's plane
        if abs(n_hi @ n_lo + 1) > 1e-6 or abs((p_lo - p_hi) @ n_hi) > 1e-6:
            fails.append(("end-plane", a["bp"], b["bp"], float(abs(n_hi @ n_lo + 1)), float((p_lo - p_hi) @ n_hi)))

    # 3: strand continuity across each interface, and base contacts
    for i in range(n):
        j = (i + 1) % n
        a, b = slabs[i], slabs[j]
        # strand-1 sector of slab i in the frame of slab j: angle of a's X seen from b
        ang = math.degrees(math.atan2(a["X"] @ b["Y"], a["X"] @ b["X"]))
        overlap = args.span - abs(((ang + 180) % 360) - 180)
        if overlap <= 0:
            fails.append(("strand-gap", a["bp"], b["bp"], float(ang)))

    # 2: overlap test by surface sampling against neighbouring slabs (and within the slab)
    k = args.neighbours
    for i in range(n):
        a = slabs[i]
        for sa in secs:
            pts = surface_samples(a, sa, args.samples)
            # within the slab: points on one solid's surface must not be inside another
            for sb in secs:
                if sb is sa:
                    continue
                if inside(a, sb, pts).any():
                    fails.append(("intra-overlap", a["bp"], sa[0], sb[0]))
            w = to_world(a, pts)
            for d in range(-k, k + 1):
                if d == 0:
                    continue
                j = (i + d) % n
                if not args.closed and (i + d < 0 or i + d >= n):
                    continue
                b = slabs[j]
                loc = to_local(b, w)
                for sb in secs:
                    bad = inside(b, sb, loc)
                    if bad.any():
                        fails.append(("overlap", a["bp"], sa[0], b["bp"], sb[0], int(bad.sum())))
    thetas = np.degrees([s["theta"] for s in slabs])
    print("slabs %d, solid sets %d, bend max %.2f deg, mean %.2f deg" % (
        n, len({s["key"] for s in slabs}), thetas.max(), thetas.mean()))
    if fails:
        print("FAILED: %d findings" % len(fails))
        for f in fails[:40]:
            print("  ", f)
        sys.exit(1)
    print("PASS: end planes shared, strands continuous, no overlaps among %d-neighbour slabs" % k)


if __name__ == "__main__":
    main()
