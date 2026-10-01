#!/usr/bin/env python3
"""Turn plasmid_mc vertex output into path files for TsPlasmidSphereDNA, with diagnostics.

For each input polygon (one "x y z" nm line per vertex, from plasmid_mc) this script
  1. passes a periodic cubic spline through the vertices,
  2. resamples it at 1 nm and rescales so the arc length is exactly N_bp * 0.34 nm,
  3. centers it at the origin and writes an xyz path (the format the shipped pBR322_a.xyz uses),
  4. reports the derived geometry: Rg, bounding box, writhe, paired fraction, superhelix
     diameter, end-loop count.

Usage: build_conformations.py --bp 4363 --out-dir conformations --prefix pBR322_ --sigma -0.06 \
           mc/pbr322_s1.xyz mc/pbr322_s2.xyz ...
"""
import argparse
import json
import os
import sys

import numpy as np
from scipy.interpolate import CubicSpline
from scipy.spatial import cKDTree

RISE = 0.34


def periodic_path(v, n_bp):
    """Resample the closed polygon v with a periodic cubic spline at ~1 nm, exact length."""
    L = n_bp * RISE
    per = np.r_[0, np.cumsum(np.linalg.norm(np.diff(np.vstack([v, v[:1]]), axis=0), axis=1))]
    cs = CubicSpline(per, np.vstack([v, v[:1]]), bc_type="periodic")
    s = np.linspace(0, per[-1], int(10 * per[-1]), endpoint=False)
    p = cs(s)
    for _ in range(3):  # spline arc length differs from the polygon perimeter; rescale to L
        seg = np.linalg.norm(np.diff(np.vstack([p, p[:1]]), axis=0), axis=1)
        p = p * (L / seg.sum())
    m = int(round(L))
    seg = np.linalg.norm(np.diff(np.vstack([p, p[:1]]), axis=0), axis=1)
    arc = np.r_[0, np.cumsum(seg)]
    t = np.linspace(0, arc[-1], m, endpoint=False)
    q = np.column_stack([np.interp(t, arc, np.append(p[:, k], p[0, k])) for k in range(3)])
    seg = np.linalg.norm(np.diff(np.vstack([q, q[:1]]), axis=0), axis=1)
    q = q * (L / seg.sum())
    return q - q.mean(axis=0)


def writhe(q):
    n = len(q)
    mid = 0.5 * (q + np.roll(q, -1, axis=0))
    d = np.roll(q, -1, axis=0) - q
    s = 0.0
    for i in range(n):
        r = mid[i] - mid
        r2 = (r * r).sum(axis=1)
        r2[i] = 1.0
        t = (np.cross(d[i], d) * r).sum(axis=1) / r2 ** 1.5
        t[i] = 0.0
        s += t.sum()
    return s / (4 * np.pi)


def diagnostics(q, n_bp, sigma, h):
    L = n_bp * RISE
    seg = np.linalg.norm(np.diff(np.vstack([q, q[:1]]), axis=0), axis=1)
    ext = q.max(axis=0) - q.min(axis=0)
    rg = np.sqrt(((q - q.mean(axis=0)) ** 2).sum(axis=1).mean())
    coarse = q[:: max(1, int(round(5.0)))]
    wr = writhe(coarse)
    dlk = sigma * n_bp / h
    tree = cKDTree(q)
    n = len(q)
    nn = np.full(n, np.inf)
    partner = np.full(n, -1)
    pairs = tree.query_ball_point(q, r=40.0)
    for i, js in enumerate(pairs):
        best, bj = np.inf, -1
        for j in js:
            c = min(abs(i - j), n - abs(i - j))
            if c < 40:
                continue
            dd = np.linalg.norm(q[i] - q[j])
            if dd < best:
                best, bj = dd, j
        nn[i], partner[i] = best, bj
    paired = nn < 25.0
    frac = float(paired.mean())
    dmed = float(np.median(nn[paired])) if paired.any() else float("nan")
    # end loops: maximal unpaired runs (cyclic) longer than 30 nm
    un = ~paired
    loops = 0
    if un.any() and not un.all():
        start = int(np.argmax(paired))
        u = np.roll(un, -start)
        run = 0
        for x in u:
            if x:
                run += 1
            else:
                if run >= 30:
                    loops += 1
                run = 0
        if run >= 30:
            loops += 1
    return {
        "n_points": int(n),
        "arc_length_nm": float(seg.sum()),
        "expected_arc_length_nm": L,
        "rg_nm": float(rg),
        "extent_xyz_nm": [float(e) for e in ext],
        "max_extent_nm": float(ext.max()),
        "writhe": float(wr),
        "delta_lk": float(dlk),
        "wr_over_dlk": float(wr / dlk) if dlk else None,
        "paired_fraction_2L_over_l": frac,
        "strand_distance_D_nm": dmed,
        "unpaired_runs_ge_30nm": int(loops),
        "branches_estimate": int(max(0, loops - 2)),
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bp", type=int, required=True)
    ap.add_argument("--sigma", type=float, required=True)
    ap.add_argument("--helical-repeat", type=float, default=10.5)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--prefix", required=True)
    ap.add_argument("--max-box-nm", type=float, default=440.0)
    ap.add_argument("inputs", nargs="+")
    a = ap.parse_args()
    os.makedirs(a.out_dir, exist_ok=True)
    records = []
    for k, f in enumerate(a.inputs, 1):
        v = np.loadtxt(f)
        q = periodic_path(v, a.bp)
        d = diagnostics(q, a.bp, a.sigma, a.helical_repeat)
        name = f"{a.prefix}{k:02d}.xyz"
        np.savetxt(os.path.join(a.out_dir, name), q, fmt="%.4f")
        d.update({"file": name, "source": os.path.basename(f), "fits_box": bool(d["max_extent_nm"] + 3.4 < a.max_box_nm)})
        records.append(d)
        print(json.dumps(d))
    with open(os.path.join(a.out_dir, f"{a.prefix}summary.json"), "w") as fh:
        json.dump(records, fh, indent=1)


if __name__ == "__main__":
    sys.exit(main())
