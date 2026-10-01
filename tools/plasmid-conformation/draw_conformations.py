#!/usr/bin/env python3
"""To-scale projections of plasmid path files: shipped blob against the lab conformations.

Every panel uses the same length scale and shows the lattice cell (box) the plasmid has to fit.
Usage: draw_conformations.py --box-nm 456 --out fig.png label=path.xyz [label=path.xyz ...]
"""
import argparse

import matplotlib.pyplot as plt
import numpy as np


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--box-nm", type=float, default=456.0)
    ap.add_argument("--out", required=True)
    ap.add_argument("items", nargs="+")
    a = ap.parse_args()
    n = len(a.items)
    fig, axes = plt.subplots(2, n, figsize=(2.6 * n, 5.4), squeeze=False)
    half = a.box_nm / 2
    for k, item in enumerate(a.items):
        label, path = item.split("=", 1)
        q = np.loadtxt(path)
        q = q - q.mean(axis=0)
        seg = np.linalg.norm(np.diff(np.vstack([q, q[:1]]), axis=0), axis=1)
        bp = seg.sum() / 0.34
        for row, (i, j, name) in enumerate([(0, 1, "xy"), (0, 2, "xz")]):
            ax = axes[row][k]
            ax.plot(np.r_[q[:, i], q[0, i]], np.r_[q[:, j], q[0, j]], lw=1.0, color="#1f4e79")
            ax.add_patch(plt.Rectangle((-half, -half), a.box_nm, a.box_nm, fill=False, ec="#999999", lw=0.8, ls="--"))
            ax.plot([-half + 20, -half + 120], [-half + 20, -half + 20], color="k", lw=2)
            ax.text(-half + 70, -half + 32, "100 nm", ha="center", fontsize=7)
            ax.set_xlim(-half - 10, half + 10)
            ax.set_ylim(-half - 10, half + 10)
            ax.set_aspect("equal")
            ax.axis("off")
            if row == 0:
                ext = q.max(axis=0) - q.min(axis=0)
                ax.set_title(f"{label}\n{bp:.0f} bp, extent {ext.max():.0f} nm", fontsize=8)
    fig.tight_layout()
    fig.savefig(a.out, dpi=200)


if __name__ == "__main__":
    main()
