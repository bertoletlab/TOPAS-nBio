# Plasmid conformation builder

Generates closed-polyline path files for `TsPlasmidSphereDNA` (and any other user of
`TsSegmentPathToBasePairs`, including the nucleus). The shipped `pBR322_a.xyz` is a compact blob
(Rg 43 nm, 145 points 10 nm apart). The files built here are supercoiled conformations sampled from a
closed wormlike-chain model, so the DNA exposure to radicals reflects a real plasmid.

## Procedure

1. Build the sampler: `c++ -O2 -std=c++17 -o plasmid_mc plasmid_mc.cpp`
2. Sample conformations (one seed per file):
   `./plasmid_mc --bp 4363 --sigma -0.06 --core-diameter 8.0 --moves 2000000 --seed 31 --out mc/pbr322_s31.xyz`
3. Convert to path files and write diagnostics:
   `python3 build_conformations.py --bp 4363 --sigma -0.06 --out-dir conformations --prefix pBR322_lab_ mc/pbr322_s*.xyz`
4. Draw them to scale against the lattice cell:
   `python3 draw_conformations.py --box-nm 456 --out fig.png "shipped blob=<pBR322_a.xyz>" "lab 01=conformations/pBR322_lab_01.xyz"`

`build_conformations.py` passes a periodic spline through the vertices, rescales the arc length to
`bp x 0.34 nm`, and reports Rg, extent, writhe, paired fraction, superhelix diameter, and end
loops. `fits_box` compares the extent with `--max-box-nm` (default 440 nm, below the 456 nm lattice
pitch for N = 8 at 50 ug/cm3).

## Model

Energy per conformation: bending `sum (P/l0)(1 - cos theta)` plus twist-writhe
`(2 pi^2 C / L)(dLk - Wr)^2`, hard-core exclusion at `d_eff`, writhe by the Klenin-Langowski
pair formula, crankshaft moves. Defaults: persistence P = 50 nm, torsional C = 49 nm, helical repeat
10.5 bp/turn, segment 5 nm. Negative supercoiling gives a right-handed interwound plectoneme (Wr < 0).

## Inputs and calibration

- P, C, h, sigma: literature values, graded in
  `research/projects/topas-nbio-recalibration/refdata/plasmid_conformation_facts.md`.
- `d_eff` (8.0 nm) was set from pUC18: the simulated superhelix diameter must match the small-angle
  neutron scattering measurement of Hammermann 1998 (13.8 nm measured, 14.1 nm simulated at 20 mM
  salt). It was not tuned to any DNA-damage curve.
- Validation of the sampler: coil analytic check, planar ellipse (Wr = 0), independent Gauss-integral
  writhe, start-independence.

## Reported, not tuned

Low-salt Wr/dLk came out at about 0.55, against 0.7 to 0.8 in the high-salt literature. The 8 files
in `conformations/` (`pBR322_lab_01..08.xyz`, summary in `pBR322_lab_summary.json`) carry this value
as sampled.

## Use in a deck

Pass the file basename to `gen_c4.py --xyz` (comma list for several; cycled over lattice sites). The
generator derives bp and extent from the path files.
