# BLAB fork of TOPAS-nBio

Fork of `topas-nbio/TOPAS-nBio`, maintained because the recalibration campaign
(`research/projects/topas-nbio-recalibration/`) keeps finding defects that block
measurements and that are worth carrying as tracked branches rather than as local edits on a
server.

`main` tracks upstream and carries no lab changes, so a diff against `upstream/main` is
always the honest statement of what we have altered.

## Branches

| Branch | What it does | Upstreamable |
|---|---|---|
| `fix/g4-11.3-species-names` | Restores hydroxyl damage on Geant4 11.3.x, which renamed radical species by prefixing a degree sign. `TsScoreDNADamageSBS` compared raw strings and silently lost every OH-mediated break; `TsIRTConfiguration` already carried the translation for the IRT path. | yes, ready |
| `fix/nucleus-parameter-names` | Three parameter-name defects in `TsNucleus`: `HydrationShellThickness` was unsettable (check and read named different parameters), `AddBackbones` was silently ignored (code asked for the singular), and `CheckOverlap` exposed the member-variable prefix. | yes, ready |
| `fix/gvalue-ratio-of-sums-estimator` | Both G-value scorers averaged `100*n/E` over events instead of dividing summed molecules by summed energy. The per-event ratio has a small, fluctuating denominator, so it is biased high and its variance is set by whichever event deposited least: eight seeds of the shipped step-by-step example gave G(H2O2) at 1 us of 0.40 to 0.48 on seven and 8.96 on the eighth. Adds ratio-of-sums columns with a delta-method uncertainty, appended so existing parsers keep working, and changes both scorers identically so an IRT-versus-step-by-step comparison never becomes a comparison of two estimators. | yes, ready |
| `fix/reaction-type-inert-on-sbs` | A deck's per-reaction `ReactionType` is silently inert on the step-by-step path and this says so. Partial diffusion control needs a per-encounter probability; `G4DNAMolecularReaction` reads only a radius, and `GetProbability()` has no consumer in the dna module outside `G4DNAMolecularDissociation`. Calling `SetReactionType(1)` there does not add partial diffusion control -- it overwrites the Smoluchowski radius `kobs/(4 pi D_sum N_A)` with a van der Waals contact radius and discards the rate, which afterwards enters nothing. Two earlier commits on this branch did exactly that and are reverted by the third; the branch now leaves the types unapplied, as before, and reports them. | yes, ready |
| `fix/background-scavenger-inert-on-g4-11.3` | A `BackgroundReaction` scavenger on the step-by-step path removed nothing on Geant4 11.3.x: `TsEmDNAChemistry` matched the scavenged species by name, and "OH" never equalled the renamed "°OH", so the first-order reaction was never registered. Silent, because the parser has already printed "Set scavenging capacity" by then; measured G(OH) at 1 us was identical to four digits from 1e5 to 3.2e6 /s and to the unscavenged run. Matches on the definition instead. Also counts the scavenged product in the molecule counter, without building a track for it, so G(Product) against capacity is readable from a chemistry-only run. | yes, ready |
| `fix/plasmid-chromosome-contents` | `TsScoreDNADamageSBS.cc:49` hardcodes `fChromosomeContents = {1}` and nothing else assigns it, so every per-DNA yield on the plasmid path is wrong by the ratio of real plasmid content to one base pair. | yes, once written |
| `feature/plasmid-sphere-hydration` | Gives `TsPlasmidSphereDNA` the nucleus's hydration shell, so quasi-direct damage can be measured on naked plasmid in the geometry the direct and indirect anchors already use. The plasmid built only bases and backbones, so the charge-transfer channel was silently zero there. The wedge solids move into the shared `TsSphereDNAPlacement` beside the placement, and the nucleus takes them from there, so there is one shell, not two that agree by hand. Same parameter names under the plasmid component (`AddHydrationShell`, `HydrationShellThickness`, `HydrationShellPhiSpan`), default off. Nucleus example and plasmid-without-shell outputs byte-identical to the previous integration build. | probably not, lab-specific |
| `feature/attack-counters-by-species` | Splits the attack counters of `feature/scavenging-counters` by attacking species. `TsScoreDNADamageSBS` lets both the hydroxyl radical and the solvated electron attack a base but only the hydroxyl attack a backbone, so the all-species backbone share `bb / (bb + base)` is not the hydroxyl-only partition Michalik 1995 describes; a 1e8 /s hydroxyl scavenger cut backbone attacks 31-fold and base attacks 1.3-fold, so most base attacks are electron attacks. Adds `Scavenged_Backbone_OH`, `Scavenged_Base_OH` and `Scavenged_Base_eaq`, appended to the ntuple and written to `<OutputFile>_scavenged_species.csv`; the existing `_scavenged.csv` and every other output stay byte-identical at fixed seed. On the unscavenged naked plasmid with both scavenging probabilities 1, the backbone share is 0.315 +/- 0.017 over all species and 0.645 +/- 0.025 hydroxyl-only, so the two are not interchangeable. | yes, with the counters branch |
| `feature/molecule-positions-scorer` | `TsScoreMoleculePositionsAtTime`, an ntuple scorer writing one row per live molecule (species, position, track and parent track, parent vertex) at a settable chemical time, default 1 ps, from the post-time-step hook. Written to measure how the electron elastic model changes spur geometry: the count of radicals within 1 nm of a hydroxyl at 1 ps ranks the five elastic models exactly as their hydroxyl decay factors do, while the electron-to-parent-cation distance does not. Rows are recorded before the zero-time contact reactions, so counts sit a few percent above the G-value scorer's 1 ps yields by exactly those reactions. | yes, as a diagnostic |
| `feature/chemistry-max-time-step` | `TsDNAMaxTimeStep`, a process on every diffusing molecule that caps the step-by-step chemistry time step. The scheduler steps each molecule to its nearest reactive partner with no ceiling (`G4Scheduler::SetMaxTimeStep` is assigned and overwritten before use), so a radical with no partner nearby jumps the whole remaining interval in one straight segment, and Brownian transport checks DNA volumes only along that segment. Pure-water G values are untouched; encounters with DNA are undercounted once radicals have to diffuse to it. Unscavenged 8-plasmid lattice at 30 us: hydroxyl attacks on DNA 48, 86, 118, 178 for no cap, 10 ns, 1 ns, 0.1 ns, with 0.03 ns at 189, so 0.1 ns is converged to about 10%. `Ch/<name>/MaximumTimeStep` sets the cap; with `MaximumTimeStepVolumePrefix` (and `MaximumTimeStepMargin`) it applies inside the named envelopes only and the step outside is d^2 / (8 D) to the nearest of them, which reproduced the global 0.1 ns result (168 against 178) at a fifth of the wall time. Off unless set; an uncapped deck's output is byte-identical to the previous integration build. Geant4's own safety-based limit (`G4DNABrownianTransportation::SpeedLevel(0)`) was tried first and lowered the attacks (4 against 48), so it is not exposed. | probably, as an option; the undercount is a Geant4-DNA transport property |
| `feature/plasmid-conformation-builder` | `tools/plasmid-conformation/`: a closed wormlike-chain sampler (`plasmid_mc`) and a builder that write supercoiled plasmid path files for `TsPlasmidSphereDNA`, which resamples any closed polyline at 0.34 nm per base pair, so these files can also feed the nucleus. The shipped `pBR322_a.xyz` is a compact blob (4350 bp, Rg 43 nm, writhe 0.92 of the linking-number change, 8.7 nm between strands), whereas supercoiled pBR322 at sigma -0.06 forms an interwound plectoneme. Eight sampled conformations (`pBR322_lab_01..08`) averaged Rg 61 ± 8 nm, extent 168 ± 21 nm, strand distance 12.1 ± 0.4 nm and writhe 0.57 ± 0.02 of the linking-number change, and one relaxed circle for bracketing ships with them; all fit the 456 nm lattice cell of eight plasmids at 50 ug/cm3. The core diameter (8.0 nm) was set from the pUC18 neutron-scattering superhelix diameter (Hammermann 1998), not from any damage curve; the writhe fraction sits below the 0.7 to 0.8 of high-salt literature and was reported, not tuned. Path files only, so no C++ change and no rebuild. At matched deck settings (8 plasmids, 3.5e5 and 1e6 /s) strand-break yield rose 1.36 ± 0.10 and 1.45 ± 0.17 times over the blob; the relaxed circle gave 1.69 ± 0.12 and 1.71 ± 0.19. | no, lab-specific |
| `feature/physical-dna-solids` | `TsPhysicalDNASolids`, a second DNA volume model for `TsPlasmidSphereDNA` (`Ge/<name>/DNAModel = "Physical"`; default `"Sphere"`, byte-identical). The sphere solids carry 0.31 of the base-pair mass with open water between them, measured against Milligan 1993 on the plasmid as a direct plateau at 0.35 of the measurement, an indirect yield that rose 1.7-fold when the solids were scaled up, and backbone attacks not converged in the chemistry time step (1.5-fold higher at 0.01 ns than at 0.1 ns) because a 0.27 nm orb is stepped over. Each base pair becomes a 0.34 nm slab of a 1.0 nm duplex: a 0.486 nm core split into two half-cylinder bases, two 86.7° backbone sectors to 1.0 nm that rotate 34.3° per base pair into continuous helical ribbons, open grooves, and 10-water hydration-shell sectors to 1.47 nm; DNA material at 1.7 g/cm3 holds the 660 Da inside the duplex with the grooves open. Slabs are `G4CutTubs` whose end faces are the planes through the leg midpoints, shared with the neighbours; one exact solid set per base pair (a binned shared table left picometre crossings). Bead paths are smoothed by a closed Catmull-Rom spline, scaled to their contour length, resampled by arc length and relaxed to at most 8° per base pair (`PhysicalMaxBendPerBasePair`); a parallel-transported frame with the closure mismatch spread as twist closes the strands without a seam. Verified: eight lab conformations pass `tools/physical-dna/check_physical_dna.py` (shared planes, continuous strands, no overlaps among three neighbours); Geant4 `CheckOverlaps` zero on closed circles of 2.6, 5.4 and 20 nm radius; 1 to 64 plasmids build linearly (11 to 480 s, 0.26 to 4.0 GB). `TsScoreDNADamageSBS` gains `IndirectPartitionMode = "Chemical"` (`ProbabilityOfScavengingOnDNA`, `HydroxylBaseFraction` 0.68): an encounter with any solid of a base pair is accepted once and the moiety drawn from the chemistry, so continuous backbones do not shield the bases; default `"Geometric"` unchanged. Volume names, copy numbers, strand materials, SDD and counters unchanged. | the solids yes, as a model option; the partition mode yes |

| `feature/chemistry-near-dna-step` | Refines `TsDNAMaxTimeStep` near the DNA solids: with `Ch/<name>/MinimumTimeStepNearDNA` set, the base-pair positions inside the named envelopes (their `Base1_` daughters) are collected once into a 2 nm grid, and a molecule whose coarse step could reach a base pair steps by (d − `NearDNARadius`, default 1.5 nm)² / (8 D), never below the minimum; every other molecule keeps the cap and the envelope rule. Why: on the sphere model at the 0.1 ns cap, hydroxyl attacks at 1e9 /s rose 1.66-fold when the cap was lowered to 0.01 ns everywhere, without flattening, at seven times the wall time. On the physical solids, minima of 0.01 and 0.003 ns gave 204 and 218 attacks at 1e8 /s and 58 and 58 at 1e9 /s (three decks, 20,000 histories), converged, at 1.4 to 1.8 times the plain-cap wall time; 1.5 and 2.1 times the sphere model's capture at identical seeds. Off unless the minimum is given. | yes, with the cap branch |

## Upstream pull requests

Each defect goes upstream as its own single-commit branch off `main`, named `pr/<slug>`, carrying
none of the lab-only files. The stack on `blab/integration` stays the lab build. Bodies and the
matching Discussion replies are drafted in `pr-drafts/`; the PI opens the pull request and posts
the reply.

| Branch | Commit | From | Status |
|---|---|---|---|
| `pr/g4-11.3-species-names` | `7fb5e3d` | `fix/g4-11.3-species-names` | prepared 2026-09-28; built clean on upstream `main`, nrtest 19/19; not opened. Answers Discussions #125, #126 (#96 related) |
| `pr/background-scavenger-g4-11.3` | `216a55a` | `fix/background-scavenger-inert-on-g4-11.3` | prepared 2026-09-28; built clean, nrtest 19/19; not opened. Answers Discussion #79 |
| `pr/nucleus-parameter-names` | `992f23d` | `fix/nucleus-parameter-names` | prepared 2026-09-28; built clean, nrtest 19/19; not opened. `TsNucleus.cc` is CRLF, so `git am` needs `--keep-cr` |

## Working rules

- One branch per defect or feature, never a combined branch: each has to be separately
  reviewable and separately upstreamable.
- **`blab/integration` is what we build from**, and it only ever merges the fix branches; no
  work is committed to it directly. A campaign needs several fixes at once while upstream
  wants them one at a time, and this is the seam between those two demands. Its merge list is
  the exact statement of which defects a given result was measured under, so a campaign
  manifest records the integration commit, never a fix-branch commit.
- Every fix commit states the measurement that revealed it, with numbers. A defect nobody
  can reproduce is a defect nobody will merge.
- Rebase on `upstream/main` before proposing anything upstream; never merge upstream into a
  fix branch, which buries the diff.
- Builds live outside this tree. `/opt/blab/topas-nbio-fork/` on blade-server is an rsync
  target, not a checkout, so the checkout stays clean and its commit is the provenance a
  campaign manifest records.
- LAB-PRIVATE campaign generators and decks never enter this repository. They live in
  `research/projects/topas-nbio-recalibration/`.

## Attribution

Commits here carry the PI as author and no other attribution. The repository is public and
upstreamable, so its history reads as lab work. Anything written with assistance is still
reviewed and owned by the named author before it lands.

## Provenance

Forked at upstream `4c95802`, which is the commit the C0 verification campaign ran on, so
every calibration result to date is reproducible from `main`.
