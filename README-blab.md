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
