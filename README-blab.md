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
| `fix/plasmid-chromosome-contents` | `TsScoreDNADamageSBS.cc:49` hardcodes `fChromosomeContents = {1}` and nothing else assigns it, so every per-DNA yield on the plasmid path is wrong by the ratio of real plasmid content to one base pair. | yes, once written |
| `feature/plasmid-sphere-hydration` | A supercoiled plasmid carrying the nucleus model's sphere volumes and a hydration shell, so an absolute anchor measured on naked plasmid can be reproduced in the geometry it was measured in. | probably not, lab-specific |

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
