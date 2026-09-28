# PR draft: pr/nucleus-parameter-names

Branch: `pr/nucleus-parameter-names` (one commit on top of `main` 4c95802)
Files: `geometry/nucleusModel/TsNucleus.cc` (+20/-3; the file keeps its CRLF line endings)

## Title

Fix three parameter-name defects in TsNucleus

## Body

### Defect

Three parameters in the `TsNucleus::Construct` block did not do what their names promised.

(i) `HydrationShellThickness`: the existence check asked for `Ge/Nucleus/HydrationShellThickness` while the value read asked for `Ge/Nucleus/fHydrationShellThickness`. A deck setting the documented name passed the check and then failed the read, so the run stopped with a parameter error. The thickness could not be set at all.

(ii) `AddBackbones`: every shipped example, including `examples/scorers/SBSDamageToDNANucleus/ExampleNucleusDNADamage.txt`, sets the plural, while the code asked for the singular `AddBackbone`. The parameter was ignored without any message, so `AddBackbones = "false"` still built backbones. The default is true and every example sets true, which is why nothing had noticed.

(iii) `CheckOverlap`: the member-variable prefix had leaked into the user-facing name, so a deck had to write `Ge/Nucleus/fCheckOverlap`. No other parameter in the class carries the prefix.

### Symptom, measured

Setup: the `SBSDamageToDNANucleus` geometry (sphere DNA model, 4.65 um nucleus, Hilbert layer 4, repeat 30) with `g4em-dna` physics only, three 5 MeV protons, and two energy-deposit scorers filtered by material: `G4_BackboneMaterial` and `G4_WATER_MODIFIED` (the hydration shell). One seed throughout, so the deposits are directly comparable.

| Deck | Before (main) | After (this PR) |
|---|---|---|
| `HydrationShellThickness = 0.5 nm` | run aborted: `Parameter name: Ge/Nucleus/fHydrationShellThickness has not been defined.` | ran; shell deposit 1.496e-3 MeV (default 0.16 nm gave 5.963e-4 MeV), backbone deposit 1.382e-4 MeV (default 7.797e-4 MeV) |
| `AddBackbones = "false"` | backbone-material deposit 7.797e-4 MeV, identical to the default deck, so backbones were still built | backbone-material deposit 0 |
| `CheckOverlap = "true"` | no overlap check ran (zero `checking CountFibers` lines in the log) | the check ran (28 `checking CountFibers` lines) |

The last row also turned up something the maintainers may want to know: with the check switched on, the stock example geometry reported `Detected Overlap, see file overlap.txt for more details.` for its fibers. That finding is outside this PR, which only makes the switch reachable.

### Fix

Read `CheckOverlap` and `AddBackbones` first and fall back to the old spellings (`fCheckOverlap`, `AddBackbone`) so that any existing deck keeps working; read `HydrationShellThickness` under the name the existence check already uses. One file, +20/-3, no behavior change for decks that set none of the three.

### Verification

- Built against Geant4 11.3.2 (G4EMLOW 8.6.1) and OpenTOPAS 4.2.3 with gcc 11.4.0 on Ubuntu 22.04, with the cmake configuration of the `main` build. The build log carried 87 unique warning lines before the patch and 87 after, all pre-existing; this patch introduced none.
- qi-topas-nbio (nrtest, latest tag, run natively against this build): all 19 tests passed, and `nrtest compare` against the unpatched build reported every output identical, so the suite does not exercise the changed path.
- Example runs: table above.

### Related

No upstream Discussion reports these three names directly (searched for `HydrationShellThickness`, `AddBackbones`, `fCheckOverlap`). Discussions #125 and #126 carry decks that set `AddBackbones`, which the code ignored until now.

## Discussion reply draft

No matching Discussion exists for this defect, so no reply is drafted. If one appears, the one-paragraph form would be: the three names in the `TsNucleus` parameter block were read under different spellings than the documented ones, so `HydrationShellThickness` aborted the run, `AddBackbones` was ignored, and `CheckOverlap` had to be written with an `f` prefix; PR [number] reads the documented names and keeps the old spellings as fallbacks.
