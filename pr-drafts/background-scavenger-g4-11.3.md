# PR draft: pr/background-scavenger-g4-11.3

Branch: `pr/background-scavenger-g4-11.3` (one commit on top of `main` 4c95802)
Files: `processes/TsEmDNAChemistry.cc`, `processes/TsDNAFirstOrderReaction.cc` (+30/-15)

## Title

Attach the background scavenger on Geant4 11.3 and count what it captures

## Body

### Defect

`TsEmDNAChemistry::ConstructProcess` attached a `BackgroundReaction` scavenger to a species by comparing the molecule definition's name with the name resolved from the deck. Geant4 11.3 renamed the radical definitions with a degree-sign prefix, so `G4OH::Definition()` is now named `°OH` while the deck, the parser, and `G4MoleculeTable` still say `OH`. The comparison therefore never matched hydroxyl, the `TsDNAFirstOrderReaction` was never registered, and a `BackgroundReaction` scavenger on the step-by-step path removed nothing. The miss was silent, since the parser had already printed `Set scavenging capacity for molecule OH` before the comparison ran.

A second, related gap sat in `TsDNAFirstOrderReaction::PostStepDoIt`: once the reaction fired, each product was handed only to `G4DNADamage`, the DNA bookkeeping, so the G-value scorer never saw it and G(Product) read zero at every time.

### Symptom, measured

Setup: `examples/scorers/SBSGetGValue/GvalueRevisedPhysicsChemistry.txt` (TsEmDNAPhysics + TsEmDNAChemistry, 1 MeV electrons, primary killed after 10 keV loss, chemistry to 1 us), 20 histories, one thread, three seeds, with and without a hydroxyl scavenger declared as

```
sv:Ch/TOPASChemistry/BackgroundReaction/Hydroxyl/Scavenger/Products = 1 "Product"
d:Ch/TOPASChemistry/BackgroundReaction/Hydroxyl/Scavenger/ScavengingCapacity = 3.2e6 /s
b:Ch/TOPASChemistry/BackgroundReaction/Hydroxyl/Scavenger/CompatibleWithStepByStep = "True"
```

Values are the scorer's `GValue` column at 999999 ps (molecules per 100 eV, mean of the per-event ratios, so a biased estimator; the numbers are quoted as the scorer wrote them).

| Build | Seed | G(°OH), no scavenger | G(°OH), 3.2e6 /s | G(PRODUCT), 3.2e6 /s |
|---|---|---|---|---|
| before (main) | 1 | 4.346 | 4.346 | 0 |
| before (main) | 2 | 5.377 | 5.377 | 0 |
| before (main) | 3 | 5.286 | 5.286 | 0 |
| after (this PR) | 1 | 4.346 | 0.299 | 4.160 |
| after (this PR) | 2 | 5.377 | 0.569 | 4.777 |
| after (this PR) | 3 | 5.286 | 0.103 | 4.917 |

Before the fix, the scavenged run reproduced the unscavenged run to every printed digit at every seed, which is the signature of a reaction that was never attached. After the fix, the scavenger removed 94% of the hydroxyl by 1 us (seed mean 5.00 to 0.32), and the captured amount appeared as G(PRODUCT). The unscavenged path did not change: its values are identical before and after.

### Fix

(i) `TsEmDNAChemistry.cc`: resolve the scavenged species through `G4MoleculeTable::Instance()->GetConfiguration(name)`, which still answers to the deck's spelling, and compare the resolved definition with the iterated one. Only the comparison was wrong; the lookup that followed it already used the table.

(ii) `TsDNAFirstOrderReaction.cc`: report each product to `G4VMoleculeCounter` at the reaction time and position, without building a track for it. The product stays inert, as it did before, and the G-value scorer can now count it. The commented-out block that built a product track was removed.

### Verification

- Built against Geant4 11.3.2 (G4EMLOW 8.6.1) and OpenTOPAS 4.2.3 with gcc 11.4.0 on Ubuntu 22.04, with the cmake configuration of the `main` build. The build log carried 87 unique warning lines before the patch and 87 after, all pre-existing; this patch introduced none.
- qi-topas-nbio (nrtest, latest tag, run natively against this build): all 19 tests passed, and `nrtest compare` against the unpatched build reported every output identical, so the suite does not exercise the changed path.
- Example runs: table above.

### Related

Discussion #79 reported scavenger concentration having no effect on the reported yields. That report used the IRT path (`IRTDamageToDNAPlasmid`), which resolves scavengers in `TsIRTConfiguration` and is not touched here; on the step-by-step path this PR is the cause and the fix.

## Discussion reply draft (#79)

We hit what looks like this symptom on the step-by-step chemistry path (`TsEmDNAChemistry`) with Geant4 11.3 and traced it to a name comparison: 11.3 renamed the radical definitions with a degree-sign prefix, so the scavenger's `OH` never matched `°OH`, the first-order reaction was never registered, and the run printed "Set scavenging capacity" and then scavenged nothing. In the shipped `SBSGetGValue` example, G(OH) at 1 us was identical to every digit with and without a 3.2e6 /s hydroxyl scavenger; after the fix in PR [number] the scavenger removed 94% of the hydroxyl and its product was counted. Your deck used the IRT path, which resolves scavengers in `TsIRTConfiguration` and is a separate code path, so if the yields stayed flat there as well the cause may differ; a quick check is whether `/Gy/Da` changes at all between 0.001 M and 10 M on Geant4 11.2, where the names did not change.
