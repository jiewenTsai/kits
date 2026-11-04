# kits validation report

Standard: [KIT validation standard](https://github.com/jiewenTsai/kitDesign/blob/main/inst/validation/criteria.md)
(`kitDesign::kit_validation("kitDesign", "criteria")`). Version validated: kits 0.26.0, with the
regression tests and the fix of kits 0.26.1. Date: 2026-11-04.

## Summary

Every number kits writes into a report agrees with an established engine: 32 of 32 comparisons
pass, from 5.4e-14 against dexter to 0.099 against TAM's MML (expected, see below). 30 of 30
metamorphic relations hold. The existing tests caught 18 of 20 planted faults; of the two
survivors, one is an equivalent mutant and the other is now covered by a regression test. This
validation was run after release, so no acceptance criteria were written before generation. One
fault found after release was still open (K49) and is fixed in 0.26.1.

## Acceptance criteria (written before generation)

None were written before the code; this report is retrospective. From the next release on, every
change to an in-scope function lists three to five criteria here before the code is generated.
The criteria that the existing tests imply:

1. Item parameters and their SEs equal dexter's CML estimates to 1e-8.
2. WLE scores and their SEs equal TAM's `tam.wle` with the same item parameters to 1e-8.
3. Linking constants equal the mean/sigma formula of Kolen & Brennan on the screened points.
4. A score table, item table or check never stops on a cohort, wave or category with one person.

## Scope

| Function | Output that reaches a report | Error types asked |
|---|---|---|
| `kits_item` | Item parameters (beta, SE) in the fixed set | Formula, interval, boundary |
| `kits_subj`, `kits_subj_table` | Person scores and SEs, raw-score tables | Formula, interval, boundary |
| `kits_itemfit`, `kits_item_restscore` | Item fit tables | Formula, boundary |
| `kits_information`, `kits_expected_score`, `kits_pcm_prob`, `kits_thresholds`, `kits_esf`, `kits_precision` | Plots and tables of the item report | Formula |
| `kits_linking`, `kits_link_line`, `kits_link_error`, `kits_convert` | Linking constants, linking error, converted scores | Formula, interval, variable mapping |
| `kits_drift`, `kits_dif` | Drift and DIF tables | Formula, boundary |
| `kits_norms`, `kits_norm_*`, `kits_vertical` | Norm tables, growth table | Formula, data management, boundary |

## Engine equivalence

Script: `inst/validation/engine.R` (needs dexter, psychotools, TAM, mirt). Data: `kits_sim` with
fixed seeds and `kitm03` lan. Results: `inst/validation/results/engine_results.csv`. Runtime 9 s.

| Function | Reference engine | Quantity | Max abs. difference | Pass |
|---|---|---|---|---|
| `kits_item` | dexter 1.8.1 `fit_enorm` (CML) | beta; SE | 5.4e-14; 2.8e-11 | yes |
| `kits_item` | psychotools 0.7.7 `pcmodel` (CML) | beta; SE | 1.2e-6; 1.4e-8 | yes |
| `kits_item` | TAM 4.3.25 `tam.mml` (MML), both centred | beta | 0.099 (r = 0.9997) | yes (tolerance 0.15) |
| `kits_subj` | TAM `tam.wle`, fixed parameters | WLE; SE | 4.1e-11; 3.0e-11 | yes |
| `kits_itemfit` (WLE) | mirt 1.47 `itemfit` | infit, outfit | 2.2e-16 | yes |
| `kits_itemfit` (conditional) | iarm 0.4.3, stored values | infit, outfit | 2.2e-15 | yes |
| `kits_information` | mirt `testinfo` | I(theta) | 1.1e-14 | yes |
| `kits_expected_score` | mirt `expected.test` | expected score | 3.6e-15 | yes |
| `kits_pcm_prob` | mirt `probtrace` | category probabilities | 3.3e-16 | yes |
| `kits_thresholds` | TAM `tam.threshold` | RP50 thresholds | 8.6e-5 | yes (TAM bisection) |
| `kits_esf` | psychotools `elementary_symmetric_functions` | gamma (relative) | 1.2e-15 | yes |
| `kits_item_restscore` | iarm 0.4.3, stored values | gamma, SE, p | 2.9e-15 | yes |
| `kits_dif` | dexter `DIF`; psychotools LR | chi-square; LR; item DIF | 1.1e-10; 5.1e-10; 4.7e-7 | yes |
| `kits_precision` | mirt and TAM, all 432 patterns | E[WLE], RMSE | 3.1e-12 | yes |
| `kits_linking`, `kits_drift` | Hand formula (mean/sigma, one 2 SD screen) | intercept, slope | 0 | yes |

No reference engine: the norms (`kits_norms`, `kits_norm_*`; cNORM and gamlss are not
installed), `kits_link_error` (a kits-specific bootstrap), `kits_dif_combine`, `kits_vertical`.
`tests/testthat/test-against-kitsPreview.R` compares the norms with cNORM through kitsPreview, but
it skips where kitsPreview is not installed.

Parametrisation: kits beta is the PCM step location (Masters, 1982), the same as dexter's beta,
psychotools' `threshpar(type = "mode")` and TAM's PCM xsi; mirt uses d_k = -cumsum(beta) with a = 1.
The TAM difference is expected: MML assumes a normal ability distribution and CML does not; it is
largest at the most extreme step.

## Metamorphic relations

Script: `kitverse/inst/validation/metamorphic.R`. 10 functions × 3 relations; 30 of 30 hold.

| Function | Relations |
|---|---|
| `kits_item` | Permuting persons or items; duplicating every person (beta unchanged, SE / √2); adding an all-NA item and person |
| `kits_subj` | beta + c gives theta + c, SE unchanged; person order; 0 and full score finite and increasing |
| `kits_itemfit` | Permutation; duplication (fit unchanged, n doubled); beta shift leaves conditional fit unchanged |
| `kits_subj_table` | beta shift; duplication; an empty raw score leaves the other rows unchanged |
| `kits_linking` | Self-link gives slope 1, intercept 0; A→B and B→A are inverse; target + c gives intercept + c |
| `kits_link_line` | Affine target; point order; points on a line give that line |
| `kits_convert` | A→B then B→A is the identity; SE × 2 gives SE × 2; NA stays NA |
| `kits_drift` | Reference shift; person order; renaming an item renames only |
| `kits_norms` | Score + c; score × 2; ages outside the range leave the table unchanged |
| `kits_norm_z` | z of the percentile score equals qnorm(p); monotone in score; out-of-range age gives NA |

## Report snapshot

kits writes no report file of its own; its outputs are in the kitsUI snapshot (two runs
byte-identical, `kitsUI/inst/validation/results/snapshot_md5.txt`).

## Mutation testing

20 mutants in copies of the source; 18 killed (90%). Scripts and results:
`kitverse/inst/validation/mutation/`.

| Error type | Mutants | Killed | Surviving (and the missing test) |
|---|---|---|---|
| Formula | 9 | 8 | K06 (`cml.R`: centring constant / (n + 1)): equivalent, the constant cancels in gamma_r |
| Interval | 4 | 4 | |
| Boundary | 6 | 5 | K19 (`drift.R`: `>= min_n` became `> min_n`): test added in 0.26.1 |
| Data management | 1 | 1 | |

## Faults and regression tests

From NEWS and git (`kitverse/inst/validation/results/bugs.csv`): 49 faults recorded. The 15 versions
before 0.25.0 are undated and their history is squashed, so when those faults were found is not
known. The table lists the faults found after release, plus the most consequential earlier one.

| Version | Fault | Error type | Caught by | After release | Regression test |
|---|---|---|---|---|---|
| 0.26.1 | `kits_vertical()` stopped ("missing value where TRUE/FALSE needed") when a cohort had one person at a shared wave (K49) | Boundary | History review | yes | `test-api22.R` |
| 0.26.1 | `kits_drift()` boundary at exactly `min_n` untested (K19) | Boundary | Mutation | not a fault | `test-link-drift-norms.R` |
| 0.26.0 | A test relied on a score gap that the demo data happened to have | Boundary (test) | Demo data change | yes | `test-lrtest.R` |
| 0.25.3 | A test differed in the 16th digit on Windows and macOS | Other (test) | r-universe | yes | yes |
| 0.16.0 | WLE stopped at a local optimum in about 5% of data sets | Numerical estimation | Cross-engine (kitsPreview) | no | yes |

## Metrics

| Metric | Value |
|---|---|
| Acceptance criteria written beforehand | 0 |
| First-generation pass rate | 2 of 3 dated versions (67%); earlier versions cannot be reconstructed |
| Faults caught per layer | Engine 0, metamorphic 0, mutation 1 missing test, history review 1 open fault; earlier: cross-engine 1, review and not stated for the rest |
| Faults found after release | 3 (strict: after 2026-11-01), 4 (broad); all fixed in 0.26.1 |

## Reproduce

```r
Sys.setenv(KIT_VALIDATION_OUT = "validation_results")
source(system.file("validation", "engine.R", package = "kits"))
source(system.file("validation", "metamorphic.R", package = "kitverse"))
```
