# kits 0.25.0

## The algorithm layer: `algo_`

* Twelve plain algorithm functions are exported: `algo_pcm()`, `algo_esf()`, `algo_cml()`,
  `algo_thresholds()`, `algo_jme()`, `algo_screen()`, `algo_mean_sigma()`, `algo_link_error()`,
  `algo_line()`, `algo_shash()`, `algo_shash_z()` and `algo_lw()`. Numbers in, numbers out: they do
  not track scales and check little. The `kits_` functions call them, so there is one implementation
  of each algorithm. They chain with the pipe, e.g.
  `kits_stats(X) |> algo_cml() |> algo_thresholds()` gives the thresholds of `kits_item(X)`.

## `kits_algorithm()` replaces `kits_code()` and `kits_equations()`

* `kits_algorithm(x)` prints, in about 30 lines, the steps behind a result (or a topic such as
  `"cml"`): one formula per step, the `kits_` and `algo_` pipelines, the key references, and where
  to read more (the tutorial file and the teaching script part).
* `kits_code()`, `kits_equations()` and the alias `kits_rtmb_code()` are removed. The teaching
  scripts stay in `system.file("scripts", package = "kits")`; the formulas (with copyable LaTeX)
  and the reproducible pipelines are in the six tutorials of `tutorial_algorithms/`.

## Person scores: the Jeffreys modal estimate

* The person scores are still called WLE. For the partial credit model the WLE equals the mode of
  the likelihood penalized by the Jeffreys prior (Firth, 1993; Magis, 2015), which kits has
  computed since 0.23.1 by maximizing log L + 1/2 log I with RTMB; the help pages now say so and
  the algorithm is `algo_jme()`. Warm's score-equation route is kept as a comparison only (part
  `"warm"` of the script `wle.R`).
* With two maxima of the penalized likelihood the higher one is taken (before: the one Newton
  reached from TAM's start value); with two equally high maxima still the smaller theta. The demo
  data have no such case.
* The SE of a person score is now 1/sqrt(I(theta)), the asymptotic SE of the WLE (Warm, 1989).
  Before it was Warm's formula sqrt(I + (I'/2I)^2) / I, which is slightly larger, most at zero
  and perfect scores. theta does not change.

## Fewer functions for the same questions

* `kits_lrtest()` is merged into `kits_dif(data, group, ...)`: `group` is a vector of groups,
  `"raw"` (the median split of Andersen's test), `"wave"` or a column name. One call gives the
  Andersen LR test (printed; the numbers in `attr(d, "lrtest")`), the Wald z of each step difference
  (`attr(d, "diff")`, `plot(d, type = "steps")`) and, with two groups, the item DIF of Bechger and
  Maris as before. The numbers are those of 0.24.1. With two groups, an item a group cannot
  calibrate is now left out of that group's DIF with a message instead of stopping.
* `kits_age_pass()` is merged into `kits_item_table(thr, norms, data = )`: with responses and ages
  in `data` the table gets the observed pass ages. Same numbers.
* `kits_drift()` is built on `kits_item()` and `kits_linking()`; same output.
* `kits_item_table(step = )`: the pass step per item, e.g. `step = c(cogD01 = 3)` for a reversed
  0-4 item whose pass is a score of 3 or more (items not named use 1). `D` is then the threshold of
  that step and the pass age counts a score of at least `step`; a `step` column says which step each
  item uses. Anchor links already took a `step` column in `anchors`.

## Norms: nlminb

* The shash norm fit (`algo_shash()`, `kits_norm_fit()`, `kits_norms()`) uses `nlminb()` with the
  RTMB gradient and Hessian, as the other fits of kits do, instead of `optim(L-BFGS-B)` with cNORM's
  settings. The model, bounds, starts and degrees are still cNORM's; on the demo data the numbers
  match 0.24.1 and cNORM to about 1e-6 (closer to cNORM than before). Convergence is judged as in
  `kits_item()`: max |gradient| of the free coefficients <= 1e-4 and their Hessian positive
  definite. A fit that stops early is restarted from where it stopped (up to 3 times), then takes
  Newton steps on the free coefficients (as `kits_item()` does) while -log L decreases; a fit that
  converges without them is unchanged. Where optim stopped before converging, the fit can now reach
  a better optimum (the by-sex example of `?kits_norms`, boys). The model gains `message` and
  `hessian_pd` and loses `grad_tol` and `restart_gain`.

## `kits_raw_equate()` removed

* The raw score of the PCM is sufficient, so `kits_subj_table()` (raw score -> WLE and SE) is the
  scoring table; with IRT ability scores no raw-to-raw relationship is needed (Kolen & Brennan,
  2014, section 6.4). A raw-score correspondence between two forms is the expected score of the
  other form at those abilities, `kits_expected_score(item, items = form_b)(tab$theta)`, and a
  subscale score on a fixed item set S is `kits_expected_score(item, items = S)(theta)`; the
  tutorials show both, with growth models and a simulation (`tutorial_algorithms/07_subscales`).

## Observed item curves: `plot(item, type = "observed")`

* `plot()` of a `kits_item()` draws, with `data =`, the observed share of children who pass each
  item against their WLE, one line per group (`group =` as in `kits_dif()`: `"wave"`, `"sex"`, a
  vector), in bins of 0.5 logits with at least 10 children; ticks at the top and bottom are the
  thetas of those who pass and fail, the dashed grey line is the model's P(X >= step) (`step =` as
  in `kits_item_table()`), and each panel gives the item's infit and outfit on `data`. It is the
  empirical curve of van Buuren and Eekhout's D-score equate groups, for the common items of KIT
  across waves. No new function; tutorial_algorithms/05_checks section 4.5.

## One wave at a time: fixed steps, new top scores

* `kits_item(fixed_params = )` fixes steps one by one: a fixed item that has a higher score in
  the new data than in `fixed_params` keeps its fixed steps, and the new steps are estimated
  (the item needs no 0 score then, but its scores above the fixed steps must be consecutive).
  Before, this stopped ("Score(s) above the highest fixed step"), so a wave where older children
  reach a score no younger child had could not be added without freeing the whole item. With it,
  each wave can be calibrated with the earlier parameters fixed, and old betas never move. The
  thresholds of an item that gains a step do move (a threshold depends on all the steps of its
  item), and so do the D-scores, pass ages and links that use them. Scoring an earlier wave again
  with the extended parameters changes its scores, so keep each wave's scores and parameters as
  they were.
* `kits_drift()`: an item with a higher score in the new data than in `ref` ("score range differs")
  is still not judged, but is now suggested for fixing (its old steps can be fixed). `ref` can also
  be a table of step parameters (`item_id`, `item_score`, `beta`; the steps of each item are 1, 2,
  ...; a table from `coef(item, rp = 0.8)` too, as beta does not depend on the RP): the thresholds
  are computed. So an item bank of `item_id`, `item_score`, `beta`, `beta_se` grows one wave at a
  time, and earlier scores are never estimated again (example in `?kits_item` and the advanced
  vignette, section 8): `kits_drift(X_new, bank)`; a drifted item gets a new id in the new wave (a
  new version, estimated freely); `kits_item(X_new, fixed_params = drift$fix)`; then the new steps
  of `coef()` are added to `bank`.

## Bug fixes

* `plot()` of the step differences with three or more groups failed when one group had no step
  among the chosen items.
* A duplicated `visit_id` in `persons` of the LR test now raises `kits_error_data`, as elsewhere.
* `kits_link_check()` reported that the same-age distributions did not converge in almost every
  run: it read the code of `nlminb()`, which often stops with "singular convergence" at the
  optimum. It now uses the gradient and the Hessian, like `kits_item()`.
* Two tests of the linking error never ran their check; they do now.

## Less code

* One `abort_if()` for the 195 one-line input checks, the helpers left over by the `algo_` layer
  removed. R code 13,804 -> 11,989 lines; every output is unchanged apart from the changes above.

# kits 0.24.1

## Shorter help pages

* Every help page keeps its title, a short description, the arguments, the value, the references it
  cites and the examples; the long details (the why, derivations, comparisons with other packages,
  version history) moved to the equations tutorial, and each page points to its section. Roxygen
  3,899 -> 2,446 lines; the package help page 329 -> 135 Rd lines.

## Teaching scripts: the other cases

* New parts in `inst/scripts/`, each checked against kits on the demo data: `cml.R` `"fixed"` (fixed
  item parameters) and `"fixed_out"` (fixed items left out of the calibration); `wle.R` `"convert"`
  (converted scores, with se_link and se_total); `link.R` `"slope"` (z of the slope against 1) and
  `"rules"` (the other screening rules); `norms.R` `"by"` (one norm per group). `kits_code(x, part = )`
  names the matching part in its header when `x` is such a case.

## Less code

* `link.R` and `lrtest.R`: the diagnostics, printing and plotting code is shorter (a shared
  `diag_flag()`); every output is unchanged (39 calls compared with `identical()`).

# kits 0.24.0

The package shows short code and formulas; the teaching moves to fixed scripts and the tutorials.

## kits_code(): an overview and teaching scripts

* `kits_code(x)` (the overview) is as before: a few lines that call kits and give exactly the numbers
  of `x`. Its pointers now name the part of the teaching script (`-> part "cml"`); `kits_norm_screen()`
  and `kits_drift()` results get an overview too.
* `kits_code(x, full = TRUE)` and `kits_code(x, part = )` return fixed teaching scripts from
  `inst/scripts/` (`cml.R`, `wle.R`, `link.R`, `norms.R`, `norm_screen.R`, `dif.R`, `drift.R`,
  `information.R`; `system.file("scripts", package = "kits")`), with a short header. Each script builds
  its inputs from the demo data, recomputes the result in base R and RTMB without kits, and ends by
  comparing with kits. Before, these scripts were generated with the numbers of `x` written in; to use
  your own data, change the inputs section. `part =` returns the header, the inputs and one section;
  every result now has parts (new: link `"points"`, screening `"screen"`, `"ages"`, drift `"points"`,
  `"line"`, `"items"`, scores `"wle"`, `"person_fit"`). The scripts no longer check the md5 of the data.
* R/code.R and R/code_overview.R: 1,892 -> 776 lines. The tests run every script on the demo data and
  compare with kits (thresholds, theta and links within 1e-6, norms within 1e-4).

## kits_equations(): short notes

* Every formula (Unicode and LaTeX) and every reference is unchanged; each note is now one phrase that
  says what the symbols are (at most about 90 characters). The longer explanations moved to the
  equations tutorial. About 27% fewer characters in the printed output.

## Removed (announced in 0.23.0)

* `kits_code(item, subj = TRUE)`: the argument `subj` of `kits_code()` for a `kits_item` is gone
  (deprecated since 0.23.0); use `kits_code(scores)` with `scores <- kits_subj(X, item)`. `kits_code()`
  of a `kits()` result keeps `subj`.
* The check of the kits 0.22.0 link names (`link$a`, `link$b`, `a_error`, ...): the `$` and `[[`
  methods of a link and the error class `kits_error_deprecated` are gone; use `intercept`, `slope` or
  `coef(link)`. No element of a link starts with a or b, so `link$a` is `NULL`.
* The warning class `kits_warning_deprecated` (nothing is deprecated now).

## Tutorials and vignettes in Quarto

* The vignettes are Quarto documents (`VignetteBuilder: quarto`; Suggests quarto). The tutorials in
  the repository are `.qmd` files rendered with `quarto render`.
* New tutorial `tutorial_models/kits_models.qmd`: each model function of kits (PCM core, CML,
  `fit_ad()`, WLE, threshold SEs, ESF, latent normal, shash, information slope, gamma SE) with its
  formula, its source read from the package at render time, what RTMB does there, and a small check.

## No derivative by hand

* The information slope I'(theta) of `kits_pcm_prob()` comes from RTMB (the information is written
  once, `information_slope()`); before, it was the sum of the third central moments. The values agree
  within 1e-12 (a test checks the identity). `pcm_item_moments()` no longer computes the fourth moment.
* The SE of the Goodman-Kruskal gamma of `kits_item_restscore()`: gamma is written as a function of
  the cell counts and RTMB gives d gamma / d n for the delta method; before, the derivative psi was
  written by hand. The SEs are unchanged (the tests compare with vcdExtra's values).

# kits 0.23.1

The model is written once and RTMB gives every derivative. No estimate of kits (calibration, WLE,
threshold SEs, latent normal of the link check, norms) rests on a derivative derived by hand any more;
the information slope I' of `kits_pcm_prob()` (sum of third central moments) and the SE of the
Goodman-Kruskal gamma of `kits_item_restscore()` (a published delta-method formula) still do. The calibration (CML), the links and the norms are unchanged to the last
bit; WLE thetas and SEs move by at most about 1e-12 (except with two maxima, below); threshold SEs by
at most about 1e-10.

## Every derivative from RTMB

* WLE (`kits_subj()`, `kits_subj_table()`, the start values of `kits_link_check()`): kits writes only
  the weighted likelihood f = log L + 1/2 log I, with log L = r theta - sum_i log S_i(theta) and the
  PCM test information I = sum_i Var(X_i | theta) from the PCM core. RTMB's Newton solver
  (`tape$newton()`) solves all item pattern x raw score combinations at once, and Warm's SE takes
  I' = dI/dtheta from RTMB. The closed-form f', f'' (third and fourth cumulants) and kits' own Newton
  steps are gone. On the demo data theta and se differ from 0.23.0 by at most 6e-13; on 1,230 and
  2,965 combinations of two simulated designs by at most 3e-13.
* WLE with two maxima (few items far apart): as TAM, the WLE is the maximum that the Newton steps reach
  from the start value, which is now TAM's, the mean beta of the items plus log((r + 0.5) / (max - r +
  1)). kits 0.16.0 to 0.23.0 searched a grid for the higher maximum (`wle_global()`); that search is
  gone. Two equally high maxima still give the lower theta (as dexter), because the start value leans
  to that side. A solution that is not a maximum of f (the Newton steps can stop at the minimum between
  two maxima when they start exactly there) is solved again from other start values.
* Threshold SEs (`kits_item()`, `coef(item, rp = )`): the delta method takes J = d thr / d beta from
  RTMB. The thresholds of an item solve sum_k (P(X >= k | theta_k) - RP)^2 = min; RTMB's Newton solver
  makes this the map beta -> thr, and its jacobian is J. Before, the partial derivatives of the implicit
  function were written by hand and computed by central differences. The thresholds themselves are
  still found by `uniroot()` (no derivative) and are unchanged; the SEs move by at most 5e-11.
* `kits_link_check()`: the latent normal of the same-age groups is fitted by RTMB and nlminb (before:
  `optim()` BFGS with numerical gradient and Hessian), and the SE of its SD comes from `sdreport()`.
  `se_method` takes the derivatives of the same-age line from RTMB. `optim()` used to stop with a
  gradient of about 1e-3, so the means and SDs move by up to 7e-6 and the same-age line values by up to
  2e-4 on the demo data; the new values have a lower -log L.
* One optimizer for every RTMB model, `fit_ad()`: nlminb on the RTMB gradient, then Newton steps with
  the RTMB Hessian (the CML and the latent normal). The CML is unchanged to the last bit.
* RTMB in the package: `ADoverload()` puts RTMB values into numeric vectors and matrices, in place of
  the `ad_zeros()` and `+ 0 * beta[1]` workarounds. The PCM core (`pcm_kernel()`, `pcm_item_moments()`)
  also takes RTMB values (log S by `logspace_add()`); numeric input takes the old path, bit for bit.

## Speed

The WLE is now solved by RTMB, which tapes every operation, so it is slower. One thread:

| Design | Combinations | 0.23.0 | 0.23.1 |
|---|---:|---:|---:|
| KIT-like: 46,900 person-waves, 24 items, 7 waves | 1,230 | 2.9 s | 3.5 s |
| 40 items, N = 3000, 10% missing at random | 2,965 | 0.9 s | 9.0 s |

The time is in building the Newton tape (about 2 ms per combination); the solve takes 0.03 s.

## kits_code()

* The scripts write only the models: one `fit_ad()` (the optimizer above) shared by every script, the
  WLE as f with RTMB's Newton solver, and the threshold SEs from RTMB. The `zeros()` helper is gone.
  Functions that put RTMB values into numeric vectors call `ADoverload()`, as the package does: in a
  session where packages such as dplyr are attached, the overloads that RTMB attaches by itself do not
  always take effect. The full script of `kits()` is about 270 lines (was 300).
* Fixtures `tests/testthat/fixtures/code0231/`.

## kits_equations()

* The WLE, threshold and same-age texts say what RTMB computes; the WLE no longer lists I'' and the
  curvature h.

## Removed

* `kits_inter()` (the interaction model) is gone, with its methods, its `kits_code()` scripts and its
  `kits_equations("inter")` text. For the item fit use the conditional infit and outfit of
  `kits_item()` (`item$itemfit`, `kits_itemfit()`). The interaction model of dexter stays in
  `kitsPreview::kits_inter()`. The shared form note of `kits_equations("subj_table")` loses its
  sentence about the interaction model, and the error of `kits_code()` for other objects no longer
  lists `kits_inter()`.

## Less duplicated code

* `cml_optimize()` and `ad_optimize()` (two copies of nlminb plus Newton steps) are one `fit_ad()`.
* `ad_zeros()`, `ad_pos()`, `ad_conv()`, `num_conv()` and the unused `poly_item()` are gone;
  `ad_clamp()` moved next to the shash norm that uses it.
* The WLE's own copy of the PCM probabilities (`wle_terms()`, `wle_item_curve()`) is gone: the WLE
  uses the PCM core.
* `copied_thresholds()` (fixed items left out) uses `step_thresholds()` and lives next to it.
* The unused `link_args_text()` is gone.
* One internal helper per job (results unchanged to the last bit on the demo data):
  * `person_sets()` and `number_booklets()` (utils.R): the item set and raw score of each person
    (7 copies) and the numbering of the response combinations (3 copies). The start values of
    `kits_link_check()` use `add_wle()`.
  * `poly_times()` (PCM core): the polynomial products of the CML and of the Lord-Wingersky recursion
    of `kits_precision()`. `betas_by_item()` moved to the PCM core and is used by `kits_esf()`,
    `kits_item_restscore()` and `simulate()`.
  * `kits_pcm_prob()` and `kits_esf()` share their setup, long table, print header and
    `as.data.frame()`; the category curves of `plot(item)` and `plot(pcm_prob)` share one helper;
    `print(kits_stats)` uses `print_rows()`.
  * Checks and messages: `some_ids()`, `has_gap()`, `can_calibrate()`, `need_params()`,
    `need_number()` (about 16 single-number arguments), `need_same_scale()` and `need_d_age()`.
    Messages about items without parameters add ", ..." after 10 items, and `kits_norm_data()` says
    the same as `kits_norm_z()` when `d` and `age_m` do not match.
  * Line arithmetic lives in line.R (`line_chain()`, `line_minus()`, `line_stretch()`). With three or
    more links, the overall slope of `kits_equations()` and `kits_precision()` is now composed link by
    link, as in `summary()`, so its last digit can move (1 ulp).
* The R code is 16,631 lines (0.23.0: 17,648); `kits_inter()` accounts for about 810 of the 1,017.

# kits 0.23.0

## Item objects no longer carry a person-wave id

* `item$design` of a `kits_item` had row names made from the first person-wave of each item pattern
  (for example `M03_011e336_M031`), so a shared item object carried one visit id per pattern.
  The row names are now 1, 2, ...; nothing else changes. A test checks that no visit or child id
  is anywhere in an item object. This was also the case in 0.22.0 and earlier.

## The a and b of a link are renamed `intercept` and `slope`

**Read this first.** A link is the line new value = intercept + slope x old value. Its elements
are now called `intercept` and `slope`; kits 0.22.0 called them `a` (the intercept) and `b` (the
slope). The letters a and b mean different things in different places (Kolen and Brennan, TAM,
plink and equateIRT write A or a for the slope and B or b for the intercept), so kits no longer
uses them. **No number changes.**

| Source | Intercept | Slope |
|---|---|---|
| kits 0.22.0 and before, the KIT link table `02_link_table.csv` | `a` | `b` |
| Kolen & Brennan (2014), TAM, plink, equateIRT | B (or b) | A (or a) |
| kits 0.23.0 | `intercept` | `slope` |

| kits 0.22.0 | kits 0.23.0 |
|---|---|
| `link$a`, `link$b` | `link$intercept`, `link$slope`, or `coef(link)` (`c(intercept = , slope = )`, the order of `lm()`) |
| `link$a_error`, `link$b_error` | `link$intercept_error`, `link$slope_error` |
| `link$a0`, `link$b0` | `link$intercept0`, `link$slope0` |
| `link$ab_cov` | `link$vcov`, a 2 x 2 matrix with dimnames `intercept`, `slope` (`ab_cov` is its off-diagonal) |
| `link$draws[, "a"]`, `[, "b"]` | `link$draws[, "intercept"]`, `[, "slope"]` |
| `new value = a + b x old value` | `new value = intercept + slope x old value` (print, `kits_equations()`, `kits_code()`) |

* `link$a`, `link$b` (and `a_error`, `b_error`, `a0`, `b0`, `ab_cov`, also with `[[`) stop with an
  error of class `kits_error_deprecated` that names the new element, instead of returning NULL by
  partial matching. No element of a link starts with "a" or "b". This check will be removed in
  kits 0.24.0.
* The same names in `as.data.frame()` of links (`intercept`, `slope`, `intercept_error`,
  `slope_error`, `cov_intercept_slope`, `intercept0`, `slope0`) and of drift checks, `confint(link)`
  (rows `intercept`, `slope`), `kits_drift()` (`intercept`, `slope`), `kits_link_check()`
  (`methods`, `lines`), the `domains` table of a link, the conversions recorded by `kits_convert()`
  (`kits_info()`), the diagnostics ("slope vs 1" is (slope - 1) / slope_error; "max |domain
  intercept - intercept|") and the scripts of `kits_code()` (`line[["intercept"]]`, `line[["slope"]]`,
  `draws[, "intercept"]`).
* `kits_equations()` writes new = intercept + slope x old; where a formula needs letters it uses B
  (intercept) and A (slope), as Kolen and Brennan (2014), and says so once.
* Every converted score, threshold and `se_link`, and every intercept, slope, linking error,
  covariance and draw is identical to kits 0.22.0 (tolerance 0, tested against numbers saved by
  0.22.0); the new scripts of `kits_code()` give the same intercept, slope and draws as the scripts of
  0.22.0, bit for bit.

## Old call -> new call

| Old call (0.22.0) | New call (0.23.0) | Status |
|---|---|---|
| `link$a`, `link$b`, `a_error`, `b_error`, `a0`, `b0`, `ab_cov` | `intercept`, `slope`, `intercept_error`, `slope_error`, `intercept0`, `slope0`, `vcov`; `coef(link)` | renamed; the old names stop with a clear error (`kits_error_deprecated`) until 0.24.0 |
| `person_id` | `visit_id` (one child at one visit, a person-wave) | renamed; `person_id` in the data still read for one version, with a message (`kits_message_data`) |
| `pid` | `subj_id` (the child) | renamed; `pid` still read for one version, with a message |
| `booklet_id`, `booklet_score`, `booklet_max` | `pattern_id`, `raw`, `raw_max` | renamed (scores, raw-score tables, `item$design`, `kits_stats()`) |
| "booklet" in prints, messages and help | "item pattern" (which items a person-wave answered); "standard form" stays | renamed |
| `item$stats$booklets`, `n_booklets` of `as.data.frame(item)` | `item$stats$patterns`, `n_patterns` | renamed |
| `kits_booklet_table` (class of `kits_subj_table(item)`) | `kits_pattern_table` | renamed |
| `kits_information(item, booklet_id = )`, `kits_expected_score(item, booklet_id = )` | `pattern_id =` | renamed argument |
| `kits_code(item, subj = TRUE)` | `kits_code(scores)`, `scores <- kits_subj(X, item)` | deprecated, removed in 0.24.0 |
| — | `kits_code(scores)`: the WLE and person fit from the item parameters the scores keep, then the conversions | new method |
| — | `kits_code(drift)`, `kits_equations(drift)`, `kits_equations("drift")` | new methods |
| `kits_age50()`, `kits_mean_sigma()`, `kits_link_se()`, `kits_link_screening()`, `kits_wle_table()`, `kits_plot_cohort()`, `kits_plot_norms()`, `kits_item_positions()`, `kits_ladder_marks()`, `kits_plot_growth()`, `kits_plot_ladder()` | see the 0.22.0 table below | removed |
| `kits_subj_table(item, data, wave_m)`, `plot(norm, d = , age_m = )` | `kits_subj(data, item, persons) \|> kits_subj_table()`, `plot(norm, scores = )` | removed (stop with an error naming the argument) |
| `s$z <- kits_norm_z(s, norm)` | `s <- s \|> kits_norm_z(norm)` | removed (stops with an error that shows the new form) |
| `dataSrc =`, `D =` (kits 0.17.0 names) | `data =`, `d =` | removed |
| `print(kits(...))` ended with the next steps of a calibration | the next steps after scoring (`plot(x)`, `x$subj \|> kits_subj_table()`, `kits_code(x)`, ...) | fixed |
| `kits_itemfit(X, item)` on the calibration data | `item$itemfit` (the same table), `plot(item, type = "fit")`; `kits_itemfit()` stays for other data | new element |

Exported functions: 55 in 0.22.0, 48 now (the eleven deprecated functions are gone; `kits_esf()`,
`kits_item_restscore()`, `kits_lrtest()` and `kits_pcm_prob()` are new).

## Item fit travels with the item object, and is conditional

* `kits_item()` now computes item fit (`x$itemfit`), the table of `kits_itemfit(data, x)` on the
  calibration data: the same code and identical numbers, columns, class and rows; only its
  `kits_info()` records the call of the calibration. Diagnostics live with the object they describe:
  person fit is in the scores (`kits_subj()`), item fit in the items. The table holds one row per
  item (`n`, `pvalue`, `rir`, `infit`, `outfit`) and nothing per person-wave, so the item object
  can still be shared. The same holds with `fixed_params =` and in `kits()` (`x$item$itemfit`).
* **Conditional infit and outfit.** `kits_itemfit(data, item, type = c("conditional", "wle"))`.
  The new default, `"conditional"`, takes the expected item score and its variance given the raw
  score of the person-wave within his or her item pattern, from the item parameters alone (the
  elementary symmetric functions of the CML; Christensen, Kreiner & Mesbah, 2013; Müller, 2020),
  with no ability estimate. Responses with a conditional variance of 0 (a zero or a perfect score,
  or an item pattern of one item) are left out. On complete data the numbers equal
  `iarm::out_infit()` (to 1e-15 with the same item parameters; to eRm's convergence, about 1e-5,
  with kits' own CML); kits does not depend on iarm or eRm (the reference values are saved in the
  tests). `type = "wle"` is the item fit of kits 0.22.0 (each person at his or her WLE, as
  Winsteps, TAM and eRm do), unchanged. The WLE version is biased downward on short forms, even
  when the PCM holds (Kreiner & Christensen, 2011): on `kitm03` lan its infit is 0.68 to 0.94
  (mean 0.86), the conditional infit 0.91 to 1.06 (mean 0.99). `x$itemfit` of `kits_item()` and
  `kits()` is the conditional version, so its numbers differ from those of the WLE version; the
  type is recorded in `attr(x, "type")` and `kits_info()`, printed in the first line of a
  `kits_itemfit` and in the title of its `plot()`.
* The printout of a `kits_item` and of a `kits()` has one more line before `Next:`, for example
  `Item fit (conditional): infit and outfit within 0.6-1.4 for all 24 items.` or
  `Item fit (conditional): 8 of 74 items outside 0.6-1.4 (6 underfit, 2 overfit: ...); ...`
  The range is the default `range` of `plot()` of a `kits_itemfit`.
* A new diagnostic "item fit" (`summary(item)`, `kits_check(item)`): a note naming the items
  outside the range, `ok` otherwise. It is a note, never a warning, because real data have some
  misfit items; so the `Checks:` line of an item can now read "1 note" where it read "OK".
* When item fit cannot be computed (for example when every person-wave has a zero or a perfect
  score), `x$itemfit` is NULL and the diagnostic says why; the calibration never
  stops for it.
* `plot(item, type = "fit", items = )` draws the plot of `plot(kits_itemfit(...))`; the default
  `type` of `plot(item)` is unchanged.
* A new column `flag` in the item-fit table (`kits_itemfit()`, `x$itemfit`, `x$item$itemfit`),
  after `outfit`, in the usual Rasch terms (Linacre; Wright & Linacre, 1994), against the same
  range 0.6 to 1.4: `"underfit"` when infit or outfit is above 1.4 (noise; it wins when both
  apply), otherwise `"overfit"` when one is below 0.6 (too predictable, usually harmless), else
  `"ok"` (the limits themselves are ok); `NA` when both are NA. The printout and the diagnostic
  count the two kinds, underfit items first; `plot()` colours the points by flag (underfit
  orange, overfit blue, ok gray). A flag is a prompt to look at the item, not a rule to remove it.
* Use `kits_itemfit()` to check the same items on other data, for example a new wave with fixed
  parameters.
* New `kits_item_restscore(data, item, min_n = 100)`: Kreiner's (2011) item-restscore test. For
  each item it compares the observed Goodman-Kruskal gamma between the item score and the rest
  score with the gamma the PCM implies given the CML item parameters (the conditional
  distribution of the item score given the total score; no theta needed), with the asymptotic SE
  of the observed gamma, z, p and Benjamini-Hochberg adjusted p over all rows. On complete data
  the numbers equal `iarm::item_restscore()` (observed gamma and SE exactly; expected gamma to
  eRm's convergence, about 1e-6; reference values saved in the tests). In incomplete designs the
  test runs within each item pattern with at least `min_n` persons (rows are item x pattern),
  because rest scores of different patterns are not comparable; smaller patterns are skipped
  with a note. `print()` lists the significant rows and the largest |observed - expected|;
  `plot()` draws observed against expected gamma with the identity line. With large samples
  small differences become significant: look at the size of the difference.
* New `kits_lrtest(data, item, split = "raw")`: the Andersen (1973) likelihood-ratio test of the PCM
  by CML. Person-waves are split by raw score (at the median within each item pattern, ties to the
  low group), by wave (`split = "wave"`), by a column of a `kits_resp()` or `persons` (for example
  `"sex"`), or by a vector; each group is calibrated on its own and LR = 2 (sum of the groups'
  conditional log-likelihoods - the joint one), with df = the groups' free step parameters minus the
  joint ones. Each group uses the steps it observes (a step nobody in the group reached, or everybody
  reached, is not estimated there and not counted in df; `x$unobserved` lists them), so a split by
  wave works on KIT's age-based design, where `eRm::LRtest()` drops every item;
  `unobserved = "drop"` follows eRm. With large N the p value is almost always tiny, so the result
  also gives each group's step parameters on a common origin, every step difference between groups
  with its Wald z, and prints the largest |difference| in logits; `plot()` draws the graphical model
  check (group against group, or each group against the joint calibration, with +/- 1.96 SE). On
  complete data the LR, df and p equal those of `eRm::LRtest()` (to 3e-8; tested against saved eRm
  1.0.10 values; kits does not depend on eRm). The groups are calibrated with the same CML code as
  `kits_item()`.
* New `kits_esf(item, items = NULL, pattern_id = NULL)`: a teaching view of the elementary
  symmetric functions gamma_r of CML for one item pattern, with the conditional distributions they
  imply, P(X_i = x | r) = e_ix gamma^(-i)_{r-x} / gamma_r, and the conditional expected score
  and variance of each item given the raw score (`$gamma`, `$prob`, `$moments`). It is not a
  second implementation: the gammas come from the same internal code as the CML likelihood
  (`cml_gamma_matrix()`, with the same centring), and a test checks that rebuilding -log L from
  `kits_esf()` gives the calibration's minimum exactly. The conditional infit/outfit of
  `kits_itemfit()` and the expected gamma of `kits_item_restscore()` use the same internal core
  (`esf_moments()`). `plot()` draws the conditional expected-score curves E(X_i | r);
  `as.data.frame()` gives the gamma table.
* `kits_equations("fit")`: the conditional residuals, the residuals at the WLE and the item flag;
  the item topic says that `x$itemfit` uses them.
* No other number changes: calibrations, scores, links and norms are identical to kits 0.23.0
  before this entry. `kits_item()` takes a little longer (under a second for `kitm03` lan).

## Shared cores and teaching functions

A core computation now lives in one internal function, and every part of kits calls it; the
exported teaching functions show that same core.

* **One PCM core.** The category probabilities P(X_i = x | theta) and their moments had six copies
  (thresholds, `kits_sim()`, item and person fit, the test characteristic curve, `kits_precision()`,
  the plots). They now come from one internal core (`pcm_probs()`, `pcm_item_moments()`): eta_x =
  sum_{g <= x} (theta - beta_ig), the row maximum subtracted before `exp()`, and **central** moments
  (Var = sum_x (x - E)^2 P, not E[x^2] - E^2, which cancels at extreme theta). The thresholds,
  `kits_sim()`, `kits_information()`, `kits_expected_score()`, `kits_raw_equate()`,
  `kits_precision()`, the person fit of `kits_subj()`, `kits_itemfit(type = "wle")` and the item
  plots use it. The WLE solver (`wle_terms()`, `wle_item_curve()`) keeps its own vectorised copy
  for now.
* **Numbers.** Calibrations (beta, thresholds and their SE), scores (theta, SE), links,
  conversions, norms and `kits_sim()` data are bit-identical to 0.23.0 before this change. The
  functions that moved from raw to central moments change in the last digits (kitm03, lan and
  motor): `kits_information()`, `kits_expected_score()` and raw-score equating by at most 6e-14 (on
  values up to 100), `kits_precision()` by 5e-14, the category probabilities by 1e-15. The WLE item
  fit and the person outfit change by at most 2e-12 and 1.4e-10 (relative 3e-12). That is where an
  item is far from a person and its variance W is about 1e-7: E[x^2] - E^2 had lost about 8 digits
  there, and the central moment has not (for a 0/1 item the old formula had a relative error up to
  7e-4 at |theta - beta| = 30; the new one has 5e-16). The conditional item fit (`item$itemfit`,
  the default of `kits_itemfit()`) does not change.
* New `kits_pcm_prob(item, theta, items = NULL)`: a teaching view of the core. `$prob` has
  P(X_i = x | theta) (long: item_id, theta, x, prob), `$moments` the expected item score, its
  variance (the item information) and the third central moment, and `$test` the expected raw
  score, the test information and its slope I'(theta), identical to
  `kits_expected_score(item)(theta)` and `kits_information(item)(theta)`. `print()` contrasts it
  with `kits_esf()`, which gives P(X_i = x | r), conditional on the raw score and free of theta.
  `plot()` draws the category curves of one item, or `type = "expected"` and `"information"`.
* `summary()` of converted scores (a `kits_subj()` result put through `kits_convert()`) shows the
  conversion: the links in the order applied (from, to, intercept, slope), the line they make
  together (new value = intercept + slope x old value; chained links compose as lines, slope =
  slope_2 x slope_1 and intercept = intercept_2 + slope_2 x intercept_1), and, when every link has
  bootstrap draws, the linking error at the mean theta, the same computation as the `se_link` of
  `kits_convert()`. The numbers are in `summary(s)$conversion` (`chain`, `intercept`, `slope`,
  `theta_mean`, `value_mean`, `se_link`). Scores without a conversion print as before.
  `kits_equations()` of a link has a new section on chains of links.
* **One line algebra.** `line_apply()`, `line_resid()`, `line_compose()`, `line_invert()` and
  `mean_sigma_moments()` are the only place where lines are applied, composed and inverted: the
  mean/sigma line, the screening residuals, `kits_convert()`, `se_link`, `summary(link, theta = )`,
  `kits_link_check()` (same-age and concurrent lines) and `kits_drift()` use them. No number
  changes.
* `kits_vertical()`: the linking error of a mean is computed from the links themselves (`se_link`
  at the mean theta) when the scores record them (scores of `kits_subj()` converted with
  `kits_convert()`); otherwise it is approximated, as before, from each child's `se_link` (a
  quadratic fitted by least squares). The two agree up to rounding (at most 1.8e-15).
* shash norms: the bounds (log sigma and epsilon in -10..10, log delta in -2..2, as in cNORM) are
  one constant, and one helper maps the coefficients to mu, sigma, epsilon and delta for both the
  fit and the norm table. Norms are bit-identical.
* `step_reach()` counts, per item step, the persons below and at or above it; the CML existence
  check, the CML statistic n_ig and `kits_drift()` use it. Each keeps its own rule; messages and
  results are unchanged. (`kits_lrtest()` keeps its own rule, on the score range of each item in
  a group, which counts cannot replace.)
* `kits_link_check()`: the latent normal of the same-age check takes log S_i(theta) from the PCM
  core. The log-likelihood changes by about 1e-16, and `optim(BFGS, reltol = 1e-12)` then stops at
  a point that differs in the 8th digit: the same-age lines move by at most 2.3e-8 (kitm03 and
  kitm36) and 6e-13 (a simulation). The concurrent-CML line and the link itself do not change.
  So `kits_link_check()` same-age lines can differ from 0.22.0 in the 8th digit (shared PCM core).

## Bug fixes

* `kits_item_table()` no longer mixes response probabilities. With a `kits_item()` as `thr`, the
  thresholds now come from `coef(item, rp = rp)`; before, they were always RP50, so `rp = 0.8`
  silently paired RP50 D positions with `age80`. A thresholds table whose `attr(, "rp")` (RP50
  when not recorded) differs from `rp` now stops with an error of class `kits_error_input`.
  At the default `rp = 0.5` the table is unchanged.
* The DIF direction labels are now decided in one place (internal `dif_direction()`), used by
  `kits_dif_combine()`, `kits_item_table()` and `plot.kits_dif()`. The labels name the actual
  groups of the `kits_dif()` result (`attr(, "groups")`): with `ref = "girl"` a positive DIF now
  reads "favors boys" (before, `kits_item_table()` and `kits_dif_combine()` always called it
  "favors girls"). `kits_dif_combine()` stops (class `kits_error_input`) when the two results use
  different `ref`. Tables without groups are read as `ref = "boy"`, as before.
* `plot()` of a `kits_item_table` colours the DIF by its side (favors the focal group: orange,
  upward triangles; favors the reference group: blue, downward triangles), so groups other than
  boys and girls (for example `ref = "M"`, "favors F") no longer stop with "Unknown DIF value";
  "favors girls" and "favors boys" keep their colours. `kits_equations("dif")` writes the
  direction for any groups (focal if dif > cut, reference if dif < -cut).
* `kits_item_table()` gains `dif_cut` (default 0.3): the cut for one `kits_dif()` given as `dif`;
  it was fixed at 0.3. A `kits_dif_combine()` result keeps the `DIF` decided by its own `cut`.
  With `ref = "boy"` and cut 0.3 all outputs are identical to before.

## Details

* `kits_code()`: the scripts call the Thurstonian thresholds `thr` and their SE `thr_se` (they were
  `tau` and `tau_se`), as `kits_equations()`; `tau` is kept for the GSED difficulty of an anchor.
  Only the names change; the numbers are identical.
* `kits_code(scores)`: the scores of `kits_subj()` keep the item parameters (summary numbers, no
  child's data) since 0.22.0; the overview scores `X` with `kits_subj(X, parms)` (identical numbers)
  and applies the conversions recorded by `kits_convert()` line by line (`intercept + slope * theta`,
  `slope * se`);
  `full = TRUE` writes the WLE in RTMB (steps 1-4), person fit and flag (5-7) and the conversions (8),
  standalone. Scores of two calibrations stacked with `rbind()` stop: write the code of each first.
* `kits_code(drift)`: the common steps (free thresholds of the new data, thresholds of `ref`, persons
  below and above each step) are written into the script, which redoes the judged steps, the
  mean/sigma line with one round of screening, the drifted items and the items suggested for
  fixing. `kits_equations(drift)` writes the rules and the line of this drift check.
* The id columns: `kits_resp()`, the data of every function, `persons =` (also of `kits()` and
  `kits_link_check()`) read `person_id` and `pid` as `visit_id` and `subj_id` for one version, with
  a message. The demo data `kitm03` and `kitm36` have the new names (values unchanged).
* Item patterns keep their ids (`bk1`, `bk2`, ... in order of appearance, as before).
* `kits_code(full = TRUE)` changed with the names, so the reference scripts of the tests were
  written again (`tests/testthat/fixtures/code023/`); the link scripts of 0.22.0 are kept in
  `fixtures/code022/`, and a test runs both and checks that the new `intercept`, `slope` and draws
  are the old `a`, `b` and draws, bit for bit (the scripts of 0.16.0 stay in `fixtures/code016/`).
* Internals: no warning of class `kits_warning_deprecated` comes from normal use; the condition class
  `kits_warning_separation` (only from `kits_age50()`) is gone; the new error class
  `kits_error_deprecated` marks the old names of a link.

# kits 0.22.0

An API cleanup: one entry point for each task, so the package is simpler to learn.
All numbers are unchanged from 0.21.0 (calibrations, links, scores, tables, norms, precision, equating and
DIF are identical on a set of simulated objects computed with both versions).
The old functions and forms still work in this version, give the same results as in 0.21.0, and warn with
class `kits_warning_deprecated` naming the replacement. They will be removed in the next version (0.23.0).
See `?kits_deprecated`.

## Old call -> new call

| Old call (0.21.0) | New call (0.22.0) | Status |
|---|---|---|
| `kits_age50(pass, age_m, wave_m)` | `kits_age_pass(pass, age_m, wave_m)$p50` | deprecated, removed in 0.23.0 |
| `kits_mean_sigma(from, to, id, item)` | `points \|> kits_link_screen() \|> kits_link_line() \|> kits_link_error()`, with `points = data.frame(id, item, from, to)` | deprecated, removed in 0.23.0 |
| `kits_mean_sigma(..., from_scale, to_scale)` | `kits_link_screen(points, from_scale = , to_scale = )` (new arguments) | deprecated, removed in 0.23.0 |
| `kits_mean_sigma(..., spread = "se", se_from, se_to)` | columns `se_from`, `se_to` of `points`, `kits_link_screen(spread = "se")` | deprecated, removed in 0.23.0 |
| `kits_link_se(co, d_link, theta = x)` | `tibble(theta = x) \|> kits_convert(co, d_link, from = "lan M36")` (its `theta` is the old `value`) | deprecated, removed in 0.23.0 |
| `kits_link_screening(link, theta)` | `summary(link, theta = )$screening` (and printed by `summary()`) | deprecated, removed in 0.23.0 |
| `kits_wle_table(item)` | `kits_subj_table(item)` | deprecated, removed in 0.23.0 |
| `kits_subj_table(data, item, wave_m)`, `kits_subj_table(item, data, wave_m)` | `kits_subj(data, item, persons) \|> kits_subj_table()` (from the scores) | form deprecated, removed in 0.23.0 |
| `kits_plot_cohort(waves, scores)` | `plot(kits_vertical(s), type = "median")` | deprecated, removed in 0.23.0 |
| `kits_plot_norms(norms, scores)` | `kits_norms(s, by = "sex") \|> plot(scores = s)` (several domains: `by = c("domain", "sex")`) | deprecated, removed in 0.23.0 |
| `kits_item_positions(items, norm)`, `kits_ladder_marks(norm)` | `kits_item_table(thr, norms, age_pass = , dif = , stems = )` | deprecated, removed in 0.23.0 |
| `kits_plot_growth(norm, items)` | `plot(kits_item_table(...))` (`type = "growth"`) | deprecated, removed in 0.23.0 |
| `kits_plot_ladder(items, marks)` | `plot(kits_item_table(...), type = "ladder", ages = )` | deprecated, removed in 0.23.0 |
| `s$z <- kits_norm_z(s, norm)` | `s <- s \|> kits_norm_z(norm)`: the scores with new columns `z` and `pct` (numbers in, numbers out as before) | form deprecated, removed in 0.23.0 |
| `plot(norm, d = s$theta, age_m = s$age_m)` | `plot(norm, scores = s)`; with a norm by group each panel draws the children of its group | form deprecated, removed in 0.23.0 |
| `ages <- kits_norm_ages(s)`, `kits_norms(s, ages, pad = 0)` | `kits_norms(s)`: without `ages` the rule of `kits_norm_ages()` (and `pad = 0`); printed | `ages` optional |
| `kits_subj(X, item, persons = persons)`, `kits_age_pass(pass, age_m, wave_m)`, `kits_dif(X, sex)` | the same, or first `resp <- kits_resp(long)`, then `kits_item(resp)`, `kits_subj(resp, item)`, `kits_age_pass(resp)`, `kits_dif(resp, by = "sex")`, `kits_itemfit(resp, item)`, `kits_inter(resp)` | new function |
| a matrix without row names, a data.frame without `person_id`: error | the rows are numbered 1, 2, 3, ... with a message | new behaviour |
| `kits_subj(..., persons)` keeps `wave_m`, `age_m`, `sex`, `cohort` | also `pid` (the child id) | new column |
| `plot(kits(X))` (Wright map) | the same, or `plot(item, type = "wright", scores = s)` | new view |
| `confint(item)`, `vcov(item)`, ... (none) | `confint()` of a `kits_item`, link and `kits_dif`; `vcov()`, `simulate()` of a `kits_item`; `logLik()`, `nobs()` (so `AIC()`, `BIC()`) and `update()` of a `kits_item` and a `kits_norms` | new methods |
| `dplyr::bind_rows()` of kits tables | `rbind()` or `purrr::list_rbind()`, which keep the records of every table; `bind_rows()` keeps the first table's and warns | new warning |
| `kits_anchors(...)` | the same, with a new `domain` column (lan, motor, cog, from the GSED ids) | new column |
| the domain of a table in a `score` column (the stacked norms of `kits_plot_norms()`, the waves of `kits_plot_cohort()`, the items and norm of `kits_item_positions()`, `kits_plot_growth()`, `kits_plot_ladder()`, `kits_ladder_marks()`) | a `domain` column: `kits_norms(by = c("domain", "sex"))`, `list_rbind(names_to = "domain")` of `kits_vertical()` results, the `domain` column of `kits_item_table()` and of the norms it takes; `score` is left for a score (as in `item_score`, `booklet_score`) | renamed; the deprecated functions still take and give `score` |
| `kits_linking(items, anchors)` with a domain without anchors | the same line (as before); now printed (`same line, no anchors: social`) and recorded in `kits_info()` and `kits_code()` | new note |
| `plot(kits_vertical(s))` (mean +/- SD) | the same, or `plot(v, type = "median")` (median, P25 to P75); one panel per domain; `cohorts`, `age_ticks`, `y_range`; a note when the median age differs from the wave | new options |
| fit one norm per sex yourself and stack | `kits_norms(s, ages, by = "sex")`, or `by = data.frame(sex = sex)` with scores as numbers | new argument |
| `kits_norm_z(s, norm)` | with a norm by group each child gets the curves of its group (`by =` with numbers) | new argument |
| `plot(norm)` | the same; a norm by group draws one panel per domain x group | new view |
| — | `plot(item, type = c("category", "expected"), items = )`: category probability curves | new method |
| — | `plot(dif, items = , cut = 0.3)`: DIF per item with its 95% CI | new method |
| — | `plot(fit, type = c("both", "infit", "outfit"), items = , range = c(0.6, 1.4))`: item fit | new method |
| — | `plot(kits_age_pass(...), items = )`: ages at which children pass, easy to hard | new method |
| — | `plot(drift, type = c("scatter", "difference"), items = )`: drift | new method |
| — | `plot(screen)`: share at the lowest score per wave | new method |
| `kits_vertical(s)` | the same, with new columns `P25` and `P75` | new columns |
| `kits_information(item)`, `kits_expected_score(item)` | the same; now S3 generics with methods for `kits_item`, `kits`, `kits_fit` and a table of parameters | unchanged interface |
| `?kits_stats`, `?kits_cml`, `?kits_thresholds` | one page, `?kits_item_steps` | help merged |
| `?kits_link_points`, `?kits_link_screen`, `?kits_link_line`, `?kits_link_error` | one page, `?kits_link_steps` | help merged |
| `?kits_norm_data`, `?kits_norm_fit`, `?kits_norm_table` | one page, `?kits_norm_steps` | help merged |

Exported functions: 53 in 0.21.0, 55 now (new `kits_resp()` and `kits_item_table()`). The eleven deprecated
functions are still exported in this version, so 44 remain after 0.23.0. The new plots and base methods are S3
methods (`plot()`, `print()`, `summary()`, `confint()`, ...), not new functions.
Every `plot()` method chooses a view with `type =` (the first value is the default) and the items with
`items =` (item ids).

## Details

* `kits_subj_table()` is the one raw-score table. Given item parameters it lists the WLE of every raw score of
  every booklet (what `kits_wle_table(item)` gave); given scores (`kits_subj()`), the standard form of each wave
  with the children of the wave (see Round 3). The forms with responses (`data`, `wave_m`, in either order) still
  work in this version and warn.
  `kits_wle_table(item, design)` with your own booklets: give the parameters of the items of the booklet,
  `kits_subj_table(dplyr::filter(parms, item_id %in% items))`.
* `kits_link_screen()` gains `from_scale` and `to_scale`, for the scales of your own points; a
  `kits_link_points()` result records its scales already, and giving them there stops with an error.
  The documentation of the screening rules, the linking error and the elements of a link moved from
  `?kits_mean_sigma` to `?kits_link_steps` and `?kits_linking`.
* `summary(link, theta = )` adds the table `screening` (the line before and after screening at `theta`, the shift,
  and with bootstrap draws `se_link` and `z`) and prints it. `summary()` was chosen over an option of
  `as.data.frame(link)`, which stays one row per link so that links can be stacked.
* The documentation of the linking error at chosen theta values (one link and a chain) moved to `?kits_convert`.
* `kits_vertical()`: `P25` and `P75` (`stats::quantile()`, type 7) follow `median`. The results of several domains
  stack with `list(lan = v_lan, motor = v_motor) |> purrr::list_rbind(names_to = "domain")`; the stacked table keeps
  its class and scale (tables on different scales stop, as with `rbind()`). `plot()` draws one panel per `domain`.
* `plot()` of a `kits_vertical()` has every feature of `kits_plot_cohort()`: the points at the actual median age,
  a note under the title when it differs from the wave by more than 1 month, the age in the hover text, `cohorts`
  (which and in which order: blue circles, orange squares), `age_ticks` and `y_range`.
* `kits_norms(..., by = )` fits one norm per group. Each group's norm is identical (tolerance 0) to `kits_norms()` on
  the children of that group alone. The result is one `kits_norms` table with the group columns in front;
  `attr(, "by")` names them, `attr(, "groups")` keeps each group's norm, the settings of `kits_info()` record `by`,
  and the diagnostics of each group are kept with the group's name. An error or a warning in one group names the
  group (an error keeps its class). With a `kits_subj()`, `by` names its columns; with scores as numbers, `by` is a
  data.frame of the group columns.
* `kits_norm_z()` with a norm by group matches each child to the curves of its group: a `kits_subj()` needs the
  group columns (a clear error otherwise); with numbers, give `by = data.frame(...)`. Children of a group without a
  norm get NA with a message. The values returned are unchanged.
* `plot()` of a norm by group: a row per domain (when `domain` is a group column), a column per other group, each
  group in its own colour; with `scores =`, each panel draws the children of its group behind the curves (the
  scores need the group columns).
* `kits_code()` of a norm by group: the overview fits each group with the three steps of the norm and runs as
  pasted (the same numbers and z); `full = TRUE` writes the whole script of each group in turn and stacks the
  norms; `part =` explains that it is for one norm.
* New `plot()` methods: a `kits_item` (category probability curves, one panel per item, or the expected item
  scores), a `kits_dif` (each item's DIF with its 95% CI from the new `attr(, "se")`, beyond `cut` in orange or
  blue), a `kits_itemfit` (infit and outfit, the acceptable range shaded), a `kits_age_pass` (easy to hard, the
  other `p` joined by a line, separated items open), a `kits_drift` (new against old thresholds with the line, or the
  distance from the line; drifted steps orange, steps suggested for fixing blue) and a `kits_norm_screen` (the share
  at the lowest score per wave, the waves left out orange, the cut dashed).
* Internals: no function of kits calls a deprecated function, so normal use never warns (`kits_drift()` and
  `kits_precision()` used `kits_mean_sigma()` and `kits_wle_table()`). The PCM category probabilities of all items
  come from one helper, used by the test characteristic curve (`kits_information()`, `kits_expected_score()`, IRT
  true-score equating) and by the Lord-Wingersky recursion (`kits_precision()`, observed-score equating);
  `kits_raw_equate(method = "true")` uses the same expected-score function as `kits_expected_score()`.
  `kits_precision()` does not use the test information (it uses the Lord-Wingersky recursion and the WLE table),
  so it shares the category probabilities rather than `kits_information()`.
  `kits_convert()` converts by an internal S3 generic (numbers, a `kits_item`, a table).
* Examples, vignettes and README use the new calls only; `kits_code(x, full = TRUE)` is byte-identical to the
  fixtures of 0.16.0, and the overviews run as pasted.
* Tests: `test-norms-by.R` (norms by group) and `test-plots-results.R` (the new plots); `test-deprecated22.R` checks that each deprecated function and form warns with
  `kits_warning_deprecated` and gives the same result as its replacement, the new forms of the plots, and that
  normal use gives no deprecation warning. The other tests use the new calls.

## Round 3: data first, models as tools

* `kits_resp(data, persons = NULL)` makes one table of responses: one row per person-wave, the columns of the
  child first (`person_id`, `pid`, `cohort`, `sex`, `wave_m`, `age_m`, those present), then one column per item;
  kits remembers which columns are items. From a long data.frame (one row per response, with the columns of the
  child) or a wide one with `persons =`. One domain at a time (a `domain` column with several stops and shows the
  `filter()`); a column of the child that varies within a person-wave stops. `print()` gives the counts per wave and
  the checks, `summary()` the standard form of each wave and the share of missing responses, `plot()` the design
  (which items each wave answered, coloured by the share passing). The functions that take responses take it, with
  the columns of the child used where they are needed; `kits_item(resp)` labels the scale `"<domain> <cohort>"`
  (with a message) when the table holds one domain and one cohort.
* Scores keep the raw score -> theta, SE of each booklet (summary numbers, no child's data) and the conversions
  applied since, so `scores |> kits_subj_table()` gives the table of each wave, identical (tolerance 0) to the old
  `kits_subj_table(item, data, wave_m)` followed by the same conversions. The counts come from the rows of the
  scores now, so a `filter()` is respected; two cohorts stacked with `rbind()` keep both records.
  `kits_subj_table(item)` (every raw score of every booklet) is unchanged.
* `kits_norm_z()` given scores returns the scores with `z` and `pct` (`100 * pnorm(z)`); given numbers, the numbers
  as before. The values of `z` are unchanged.
* `kits_norms(s)` without `ages` uses `kits_norm_ages(s)` and `pad = 0`, as the tutorials did; `print()` states the
  range and where it came from, notes waves that `kits_norm_screen()` would leave out, and gives the next steps.
  The recorded call keeps `ages` as a setting, so `update(norm, degrees = )` refits it.
* `kits_item_table(thr, norms, age_pass, dif, stems, rp)` is one table per item (D, the age at which a share `rp`
  pass, the age at which the norm median reaches the D, the DIF, the stem) that keeps the norm; `plot()` draws the
  growth chart with the items (`type = "growth"`) or the developmental ladder (`type = "ladder"`), the same plots as
  `kits_plot_growth()` and `kits_plot_ladder()`.
* `print()`: scores show their main columns and name the others, say what else they keep, and end with `Next:`;
  so do item parameters, links and norms. `kits_subj_table()` names the kind of table in its header (every raw score
  of every booklet, or the standard form of each wave), and `kits_age_pass()` the definition of a pass.
* Printouts after a beginner tutorial: a diagnostics table after `select()` prints only its columns (no
  stray `message` column); a `kits_subj` whose columns were selected or reordered prints them as given, without
  the `Other columns`, `Also kept` and `Next` lines (row filters keep them); a norm prints the age range it was
  fitted on (kept in `attr(, "fitted")`), and after `filter()` how many rows are shown; the note that the norms
  are not on D comes once, from `kits_norms()` (no longer from `kits_norm_screen()`), and says what to do;
  `kits_item()` reports an early nlminb stop that the Newton steps finished as converged, in plain words.
  One word for one thing: a booklet is any combination of items a person-wave answered, the standard form the
  most common booklet of a wave (prints, messages and help); `summary(resp)` prints the first item ids of each
  standard form and, for the other booklets, the items they lack or add (the full lists stay in
  `summary(resp)$forms`).
* New base methods (see `?kits_methods`): `confint()`, `vcov()`, `logLik()`, `nobs()`, `update()`, `simulate()`.
* Stacking: `rbind()` and `purrr::list_rbind()` check the scales and combine the records of all tables (the links
  without SE, the raw score tables); `dplyr::bind_rows()` cannot, so it warns (class `kits_warning_scale`) when the
  first table has records.
* `kits_linking(domains, anchors)`: domains without anchors (for example social) get the same line as the domains
  with anchors, as before; this is now printed, recorded in the settings of `kits_info()` (`no_anchors`) and written
  by `kits_code()`. KIT scale scores are in D units, not GSED D. The numbers of the domains with anchors are
  unchanged, and a link with no domain of at least 3 anchors still stops.
* Tests: `test-round3.R` (the response table, ids, the tables from scores against the old ones, z, the item table
  against the old functions, the base methods, stacking, anchors, domains without anchors).

## Data

* New data `kitm03` and `kitm36`: simulated KIT data (language and motor, 500 children per cohort), each a list of
  `responses` and the `anchors` key for linking to D; no row comes from a real child.

# kits 0.21.0

The WLE is computed in base R, without RTMB; the CML stays in RTMB.
The numbers of 0.20.0 are unchanged: the item parameters are identical, and theta and se differ by at most 2e-11.
On the real KIT data (five domains, two cohorts, 296,343 person-waves) the largest difference is 1.9e-11, and
scoring all ten calibrations takes 20.2 s instead of 26.0 s.

## WLE in base R

* Every WLE of kits (`kits_subj()`, `kits_wle_table()`, and through them `kits_itemfit()`, `kits_subj_table()`,
  `kits_precision()`, `kits_link_check()` and `kits()`) is now found by Newton steps in base R, for all item set x
  raw score combinations at once. Both derivatives of the weighted likelihood f = log L + 1/2 log I have closed forms
  from the PCM category probabilities: f' = r - E + I' / 2I and f'' = -I + (I''/I - (I'/I)^2) / 2.
  `RTMB::MakeADFun()` and `nlminb()` are no longer used for the WLE.
* Kept from 0.16.0 and 0.18.0: the grid search that keeps the higher of two maxima of the weighted likelihood
  (the lower theta when the two are equally high, as in dexter), the finite WLE for zero and perfect scores, and
  Warm's SE sqrt(I + (I'/2I)^2) / I.
* The grid search now computes log S and I once per item set on one theta grid (they do not depend on the raw
  score) and sums the items with a matrix product; this was most of the time of the WLE. The grid covers the same
  range (every 0.1, within 15 of the start value of each combination), on points shared by all combinations.
  A combination with two valleys on the grid now compares its current maximum with both, so a tie gives the lower
  theta whatever maximum the Newton steps reached first (in 0.20.0 this depended on the start value).
* Equal to 0.20.0 and to dexter: on the four benchmark designs below, theta and se differ from 0.20.0 by at most
  3.4e-13 and 3.9e-13 (infit and outfit by 2.4e-13; the flags are identical); on 120 random item sets (2 to 8
  items with 1 to 3 steps, a quarter of them with items far apart), by at most 1.2e-11 (at a zero score, where
  0.20.0 was the less precise), and the largest difference from `dexter::ability_tables(method = "WLE")` went from
  1.2e-11 to 2.2e-12.
* New tests (`test-wle21.R`): equality with the WLE tables of 0.20.0 (stored in `fixtures/wle020.rds`) on several
  booklets, items far apart (two maxima), fixed parameters with a left-out item and betas shifted by 60; the score
  r - E + I'/2I is 0 at every estimate, including zero and perfect scores; no point of a fine grid has a higher
  weighted likelihood; equal maxima give the lower theta; dexter on 30 random item sets (within 5e-8).

## CML stays in RTMB

* An analytic gradient of the conditional likelihood in R was written and timed (the gradient of each step is
  -n_ig + sum_b sum_r n_br P(X_i >= g | b, r), with P from the ratios of the elementary symmetric functions
  gamma^(i)_{r-k} / gamma_r, carried through the items by one forward and one backward pass over the booklets;
  then the exact Hessian of RTMB at the end). It reached the same estimates (beta within 3e-14, SE within 8e-14),
  but it was slower on every design: on the expert designs one evaluation of the R gradient takes 13 to 17 times
  as long as the RTMB gradient of the taped likelihood (0.03 s and 0.20 s), and nlminb needs 180 to 280 of them:

  | Design | Booklets | CML, 0.20.0 (RTMB) | CML, R gradient |
  |---|---:|---:|---:|
  | KIT-like, 46,900 person-waves, 24 items | 59 | 0.04 s | 1.2 s |
  | KIT-like, 46,900 person-waves, 24 items | 132 | 0.07 s | 0.6 s |
  | 40 items, N = 3000, 2% missing | 514 | 2.5 s | 7.8 s |
  | 40 items, N = 3000, 10% missing | 2,689 | 12.2 s | 51.8 s |

  So `kits_cml()` and every CML in kits are unchanged: `nlminb()` on the RTMB gradient, then Newton steps with the
  exact RTMB Hessian.
* A new test checks the gradient of the CML likelihood three ways on a small design with two booklets: the RTMB
  gradient, a numerical gradient, and the formula above with P(X_i >= g | b, r) from the enumeration of all
  response patterns.
* RTMB stays in Imports: it computes the CML (`kits_item()`, `kits_cml()`, `kits_drift()`, `kits_dif()`,
  `kits_link_check()`), the shash norms, the interaction model, and it is used by the scripts of `kits_code()`.

## Time (one thread, median of 3 runs)

Simulated data (`kits_sim()`). KIT-like: 7 waves of 6,700 person-waves, each with a targeted form of 7 or 8 of
24 items (1 to 3 steps), plus random missing responses. Expert review: 40 items (20 with 2 steps), N = 3000,
2% or 10% of the responses missing at random. dexter 1.8.1: `fit_enorm()` and `ability(method = "WLE")` on the
long data.frame (kits takes the matrix and converts it inside the timing).

| Design | Booklets | `kits_item()` 0.20.0 | 0.21.0 | dexter `fit_enorm()` | `kits_subj()` 0.20.0 | 0.21.0 | dexter `ability()` |
|---|---:|---:|---:|---:|---:|---:|---:|
| KIT-like | 59 | 1.66 s | 1.65 s | 0.13 s | 2.12 s | 1.97 s | 0.06 s |
| KIT-like | 132 | 1.71 s | 1.69 s | 0.19 s | 2.72 s | 1.96 s | 0.06 s |
| 40 items, 2% missing | 514 | 2.90 s | 2.91 s | 1.27 s | 4.83 s | 0.47 s | 0.07 s |
| 40 items, 10% missing | 2,689 | 12.81 s | 12.76 s | 5.65 s | 10.81 s | 0.75 s | 0.24 s |

The WLE alone (all combinations of a design): 0.17 -> 0.03 s, 0.65 -> 0.06 s, 4.77 -> 0.22 s and 9.43 -> 0.49 s.
In the KIT-like designs the rest of `kits_subj()` (about 1.9 s) is the data handling: the item set of each
person-wave (0.8 s) and the person fit (0.9 s).

## The scripts of `kits_code()`

* The whole scripts (`full = TRUE`) are unchanged and still use RTMB for the WLE (`nlminb()` on the RTMB gradient,
  then Newton steps); the fixtures of 0.16.0 are identical, header md5 lines included. The script and the package
  now take different routes to the WLE and reach the same maximum.
* Documentation: `?kits_subj`, `?kits_code`, `?kits_equations`, the package page, the README, the flow vignette and
  the Description say that the WLE is computed in base R, and that RTMB computes the CML, the norms and the
  interaction model (and the scripts). `kits_equations()`: the WLE sections say "in base R" and cite R instead of
  RTMB.

# kits 0.20.0

kits and kitsPreview now share one version number, starting at 0.20.0 (kits skips 0.18.1 and 0.19; kitsPreview jumps from 0.11).

A fix of the linking error in a pipe of `kits_convert()`, and documentation.
The numbers of 0.18.0 are unchanged.

## Linking error in a pipe

* `x |> kits_convert(co) |> kits_convert(d_link)`, with `co` made with `n_boot = 0` and `d_link` with draws, added
  `se_link` (and `se_total`) that held the error of `d_link` only, so the linking error was understated.
  `kits_convert(x, co, d_link)` in one call rightly added no column.
  Now a converted table records the links whose error it does not carry (`attr(, "links_without_se")`: links without
  draws, and every link applied to a `kits_item()`, recorded on its thresholds so that `coef()` keeps it).
  A later `kits_convert()` that sees the record adds no `se_link` either, and gives the message of the one-call form
  (class `kits_message_scale`), naming the earlier links.
  A table that already has `se_link` loses it, with the message, when the next link has no draws (as before).
  The record is kept by `[`, the dplyr verbs and `rbind()`; `dplyr::bind_rows()` keeps only the record of its first
  table (dplyr rebuilds the result from the first table), so stack converted tables with `rbind()`, or convert each
  table with all its links before stacking.
  New tests check that the pipe and the one call agree on whether `se_link` and `se_total` exist, and on their values,
  for scores, thresholds and numeric vectors, in both orders.

## Documentation

* Help pages brought in line with the code:
  * `kits_drift()`: items whose new data have no 0 score or one score only are also suggested for fixing;
    `points` has the column `judged`.
  * `kits_item()`: the result is the same as `kitsPreview::kits_item()` (its `dataSrc` is `data` here); `fit` also has
    `inputs$left_out` and `optim`; the Details name the stops for steps without a CML estimate (`kits_error_items`)
    and for booklets in separate groups (`kits_error_data`), and the diagnostic "extreme estimates".
  * `kits_check()`: the list of checks adds "extreme estimates", the age-separated items of `kits_age_pass()` and the
    links' domains without kept anchors.
  * `kits_norms()`, `kits_norm_fit()`: the fit stops (`kits_error_data`) with fewer than `10 * sum(degrees + 1)`
    person-waves inside the ages. `kits_norm_data()`: `attr(, "norm")` also holds `age_x_text` and `scale`.
  * `kits_linking()`, `kits_mean_sigma()`: `settings$key` for anchors from `kits_anchors()`; with `spread = "se"`,
    `points` has `se_from` and `se_to`; `seed = NULL` gives a run that cannot be repeated (also `kits_link_error()`).
  * `kits_link_check()`: repeated ids in `persons` stop; the result has `settings` (`tol`); the columns of `waves`.
  * `kits_age_pass()`: NA also for items answered at fewer than 2 distinct ages.
  * `kits_inter()`: every element of `est` and `inputs`.
  * `kits_subj()`, `kits_wle_table()`, `kits_precision()`, `kits_raw_equate()`: `item` may be a `kits()`;
    `kits_subj()`: other columns of `persons` are dropped with a message.
  * `kits_vertical()`: the shrinkage check uses `sd` when the scores have no `se`; waves whose `sd_true` is 0 are noted.
  * `AGE_TICKS` and `NORM_PCT` say which defaults they are.
  * The package page: `kits_warning_anchors` and the `kits_convert()` message in the list of condition classes;
    `plot()` of `kits_vertical` and `kits_precision`; the person information of `kits_link_check()` and
    `kits_norm_screen()`.
* README and vignettes: the link steps equal `kits_linking()` except `info`; the three norm steps give the same numbers
  as `kits_norms()` but record no scale, so `kits_norm_z()` cannot check it; the norm overview of `kits_code()` ends
  with `kits_norm_z()`; the link overview has the fourth line `kits_link_error()` only for a link with draws;
  `kits_norm_data()` in the table of renamed arguments; one line each on `AGE_TICKS` and `NORM_PCT`.

# kits 0.18.0

Links and norms in named steps, one-row summaries with `as.data.frame()`, and the design items of the expert review
(Part B): overlapping person ids in `kits_link_check()`, a faster CML, a flag for perfect separation, consistent
argument names, and shorter code lines.
The numbers of 0.17.0 are unchanged, except the last digits of the CML (below 1e-10; see "CML with many booklets").

## Links and norms in named steps

* `kits_linking()` is four exported steps, each with an example:
  `kits_link_points()` (one point per common step, or per anchor) -> `kits_link_screen()` (the line of all points and
  one round of screening) -> `kits_link_line()` (mean/sigma on the kept points, with the domain check) ->
  `kits_link_error()` (the item bootstrap). `kits_linking()` and `kits_mean_sigma()` run the same steps, so the four
  steps give an identical link (every element but `info`).
* `kits_link_screen()` also takes a plain data.frame of points (`id`, `item`, `from`, `to`), so points from elsewhere
  go through the same steps. `kits_link_error(link, n_boot = 0)` removes the linking error.
* `kits_norms()` is three exported steps: `kits_norm_data()` (the data inside the ages) -> `kits_norm_fit()` (the shash
  by maximum likelihood, the best start kept) -> `kits_norm_table()` (the parameters and percentiles at each age,
  with the diagnostics). The three steps give the same norm table as `kits_norms()`.
* The overviews of `kits_code()` for links and norms call these steps, one step per line, and run as pasted
  (the tests check that they give exactly the numbers of the result).
* Shorter parts of `kits_code()`:
  * `kits_dif()`, `part = "cml"`: the CML of one group, the same script as `kits_code(item, part = "cml")`
    (86 lines; was 127).
  * `kits_inter()`: `"form"` (31 lines), the new `"stats"` (the sufficient statistics of one wave, 55 lines) and
    `"model"` (the interaction model of that wave, 81 lines; was 176). The `"model"` part has the sufficient
    statistics written in, so it needs no data.
  * `kits_norms()`: `"fit"` and `"table"` take the steps before from kits (`"table"`: 46 lines; was 91).
* The full scripts (`full = TRUE`) are unchanged. In three link fixtures of the tests, the data md5 in the header
  changed, because the thresholds moved in their last digits (see below); every other line is identical.

## `as.data.frame()`: one-row summaries

* New `as.data.frame()` methods give a plain data.frame with one row (or one row per wave, method or booklet),
  ready for `dplyr::bind_rows()`; `checks` is the most serious status of the diagnostics.
  The columns are documented in `?as.data.frame.kits`.
* A link (`kits_linking()`, `kits_mean_sigma()`): `from`, `to`, `kind`, `a`, `b`, `a_error`, `b_error`, `cov_ab`,
  `a0`, `b0`, `n_points`, `n_dropped`, `n_items`, `n_failed`, `cut`, `spread`, `unit`, `step`, `n_boot`, `seed`, `key`,
  `checks`.
* `kits_drift()` (one row), `kits_inter()` (one row per wave), `kits_link_check()` (one row per method), `kits_info()`
  (the provenance of any result, one row), `kits_item()` and `kits()` (one row), `kits_norm_screen()` (its waves),
  `kits_norm_fit()` (one row), `kits_stats()` (one row per booklet).

## B1: overlapping person ids in `kits_link_check()`

* New arguments `persons_from` and `persons_to`, one person table per cohort, so the two cohorts may use the same
  person ids. `persons` (one table) still works when the ids do not overlap; with repeated ids it now stops
  (`kits_error_input`) and points to `persons_from` and `persons_to`.

## B2: CML with many booklets

* `kits_item()`, `kits_cml()` and every CML inside kits (`kits_drift()`, `kits_dif()`, `kits_link_check()`):
  * the conditional likelihood is taped for all booklets at once (one matrix row per booklet), so RTMB records it
    about five times faster; it is the same formula as the per-booklet loop of the scripts;
  * `nlminb()` uses the RTMB gradient only (quasi-Newton), and the exact Hessian is computed for the Newton steps
    at the end (two or three times instead of at every iteration).
* Time of `kits_item()` (one thread):

  | Design | Booklets | Parameters | 0.17.0 | 0.18.0 |
  |---|---:|---:|---:|---:|
  | 40 items, N = 3000, 2% missing (expert review) | 511 | 79 | 20.8 s | 6.4 s |
  | 40 items, N = 3000, 10% missing (expert review) | 2,662 | 79 | 90.0 s | 26.6 s |
  | KIT-like: 24 items, 47,000 person-waves | 56 | 33 | 3.1 s | 2.8 s |
  | KIT-like: 24 items, 47,000 person-waves | 139 | 33 | 3.2 s | 2.6 s |

* The estimates agree with 0.17.0 to the tolerance of the optimizer: on these four designs the largest
  |difference| is 4e-14 in `beta` and `thr`, 9e-12 in `SE_beta` and `thr_se`; on the small data sets of the
  tests, at most 2e-15 in `beta` and 8e-12 in the SEs.
* The diagnostic "optimizer convergence" reports the code of the quasi-Newton run, which is often 1
  ("false convergence (8)" or "singular convergence (7)") at the very strict `rel.tol = 1e-14`; as before, the
  status is "ok" when the gradient is below 1e-4 and the Hessian is positive definite after the Newton steps.

## B4: perfect separation in `kits_age_pass()` and `kits_age50()`

* `kits_age_pass()` returns a table (class `kits_age_pass`, a data.frame) with one row per item: `item_id`, one
  column per `p` (`p50`, ...), and the new columns `separated` (`glm()` did not converge, or its fitted
  probabilities reached 0 or 1), `age_last_fail` and `age_first_pass` (the bracket where the age is not
  identified). The ages are exactly those of 0.17.0. Its diagnostics count the separated items (a note).
  Before, it returned a matrix (a named vector for one item).
* `kits_age50()` still returns the named vector of ages, and warns (`kits_warning_separation`) when items are
  separated.
* `glm()` is no longer wrapped in `suppressWarnings()`: its two separation warnings become the flag, and any other
  warning is shown.

## B5: argument names and S3 dispatch

* Every argument name is lower snake_case. The old names still work in 0.18.0, with a warning of class
  `kits_warning_deprecated` that names the new argument, and give identical results. Giving both names stops.

  | Old | New | Functions |
  |---|---|---|
  | `dataSrc` | `data` | `kits()`, `kits_item()`, `kits_stats()`, `kits_subj()`, `kits_subj_table()`, `kits_itemfit()`, `kits_dif()`, `kits_inter()`, `kits_drift()` |
  | `D` | `d` | `kits_norms()`, `kits_norm_z()`, `kits_norm_screen()`, `plot()` of a `kits_norms` (`plot.kits_norms()`) |

* The old names sit after `...`, so R matches them exactly and never partially (no clash between `D` and `d`, or
  between `dataSrc` and `data`). The functions that gained `...` for this stop on an argument they do not have
  (`kits_error_input`) instead of ignoring it.
* `kits_norms()` and `kits_norm_z()` dispatch on their first argument (S3): `kits_norms(d, age_m, ages, ...)` and
  `kits_norm_z(d, age_m, norm)` for scores; `kits_norms(s, ages, ...)` and `kits_norm_z(s, norm)` for a
  `kits_subj()`. The second argument no longer shifts: `kits_norms(s, age_m = x)` stops with a clear error
  (it used to treat `x` as `ages`).
* In the scripts of `kits_code()`, the scores keep the name `D` of the formulas (`D` and `age_m`, as before).
* Columns named `D` in tables (for example the items of `kits_plot_ladder()`) are data, not arguments, and keep
  their name.

## B6: style

* No code line is longer than 130 characters, except the first line of six script templates written by
  `kits_code(full = TRUE)` (they must stay byte-identical); long strings are split with `paste0()`.
* No `if`/`else` inside the arguments of a diagnostic row: the status and the message are computed first, one
  concept per line (37 places).

## Other changes

* A few errors on wrong input now come in another order (for example a wrong `cut` before a missing anchor), and
  `kits_norm_data()` stops when fewer than 2 person-waves are inside the ages.

# kits 0.17.0

Three changes: `kits_convert()` is taught with the pipe, the D-link anchors are package data, and `kits_code()` has three tiers.
No number of 0.16.0 changed (every object compared in the tests).

## `kits_convert()` with the pipe, one link per line

* The docs, the README and both vignettes use the pipe: `s36 |> kits_convert(co) |> kits_convert(d_link)`.
* `kits_convert(s36, co, d_link)` is the same thing in one call, with the same result (`theta`, `se`, `se_link`, `se_total`, scale).
* A converted numeric vector records its new scale (`attr(, "scale")`), so the next link in a pipe checks it; before, the second link warned "scale not recorded".
* The wrong order still stops in a pipe (`kits_error_scale`).
* The name stays: Kolen and Brennan write "Scale I parameter estimates are converted to Scale J"; `kits_linking()` finds the line, `kits_convert()` applies it.

## `kits_anchors()`: the anchors of the D link

* New function `kits_anchors(key = c("gsed2406", "gsed2510"))` returns the anchor table of one dscore key: `item_id`, `gsed_item`, `key`, `tau`, `note`.
* gsed2406 (the default) has 29 items; gsed2510 has 8. The table is stored in the package (`R/sysdata.rda`).
* `tau` comes from `dscore::builtin_itembank` (dscore 2.0.0, Apache License 2.0); matching KIT items to GSED items is KIT's own decision, recorded in `note`.
* gsed2212 has exactly the same `tau` as gsed2406; dscore marks gsed2212 as retired, so the key is called gsed2406.
* `key = "gsed2510"` warns (`kits_warning_anchors`): only 8 items have a `tau` (2 language, 6 motor, 0 cognitive), and its D metric differs, so D values of the two keys are not interchangeable.
* An anchor table with a duplicated `item_id` (for example two keys bound together) stops `kits_linking()` (`kits_error_input`) and asks to choose one key.
* The link records the key: `link$settings$key`, `kits_info()`, and the printed "points from:" line.

## `kits_code()` in three tiers

* Overview (the default): about 30 lines, one line per step; each comment points to the detail of its step, for example `-> kits_code(item$stats)` or `-> kits_code(item, part = "cml")`.
* The overview uses `library(kits)` and exported functions only, so a student can paste it and get exactly the same numbers; its header names the inputs needed (for example `X`).
* One step (`part =`): the detail of one step, runnable alone from the inputs in its header.
* Item parts: `"stats"`, `"cml"`, `"thr"`; scoring parts: `"wle"`, `"person_fit"`.
* Links: `"line"`, `"error"`, `"domains"`; norms: `"data"`, `"fit"`, `"table"`, `"z"`; DIF: `"cml"`, `"dif"`, `"overall"`; interaction model: `"form"`, `"model"`.
* The full script (`full = TRUE`) is the `kits_code()` of 0.16.0 and needs no kits; no section header changed.
* The script of `kits_norm_screen()` is short, so its overview is the full script.
* `kits_rtmb_code()` stays an alias with the same arguments.
* For the full script of 0.16.0, add `full = TRUE`: `kits_code(item, subj = TRUE, full = TRUE)`.

## Calibration in three steps; the result packs key objects

* `kits_item()` reads as three steps, one line each: `kits_stats()` (sufficient statistics), `kits_cml()` (CML) and `kits_thresholds()` (thresholds); all three are exported, with examples.
* `item$stats` (class `kits_stats`): `n_ig` of each step, `n_br` of each booklet and raw score, and the booklets; it has a print method and holds counts only, no person ids.
* `kits_code(item$stats)` writes how they are counted (no kits needed).
* `kits_cml(item$stats)` is identical to `item$fit`; `item$fit$inputs$left_out` records the fixed items left out of the calibration.
* New functions `kits_information(item)` and `kits_expected_score(item)`, like dexter's `information()` and `expected_score()`, return a function of theta; give `items` or `booklet_id` to choose the items.
* Each has its own script (`kits_code(information)`), with the item parameters written into it.
* The sufficient statistics of the CML are counted with dplyr (before: `tapply()`, `mapply()`); the numbers are identical.

# kits 0.16.0

Fixes from three reviews (an expert review, a beginner's review and the privacy findings of the Shiny plan).
The numbers of 0.15.0 are unchanged (to about 1e-13), except where a fix below says so.

## Numbers that change on purpose

* `kits_inter()`: `SE_sigma` is the SE from the full Hessian (all parameters, with the mean-0 normalization of sigma),
  no longer dexter's block SE of each item, which was about 16% too large; `fit_IM = sigma / SE_sigma` grows by the
  same factor. In a simulation under the PCM (6 items, 500 children, 150 replicates) the SD of sigma was 0.059, the
  new SE 0.060 (block SE 0.070), and |fit_IM| > 2 flagged 4.3% of the items (1.7% before; 5.0% at 1.96).
  The block SE is kept in `fits[[w]]$est$se_sigma_block`. `sigma` is unchanged.
* `kits_dif()`: the overall test centres the compared steps of each group, so it compares the groups on a common
  origin. With an item left out (answered in one group only, or with different scores), 0.15.0 counted the
  difference of the two origins as DIF and rejected almost always (50 of 50 replicates with no DIF); now the type I
  error at 5% is 0.055 (200 replicates; 0.050 with no item left out). Without a left-out item the test is unchanged
  (and equals dexter's `DIF_overall`); with one it now differs from dexter's `DIF()`, which keeps the uncentred test.
  The size `dif` is unchanged in every case.
* `kits_norms()`: the norm table, its diagnostics and `kits_norm_z()` use the clamps of the fit (log sigma and epsilon
  within -10..10, log delta within -2..2). Where a polynomial goes beyond its clamp, 0.15.0 reported a distribution
  other than the fitted one (in the reviewer's example delta 12.2 instead of 7.39, and 17% of the children aged over 30
  outside P3-P97 instead of 6%). Elsewhere nothing changes. This departs from cNORM, whose
  `predictCoefficients_shash()` does not clamp.
* `kits_norms()`: when a restart of the optimizer lowers -log L by 1e-6 or more, its better solution is kept (0.15.0
  flagged the fit and kept the worse one). Fits that converged are unchanged.
* Booklet ids: `kits_item()$design` and `kits_wle_table(item)` use the ids of `kits_subj()` (`bk1`, `bk2`, ..., persons
  sorted in the C locale, booklets in order of appearance) instead of `set1`, `set2`, ..., so the tables can be joined
  on `booklet_id`. The rows of `kits_wle_table(item)` follow the new ids; every value per booklet is unchanged.
* `kits_sim()`: the default person ids get a prefix made from the simulated data (`"p3fa2_0001"`), so two calls no
  longer give the same ids. The scores are unchanged; `prefix = "p"` gives the ids of 0.15.0.
* md5 fingerprints of data with non-ASCII ids now ignore the encoding mark of the strings; ASCII data keep their md5.
* WLE (`kits_subj()`, `kits_wle_table()`, `kits_subj_table()`, `kits_precision()`): when the weighted likelihood of an
  item set and raw score has two maxima (few items far apart, a middle raw score), 0.15.0 could stop at the minimum
  between them (its start value), for example theta = -0.05 instead of -0.57 or 0.47 for two dichotomous items 4.3
  logits apart and a raw score of 1. Now the higher maximum is kept, and with a tie the lower theta, as in dexter.
  This was found by the agreement test with kitsPreview, which failed in about 5% of its random data sets
  (incomplete designs with two-item booklets); kits now agrees with dexter to 5e-8 on 120 such data sets.
  Item sets with one maximum are unchanged.

## Privacy

* `kits_info()$call` never holds data. The call keeps the names of the objects and the settings; an argument given
  by value (`do.call()`, `purrr::map()`) and the constants inside an expression given as data are written `<data>`.
  `kits_check()` labels an object given by value "object 1", ... (it deparsed the object before).
* The default scale label is the name of the data object; when the data are not a named object (an expression such as
  `X[keep, ]`, or a value), the label is `"scale 1"`, with a message (class `kits_message_scale`). Give `scale`.
* Messages and errors give the number of person-waves concerned, not their ids (dropped person-waves, repeated ids,
  missing rows in `persons`, values of `wave_m` or `sex` by person id). `options(kits.show_ids = TRUE)` shows ids again.
* `plot(norm, D = , age_m = )` draws the children binned by default (about 40 x 40 cells; cells with fewer than
  `min_cell = 5` children are not drawn), so the plot and its html hold no single child. `persons = "points"` draws
  every child, as before.

## Stacking tables

* `dplyr::bind_rows()` of kits tables checks the scales like `rbind()`: tables on different scales stop
  (`kits_error_scale`); tables on one scale keep the class and the scale; a plain data.frame or tibble mixed in gives
  `kits_warning_scale`. 0.15.0 silently labelled the result with the scale of the first table, so `kits_convert()`
  then converted half of the rows wrongly. (vctrs methods; vctrs is a new import, already a dependency of dplyr.)
* After `rbind()` or `bind_rows()`, `summary()` and `kits_check()` work: a stacked `kits_subj` gets its diagnostics
  recomputed (with the reliability of each wave); other tables stacked with `rbind()` keep the diagnostics of each
  table ("table 1: ...").
* `kits_check()` no longer stops on an object without diagnostics: a table that lost them gets a note, an object that
  is not a kits result a warning row, and the other objects are still checked (no purrr backtrace).

## CML

* The elementary symmetric functions are computed on the betas centred within each booklet, so long booklets and
  fixed items far from the origin no longer underflow (a 40-item booklet about 6 logits from the origin stopped with
  "infinite or missing values in 'x'"). The same in the scripts of `kits_code()`.
* An item step with no finite CML estimate (reached, or missed, by none of the person-waves whose raw score is neither 0
  nor the top score, for example an item passed only by children with a perfect score) stops with `kits_error_items`;
  0.15.0 returned betas of -7 to 23 with SEs of 23,000 and "Checks: OK". A new diagnostic, "extreme estimates", warns
  when a free step has an SE above 10 or |beta| above 15.
* A design that CML cannot identify stops with `kits_error_data` and says why: item sets that form separate groups with
  no common item (or, with fixed items, a group without a fixed item), or one item per person-wave. The Newton steps
  check the Hessian before solving.

## Other fixes (expert review)

* `kits_convert()` of a `kits_item` records `beta_scale` on its thresholds too, so `kits_subj(X, coef(item_D))` labels
  the scores with the calibration scale (0.15.0 labelled logits "D").
* `kits_linking()` with domains no longer crashes when a domain loses all its anchors and `n_boot > 0`.
* `kits_equations()`: the `sd_true` note of `kits_vertical()` follows `extremes`; the chain of links no longer says
  "M36 -> M03"; the KIT threshold is written thr_ik (tau is the GSED difficulty) and the true score T_F; the DIF and
  interaction-model formulas follow the fixes above.
* Input checks: `n_boot` must be a whole number; `kits_norm_z()` needs positive `sigma` and `delta`; `kits_norms()` stops
  when all scores are the same; `kits_raw_equate(method = )` gives a classed error; a factor of scores is named;
  `kits_vertical()` with one child per wave gives no base warning.
* `kits_norms()` documents that the table runs from `floor()` to `ceiling()` of the observed ages.
* `pcm_probs()` is defined once.
* The docs of `kits_link_se()` and `kits_convert()` say that the links' resamples are treated as independent, which is
  an approximation when links share items, and that the chain variance is the delta method centred at the estimates;
  `kits_link_se()` names shared items in a message.
* The bootstrap uses fixed RNG kinds (Mersenne-Twister, Inversion, Rejection), so a `seed` gives the same draws in any
  session; the caller's RNG kind and state are restored. With the default kinds nothing changes.
* DESCRIPTION: `kits_code()` covers the main estimates, not every result; its error lists the classes it covers.
* `kits_drift()` diagnostics give the expected number of steps flagged by chance.

## Messages and defaults (beginner review)

* Anchors with other column names: the error says that an anchor table needs `item_id` and `tau`, and a thresholds
  table `item_id`, `item_score`, `thr`.
* Scores coded 1, 2, 3: the error says to subtract 1 so that scores start at 0 (not to leave the items out).
* `persons` without the recognized columns lists the columns it has; unused columns give a message.
* `kits_dif()`: `ref` not in `sex` lists the values of `sex`.
* `kits_dif_combine(cohorts = NULL)`: the names default to the argument names (or "1", "2"), not M03 and M36.
  `kits_plot_cohort(cohorts = NULL)` and `kits_plot_norms(sexes = NULL)` take them from the data.
* `kits_sim()`: `prefix` is documented (give each cohort its own); the default ids no longer collide.
* The norm warnings name the `degrees` argument and suggest `c(3, 2, 1, 1)` for few waves.
* `kits_save_plot()` checks for pandoc first and says how to set `RSTUDIO_PANDOC` (class `kits_error_input`).
* `print(item, n = 20)` prints 20 steps; fixed steps are marked (`fixed` column and a header line).
* The convergence check of `kits_item()` says "converged" when nlminb ends with code 1 at its strict tolerance.
* A negative or low reliability of a wave says why (items too easy or too hard) and what to look at.
* `kits_convert()` says which links have no draws when `se_link` cannot be added.
* `kits_info()` prints "person ids md5" and "Persons table".
* `?kits_age50` calls the function in its last example; `?kits_item` shows `fixed_params` and uses 400 persons;
  the help titles of the report plots no longer carry the report letters "Plot A" to "Plot D".

## Documentation

* README reordered for beginners: what kits is, installation (all dependencies, pandoc, webshot2, a private library),
  a 10-line first example, a cheat sheet, the workflow (steps 0 to 10, as in the vignettes), a glossary, common errors,
  "my data are not KIT", and the advanced topics (Kolen & Brennan) last. `?kits-package` has a glossary and a runnable
  example.
* The vignette is split: `vignette("kits-flow")` (steps 0-6, simple data without DIF) and the new
  `vignette("kits-advanced")` (steps 7-10: items and DIF, a new wave, formulas and code, Kolen & Brennan).
* Terms: booklet (作答組合) is the item set a child answered, the standard form (題本) is the item set of most children of
  a wave.

# kits 0.15.0

Code style: tidyverse internals; dependencies added. The numbers of 0.14.0 are unchanged, apart from the fixes below.

## Code style

* The data handling of the package (long data, booklets, sufficient statistics, joins, summaries, the tables of the
  results) is written with dplyr, tidyr and purrr and the native pipe `|>`, in small functions with one step per line.
  `aggregate()`, `tapply()`, `merge()`, `ave()`, `split()`/`lapply()`/`do.call(rbind, ...)` and manual `match()` joins
  are gone from the data handling. dplyr, tidyr and purrr are new dependencies (Imports).
* The numerical cores stay in base R and RTMB, as they mirror the formulas or run in inner loops: the CML with the
  elementary symmetric functions, the WLE and its Newton steps, the Lord-Wingersky recursion, the shash likelihood,
  the interaction-model likelihood, the bootstrap of a link, the conversion of a response matrix to long data and the
  matrix algebra.
* Tables that are sorted by item, person or domain ids use the order of `sort()` (the collation of the session), as
  before, so every table and every bootstrap draw keeps its order.
* The results are still data.frames (`kits_table` and the other classes), not tibbles. Every data, `persons` and
  table input accepts a tibble.
* dplyr verbs on kits tables keep the class and the scale (a new `dplyr_reconstruct()` method), and drop the
  diagnostics and the provenance, as `[` does. For example `filter(s, flag == "ok")` is still a `kits_subj` on its scale,
  and `filter(coef(item), item_id %in% fix)` can be given as `fixed_params`.
* The scripts of `kits_code()` handle their data the same way (`library(dplyr)`, `tidyr`, `purrr` in the header); the
  estimation stays in RTMB and base R with the names of `kits_equations()`, and the scripts reproduce the package.
* The examples and the vignette use dplyr, tidyr and purrr for the data handling; kits tables are still stacked with
  `rbind()`, which checks their scales.

## Fixes

* `kits_vertical()`: new argument `extremes` (default `TRUE`, every child, as in 0.14.0). `extremes = FALSE` leaves the
  zero and perfect scores out of `sd_true`. In simulations both rules underestimate the true SD when a wave has a floor
  or a ceiling, and leaving the extremes out underestimates it more, so the default stays; the docs say to read
  `sd_true` of such a wave as a lower bound.
* `kits_precision()`: the default `criterion` is 1 logit (times the slopes of the links) instead of 0.5. Forms of 6 to
  8 items, as the standard forms of KIT, have a smallest RMSE of about 0.6 to 0.7 logit and never reached 0.5. An RMSE
  of 1 logit is a reliability of 0.75 among children with an SD of 2 logits; give `criterion` for another rule.
* `kits_linking()`, `kits_mean_sigma()`: the documentation says that `spread = "se"` suits small samples only (in KIT,
  `thr_se` is about 0.015 while the residuals are about 0.2, so it dropped 39 of 54 steps), and that `"sd"` (the
  default) or `"mad"` is the rule at KIT sizes. A new diagnostic warns when the screening keeps fewer than half of the
  points.
* `kits_equations("link", format = "latex")`: the note on a chain of links ends in proper math,
  `$\mathrm{Var}(b_{k})$ $v_{k-1}$`, instead of escaped text.
* `kits_equations("dif", format = "latex")`: the overall test says "without its first step", as the Unicode version and
  the code do.

# kits 0.14.0

The second batch of methods from Kolen and Brennan (2014), *Test equating, scaling, and linking* (3rd ed.).
Page numbers below are the printed pages of the book. The defaults are unchanged, so every result of 0.13.0
(and the agreement with kitsPreview) stays the same.

## New features

* New `kits_link_check(link, item_from, item_to, X_from, X_to, persons, theta, tol)`: independent checks of a link
  between two calibrations (pp. 190, 238, 295-296, 431-432).
  * Same-age groups: at the waves both cohorts share (36 and 48 months in KIT), the ability distribution of each
    cohort is estimated on its own calibration as a normal distribution by marginal maximum likelihood with the item
    parameters fixed (the raw score is sufficient; 121 quadrature points; `optim()`), with the SEs of the mean and SD.
    The line matches the means and SDs (b = sqrt(mean sigma_to^2 / mean sigma_from^2), every shared wave weighted
    equally), one line per shared wave and one for all of them. It needs no common items, but assumes no cohort effect.
  * Concurrent CML: the stacked data of both cohorts are calibrated together with `kits_item()`; each separate
    calibration is put on the concurrent scale by mean/sigma of all its thresholds, and the implied link is their
    composition (b = b_f / b_t, a = (a_f - a_t) / b_t).
  * `at` compares the lines at chosen theta (by default the wave medians of `from`): the value, the difference from the
    link, `se_link` and `z` = diff / se_link when the link has draws, and for the same-age lines their own sampling
    error (`se_method`, delta method) and `z_both`. Notes when a method differs by more than 2 (or by `tol` without
    draws). `print()`, `summary()`, `kits_check()`, `kits_info()` and `kits_equations()` (new topic `"link_check"`).
* Screening rules (pp. 190, 202-203, 288-289): `kits_linking()` and `kits_mean_sigma()` have `spread = c("sd", "mad", "se")`
  and `unit = c("step", "item")`.
  `"mad"` uses cut x 1.4826 x MAD(e); `"se"` standardizes each residual by sqrt(se_to^2 + b0^2 se_from^2)
  (`thr_se` of both calibrations; anchors have no SE, so only `from`; `kits_mean_sigma()` takes `se_from`, `se_to`);
  `unit = "item"` drops every step of an item with a flagged step.
  The link records the rule (`settings$spread`, `settings$unit`), `points` gains `z` and `flagged`,
  the print shows the rule, `plot()` draws the band of the rule, and `kits_equations()` and `kits_code()` follow it.
* New `kits_link_screening(link, theta)`: the line before screening (a0, b0) against the final line at chosen theta,
  the shift, and the shift in linking errors when the link has draws.
* Domain balance (p. 239): a link with domains prints the kept points per domain and notes when one domain has more
  than 2/3 of them (the common items should be a miniature of the test; the KIT D link keeps 22 motor of 25 anchors).
* Response probability (pp. 414-419): `coef(item, rp = 0.8)` gives the thresholds at which P(X >= k) = rp, with the
  delta-method SE, and records `attr(, "rp")` (kept by `[`, checked by `rbind()`). The default `rp = 0.5` returns the
  thresholds as before. Links stop on a table with another RP: the D link stays at RP50, as the GSED tau is an RP50 value.
* New `kits_age_pass(pass, age_m, wave_m, step, p)`: the age at which a share `p` of the children pass, for one or
  more `p` (a matrix, one column per `p`). `kits_age50()` is now a thin wrapper with `p = 0.5` and gives the same values.
  `kits_item_positions(items, norm, rp)` uses the column `age80` (and so on) for `rp = 0.8`, and
  `kits_plot_growth()` and `kits_plot_ladder()` label the plots with the RP.
* Quality-control notes (Table 8.8 p. 336, pp. 288, 301, 304, 333):
  * `kits_item()`: "smallest main booklet", a note when a booklet with 5% or more of the person-waves has fewer than
    400 (the rule of thumb for the Rasch model, p. 304).
  * Links between calibrations: "common items / items", a note when the kept common items are less than 20% of the
    items of either calibration (p. 288); `link$items_total` holds the item counts.
  * `kits_vertical()`: the cohorts at a shared wave also get `z_link` (the mean difference over the linking error,
    note beyond 2) and `sd_ratio` (note outside 0.8 to 1.2, p. 301).
  * `?kits_check` and the vignette have a QC checklist for new waves: drift check before fixing, fixing to the original
    calibration, double links, booklet n, common-item share, consistency with the other cohort, md5 provenance.

## Other changes

* `kits_vertical()`: when only one cohort at a shared wave has `se_link` (for example the reference cohort, which was
  not converted), the other now counts as 0 linking error in `z` and `z_link`; before, the linking error was left out
  of `z`. On D, both cohorts carry the error of the link to D, which cancels in their difference, so compare cohorts on
  the reference scale to judge the cohort link alone.
* `kits_equations()`: the `"item"` topic has the RP threshold; the `"link"` topic notes that the default screening is the
  rule of kitsPreview, and a link to anchors notes the RP50 of the GSED tau; "This link" gives the screening rule, the
  line before screening, the common-item share and the domain balance; `"vertical"` has `z_link` and `sd_ratio`.

# kits 0.13.0

The first batch of methods from Kolen and Brennan (2014), *Test equating, scaling, and linking* (3rd ed.).
Page numbers below are the printed pages of the book.

## New features

* Linking error at each theta and through a chain of links (pp. 248-258, 269-270, 276-277).
  * With `n_boot > 0`, a link keeps its bootstrap draws (`link$draws`: the a and b of the final screened line in each
    resample that could be computed) and their covariance (`link$ab_cov`). With `n_boot = 0` (the default) nothing changes.
  * New `kits_link_se(..., theta)`: the converted value and its linking error at chosen theta values (for example the
    wave medians), for one link or a chain. For one link it is the SD of a + b theta over the draws,
    sqrt(Var(a) + theta^2 Var(b) + 2 theta Cov(a, b)); at theta = 0 it is `a_error`.
    For a chain, the draws of different links are independent (different items) and every draw is paired with every
    draw of the next link; the variance is computed exactly from the moments of the draws, centred at the estimated lines
    (the Braun-Holland delta method of p. 269 plus the product term), so no random pairing and no seed.
  * `kits_convert()` adds `se_link` and `se_total = sqrt(se^2 + se_link^2)` to a `kits_subj()`, a `kits_subj_table()` or a
    thresholds table (with `thr_se`) when every link of the chain has draws. Otherwise no column is added and no warning
    is given. Converting in two steps gives the same `se_link` as one step.
    `se_link` is the error from sampling the linking items only.
  * `rbind()` of kits tables fills `se_link` and `se_total` with NA for the tables that do not have them.
* Slope against 1 (p. 229; Masters, 1984): a link between two calibrations made with `n_boot > 0` has the new diagnostic
  "slope vs 1", (b - 1) / b_error, with a note when it is beyond 2 (the scale stretches or shrinks between the two
  calibrations: the items discriminate differently). The free slope is kept; links to anchors skip the check.
* New `kits_vertical(s)`: the statistics of a vertical scale per cohort x wave (pp. 461-469): n, mean, SD, median,
  `sd_true = sqrt(max(SD^2 - mean(se^2), 0))`, the growth from the previous wave, the effect size of Eq. 9.47
  (Yen, 1986; the square root of the unweighted average of the two variances) and `p_sup`, the probability that a child of
  this wave scores above a child of the previous wave (based on order only, so unchanged by monotone transformations, p. 469).
  Cohorts at a shared wave (36 and 48 months in KIT) are compared (`attr(, "cohorts")`), with the linking error of each
  mean when the scores have `se_link`. The diagnostics note scale shrinkage (p. 468), means below the previous wave and
  cohorts that differ. `print()`, `summary()`, `plot()` (means +/- SD by age, blue and orange) and `kits_equations()`.
* New `kits_precision(item, forms, theta, link, criterion)`: the precision of each form across the scale
  (pp. 232, 405, 413, 463-465). The raw-score distribution P(r | theta) comes from the Lord-Wingersky recursion for
  polytomous items (Eq. 6.41), each raw score is mapped to the WLE of `kits_subj_table()`, and the result gives
  E[theta-hat | theta], the bias, the SD (conditional SEM) and the RMSE, optionally converted with links (for example to D),
  and the usable range of each form (RMSE < `criterion`, by default 0.5 logit). `forms` can be a `kits_subj_table()`
  (the standard form of each wave) or a named list of item sets. `plot()` draws the RMSE curves.
* New `kits_raw_equate(item, from_items, to_items, method = c("true", "observed"), theta_dist)`: raw-to-raw tables of two
  forms calibrated on one scale. True-score equating (pp. 192-197): tau_from(theta) = r, then tau_to(theta);
  0 goes to 0 and the maximum to the maximum (Eq. 6.24 with no guessing). Observed-score equating (pp. 197-200, 232):
  Lord-Wingersky distributions integrated over a theta distribution (a normal fitted to a theta sample, or given
  quadrature points), then the equipercentile method of chapter 2 (Eq. 2.14-2.18).
* `kits_equations()`: the `"link"` topic has the linking error at theta and through a chain, the slope against 1
  (links between calibrations) and, for links to anchors, a note on wording: linking to D is a calibration
  (scale aligning), not an equating (pp. 493-500). New topics `"vertical"`, `"precision"` and `"raw_equate"`, with methods
  for their results. A converted `kits_subj()` or `kits_subj_table()` with `se_link` says so. New references:
  Masters (1984), Thissen et al. (1995), Lord and Wingersky (1984), Yen (1986).
* `kits_code()` of a link made with `n_boot > 0` also computes the covariance of a and b, the linking error at theta
  (`se_link()`), and for links between calibrations the slope z. `kits_code()` does not cover the new
  `kits_vertical()`, `kits_precision()` and `kits_raw_equate()`.

## Other changes

* `?kits_linking` says that the link to D through anchor items aligns scales (a calibration), and how the
  slope against 1 and the linking error at theta are reported.

# kits 0.12.0

## Breaking changes

* `kits_data()` and its class are removed. The data go in directly, as in mirt and TAM,
  because the model is fixed (the PCM, which is dexter's ENORM with scores 0, ..., m).
  Every `dataSrc` takes a matrix (rows are person-waves with the person ids as row names, columns are items,
  NA means not administered), a wide data.frame of the same shape (person ids in a `person_id` column or the row names;
  every other column is an item; non-numeric item columns stop with `kits_error_input`),
  or the long data.frame with `person_id`, `item_id`, `item_score`.
  One converter checks them as `kits_data()` did, with the same classed errors (`kits_error_data`):
  unique person and item ids, whole scores of at least 0, no gaps in the observed scores, and a message
  (`kits_message_data`) for dropped items and person-waves without responses.
  The gap check applies where a model is fitted on the data (`kits_item()` without `fixed_params`, `kits()`,
  `kits_inter()`, `kits_dif()`); scoring with given parameters (`kits_subj()`, `kits_subj_table()`, `kits_itemfit()`)
  and `kits_drift()` (which reports a `"score gap"`) do not need every score to occur.
  Person information is now given as plain arguments. To migrate:
  * `d <- kits_data(X, persons, scale = "lan M03"); item <- kits_item(d)` becomes
    `item <- kits_item(X, scale = "lan M03")` (by default the scale is the name of the data, as before).
  * `kits_subj(d, item)` becomes `kits_subj(X, item, persons = persons)`: `persons` is a data.frame with `person_id`
    and any of `wave_m`, `age_m`, `sex`, `cohort`, joined on `person_id` (extra rows are fine; missing ones stop).
    The norm functions take the result as before (`kits_norm_screen(s)`, `kits_norms(s, ages)`, `kits_norm_z(s, norm)`,
    `kits_norm_ages(s)`), and the diagnostics report the reliability of each wave. `kits()` also takes `persons`.
  * `kits_subj_table(d, item)` and `kits_inter(d)` become `kits_subj_table(X, item, wave_m = )` and `kits_inter(X, wave_m = )`,
    with one value per row (or a vector named by person id). Without `wave_m`, all rows are now one group
    (`wave_m` is NA, and the fit of `kits_inter()` is named `"all"`); before, it was an error.
  * `kits_age50(d, step = 1)` becomes `kits_age50(X, age_m, wave_m, step = 1)`: numbers are scores and a score of at least
    `step` is a pass (so 0/1 codes work as before), logical values are used as they are.
    `age_m` and `wave_m` are required.
  * `kits_itemfit(d, item)`, `kits_drift(d, ref)` become `kits_itemfit(X, item)`, `kits_drift(X, ref)`.
* `kits_linking()` no longer computes the linking error by default: `n_boot = 0`, as in `kits_mean_sigma()`.
  Give `n_boot` (for example 200, or 2000 for a report) to get `a_error` and `b_error`.
  Without it the print shows no linking error and the diagnostics no longer note that it was not computed.

## New features

* `kits_dif()` and `kits_dif_combine()` (sex DIF) are now in kits, with the interface and output of kitsPreview,
  but computed without dexter: each sex is calibrated by the CML of kits, and the relative item parameters are compared
  as in `dexter::DIF()` (Bechger & Maris, 2015). `dif` agrees with kitsPreview to about 1e-15 on the test data.
  `kits_dif()` also keeps dexter's overall test (`attr(, "overall")`: chi-square, df, p), the groups, the items left out
  (items in one group only, with a message; items with different scores in the two groups, with a warning of the
  new class `kits_warning_data`), diagnostics and `kits_info()`.
  `kits_equations()` has a `"dif"` topic and a method; `kits_code()` writes a standalone RTMB script for a `kits_dif()`.
* The `DIF` column of `kits_dif_combine()` colours the items of `kits_plot_growth()` and `kits_plot_ladder()`
  (none gray, favors girls orange, favors boys blue, inconsistent black); the examples and the vignette show all four.
* `kits_info()` records a fingerprint of the person information (`persons`, `wave_m`, `sex`), and the form of a
  data.frame input (`"wide data.frame"` or `"long data.frame"`); `kits_code()` writes how to build `X` from it.

## Other changes

* The Wright map (`plot()` of a `kits`) has much less space between the persons panel and the items panel.
* `plot()` of a `kits_inter` with one fit needs no `wave`.

# kits 0.11.0

## Breaking changes

* `kits_linking_d()` is removed. `kits_linking()` now makes both kinds of link:
  `kits_linking(from, to)` with `to` a calibration links two calibrations (for example two cohorts),
  and with `to` a table of anchor items (`item_id`, `tau`) links to D.
  Replace `kits_linking_d(steps, anchors)` by `kits_linking(steps, anchors)`.
  The new arguments are `step` (all common steps, one step, or the pass step of the anchors)
  and `to_scale` (by default the scale of `to`, or `"D"` for anchors).
  Domains (a named list in `from`) now also work for links between calibrations.
  `link$settings` and `kits_info()` record `kind` and `step`; anchor points record their `step`.

## New features

* `kits_equations()`: the formulas behind a result, in Unicode or LaTeX, with references.
  They describe the kits engine (CML, WLE and the interaction model in RTMB, the shash norm with `optim()`),
  and end with the scales and numbers of the result, including the links applied by `kits_convert()`.
* `kits_code()` and its alias `kits_rtmb_code()`: a standalone RTMB and base R script that recomputes a result
  (`kits_item()`, `kits()`, `kits_inter()`, links, `kits_norms()`, `kits_norm_screen()`).
  The header gives the kits version, the scale and the data md5; the script checks the md5.
* `rbind()` of kits tables (for example two `kits_subj()` results) checks that their scales match.

## Other changes

* Tables longer than 30 rows print their first 10 rows (`print(x, n = )` for more).
* `kits_info()` records the name of a `kits_data()` or `kits_subj()` input (`input`),
  and the `a` and `b` of the links applied by `kits_convert()` (`links`).
* `kits_norm_ages()` also takes a `kits_subj()` result and uses its `age_m` and `wave_m` columns,
  like the other norm functions: `kits_norm_ages(s[screen$normed, ])`.
* The print of a link with domains labels the `domains` table as a domain check: one slope shared by the domains,
  fitted with domain intercepts on the kept points; it is not the slope `b` of the link. The numbers are unchanged.
  `kits_equations()` and `kits_code()` use the same label.
* `print()` of a `kits_drift()` result no longer writes `kits_item(X_new, ...)`: it names the recorded `kits_data()` input,
  or `<new data>`. `kits_drift()` now records its `kits_data()` input in `info$input`.

## Bug fixes

* `kits_norms()` no longer warns "stopped before converging (optim code 52)" when the fit has converged.
  With `factr = 1e-8` (the setting of cNORM), L-BFGS-B often ends at the optimum with code 52
  (its line search finds no lower value, because the requested precision is finer than the rounding of -log L).
  Code 52 now counts as converged when the largest absolute gradient from RTMB (over the coefficients not held at a limit)
  is at most 1e-6 x the number of person-waves and a restart from the solution does not lower -log L by 1e-6;
  otherwise it still warns (class `kits_warning_convergence`, now with the gradient in the message).
  The fitted norms are unchanged, so they still agree with kitsPreview and cNORM.
  `summary()` of a norm shows the new check "max |gradient|".
* `kits_equations()`:
  * Links cite the GSED D-score paper (Weber et al., 2019) only when they link to anchor items, not between calibrations.
  * A `kits_norm_screen()` result has its own topic, `"screen"`: the floor share of each wave (rounded to 1 decimal
    before it is compared with `cut`, as the code does), the cut, and the age range rule of `kits_norm_ages()`,
    without the shash norm and its references. The `"norms"` topic now covers only the norm fit,
    and adds the convergence rule. "This screen" lists `floor_pct` by wave.
  * The `"fit"` topic gives the formula of `pvalue` in `kits_itemfit()`.
  * Distinct symbols: the four polynomials of the norm have their own coefficients
    (c^(mu), c^(sigma), c^(epsilon), c^(delta); `c_mu`, ... in the scripts), and the free item parameters are eta
    (beta = beta0 + A eta), so `b` is only the booklet.
  * No double space in lists of waves ("waves 6, 12, 24"); the same fix in the diagnostics of `kits_norm_screen()`.
  * `as.character()` returns an unnamed character vector.
* `kits_code()` scripts no longer reuse names of the equations for other things:
  the threshold step keeps `beta` (all steps) and uses `beta_i` for one item, the WLE step counts item slots in `n_slot`
  (`J` stays the threshold Jacobian), and the step size, the exp matrix and local `mean`, `c`, `h`, `e` are renamed;
  `eta` holds the free item parameters and `h` the WLE curvature.
* The README pointed to `?kits` for the condition classes, which opens the help of `kits()`;
  it now points to `?"kits-package"`.
* `kits_info()` prints settings with several values without padding.
