#' kits: building the KIT developmental scale (kit scale)
#'
#' kits puts children of different ages, who answered different items, on one scale, and builds
#' same-age norms; KIT is Kids in Taiwan: National Longitudinal Study of Child Development and Care.
#' The workflow functions start with `kits_` (kit scale); under them, the `algo_` functions hold
#' the algorithms alone (numbers in, numbers out; the last item of the workflow below). The data go
#' in directly, as a response matrix, a wide or a long data.frame (see [kits_item()]); scores
#' start at 0.
#' The tutorials: `vignette("kits-flow", package = "kits")` (steps 0 to 6),
#' `vignette("kits-advanced", package = "kits")` (steps 7 to 10); the algorithm behind a result:
#' [kits_algorithm()], with the tutorials in `tutorial_algorithms/` and the teaching scripts in
#' `inst/scripts`.
#'
#' @section Main workflow:
#' * Step 0, data: the responses of each cohort x domain and a `persons` table, or one table of
#'   responses with the columns of the child, [kits_resp()]; [kits_sim()] simulates responses.
#' * Step 1, calibrate each cohort x domain: [kits_item()], the PCM (Masters, 1982) by CML in RTMB
#'   (Kristensen et al., 2016; Kristensen, 2026). Always give `scale`, e.g.
#'   `kits_item(X03, scale = "lan M03")`. Without links, [kits()] calibrates and scores in one step.
#' * Step 2, cohort link (M36 -> M03) through the common items: `kits_linking(item36, item03)`.
#' * Step 3, link to D (Weber et al., 2019) through the anchor items of [kits_anchors()]:
#'   `kits_linking(item03, anchors)`. A link in four steps: [kits_link_steps].
#' * Step 4, score and convert: [kits_subj()] (the WLE, Warm, 1989), then
#'   `s36 |> kits_convert(co) |> kits_convert(d_link)`; raw-score tables: [kits_subj_table()].
#' * Step 5, norms: [kits_norm_screen()], [kits_norm_ages()], [kits_norms()] (shash, Jones &
#'   Pewsey, 2009, with the model of Lenhard & Lenhard, 2021; in three steps: [kits_norm_steps])
#'   and the same-age z [kits_norm_z()].
#' * Step 6, checks: `summary()` and [kits_check()]; [kits_itemfit()], [kits_item_restscore()],
#'   the LR test of [kits_dif()], [kits_esf()], [kits_info()]; R's generics [kits_methods]; one-row summaries
#'   [as.data.frame.kits].
#' * Step 7, items and report plots: DIF [kits_dif()] (Bechger & Maris,
#'   2015) and [kits_dif_combine()], [kits_item_table()]; `plot()` of most results (plotly),
#'   [kits_save_plot()].
#' * Step 8, a new wave: [kits_drift()], then `kits_item(X_new, fixed_params = drift$fix)`.
#' * Step 9, the algorithms: [kits_algorithm()]; the calibration in three steps
#'   [kits_item_steps]; [kits_information()], [kits_expected_score()], [kits_pcm_prob()].
#' * Step 10, properties of the scale (Kolen & Brennan, 2014): `se_link` of [kits_convert()],
#'   [kits_vertical()], [kits_precision()], [kits_link_check()]; raw-score correspondence
#'   with [kits_subj_table()] and [kits_expected_score()].
#' * The algorithms alone (numbers in, numbers out, no scale labels), which the `kits_` functions
#'   call: [algo_pcm()], [algo_esf()], [algo_cml()], [algo_thresholds()], [algo_jme()],
#'   [algo_screen()], [algo_mean_sigma()], [algo_link_error()], [algo_line()], [algo_shash()],
#'   [algo_shash_z()], [algo_lw()]; for example `kits_stats(X) |> algo_cml() |> algo_thresholds()`.
#'
#' ```
#' # 0 data: two cohorts; items i3-i5 are common to both
#' set.seed(1)
#' betas <- list(i1 = -2, i2 = -1, i3 = c(-0.5, 0.5), i4 = 0, i5 = c(0.5, 1.5), i6 = 2, i7 = 3)
#' X03 <- kits_sim(rnorm(500, 0), betas[1:5], prefix = "a")
#' X36 <- kits_sim(rnorm(500, 2), betas[3:7], prefix = "b")
#' persons <- data.frame(visit_id = c(rownames(X03), rownames(X36)),
#'                       wave_m = rep(c(6, 36), each = 500))
#' persons$age_m <- persons$wave_m + runif(1000, -1, 1)
#' anchors <- data.frame(item_id = c("i1", "i2", "i3", "i4", "i5"), tau = c(46, 53, 56, 60, 63))
#' item03 <- kits_item(X03, scale = "lan M03")                                  # 1 calibrate
#' item36 <- kits_item(X36, scale = "lan M36")
#' co <- kits_linking(item36, item03)                                           # 2 cohort link
#' d_link <- kits_linking(item03, anchors)                                      # 3 link to D
#' s36 <- kits_subj(X36, item36, persons = persons) |>                         # 4 score
#'   kits_convert(co) |>                                                        #   M36 -> M03
#'   kits_convert(d_link)                                                       #   M03 -> D
#' kits_check(item03, item36, co, d_link, s36)                                  # 6 checks
#' ```
#'
#' @section Errors, warnings and messages:
#' Every error, warning and message of kits has a class, so it can be caught by class, for example
#' `tryCatch(kits_convert(s, link), kits_error_scale = function(e) ...)`.
#' All errors also have the class `kits_error`, all warnings `kits_warning`, all messages `kits_message`.
#' * `kits_error_input`: a wrong argument (type, length, missing column, unknown option).
#' * `kits_error_data`: a problem in the data (row names, ids, scores that are not whole numbers, score gaps,
#'   item patterns that form separate groups with no common items).
#' * `kits_error_items`: a problem in the item parameters or common items (missing parameters, too few common items or anchors).
#' * `kits_error_scale`: scales that do not match (a link applied to the wrong scale, a link between two equal scales,
#'   `fixed_params` on another scale than `scale`, a norm on another scale than the scores,
#'   tables of different scales stacked by `rbind()` or `bind_rows()`).
#' * `kits_error_fit`: a model that cannot be fitted (nothing to calibrate, too few points kept, no wave fitted).
#' * `kits_warning_convergence`: an optimizer that stopped before converging (the norms).
#' * `kits_warning_scale`: a conversion or link to check (scale not recorded, points negatively related,
#'   no linking error because all points belong to one item, a domain without kept anchors,
#'   a kits table stacked with a plain data.frame).
#' * `kits_warning_norms`: a norm to check (coefficients at their limits, a large delta).
#' * `kits_warning_data`: data left out of a comparison (items with different scores in the two groups of [kits_dif()]).
#' * `kits_warning_anchors`: anchors to use with care ([kits_anchors()] with `key = "gsed2510"`, another D metric).
#' * `kits_message_data`: data dropped or worth knowing (empty rows or items, items not drawn, items in one DIF group only).
#' * `kits_message_scale`: scores that are not on the D scale given to the norm functions, a scale labelled
#'   "scale 1" by default, links that share kept items, and a [kits_convert()] that adds no `se_link` because a link
#'   of the chain (in this call or an earlier one) has no bootstrap draws.
#' * `kits_message_fit`: items or waves left out of a fit.
#'
#' Messages and errors give counts, not person ids; set `options(kits.show_ids = TRUE)` to show the ids.
#'
#' To cite kits, see `citation("kits")`; please cite the methods below with it.
#'
#' @references
#' Bechger, T. M., & Maris, G. (2015). A statistical test for differential item pair functioning.
#' *Psychometrika, 80*(2), 317-340. \doi{10.1007/s11336-014-9408-y}
#'
#' Jones, M. C., & Pewsey, A. (2009). Sinh-arcsinh distributions. *Biometrika, 96*(4), 761-780.
#' \doi{10.1093/biomet/asp053}
#'
#' Kolen, M. J., & Brennan, R. L. (2014). *Test equating, scaling, and linking* (3rd ed.). Springer.
#' \doi{10.1007/978-1-4939-0317-7}
#'
#' Kristensen, K. (2026). RTMB: R bindings for TMB (R package version 2.0). \doi{10.32614/CRAN.package.RTMB}
#'
#' Kristensen, K., Nielsen, A., Berg, C. W., Skaug, H., & Bell, B. M. (2016). TMB: Automatic differentiation and
#' Laplace approximation. *Journal of Statistical Software, 70*(5), 1-21. \doi{10.18637/jss.v070.i05}
#'
#' Lenhard, W., & Lenhard, A. (2021). Improvement of norm score quality via regression-based continuous norming.
#' *Educational and Psychological Measurement, 81*(2), 229-261. \doi{10.1177/0013164420928457}
#'
#' Masters, G. N. (1982). A Rasch model for partial credit scoring. *Psychometrika, 47*(2), 149-174.
#' \doi{10.1007/BF02296272}
#'
#' Warm, T. A. (1989). Weighted likelihood estimation of ability in item response theory.
#' *Psychometrika, 54*(3), 427-450. \doi{10.1007/BF02294627}
#'
#' Weber, A. M., Rubio-Codina, M., Walker, S. P., van Buuren, S., Eekhout, I., et al. (2019).
#' The D-score: a metric for interpreting the early development of infants and toddlers across global settings.
#' *BMJ Global Health, 4*(6), e001724. \doi{10.1136/bmjgh-2019-001724}
#'
#' @keywords internal
"_PACKAGE"

## 整理資料用的 tidyverse 動詞
#' @importFrom dplyr filter mutate transmute summarise reframe arrange select rename rename_with relocate pull
#' @importFrom dplyr left_join inner_join anti_join bind_rows bind_cols distinct count n n_distinct row_number
#' @importFrom dplyr if_else case_when coalesce lag first desc across any_of all_of everything pick slice .data
#' @importFrom dplyr dplyr_reconstruct
#' @importFrom purrr map map_chr map_dbl map_int map_lgl map2 map2_dbl pmap pmap_chr imap walk walk2
#' @importFrom purrr list_rbind list_flatten compact discard
#' @importFrom tidyr pivot_wider expand_grid
NULL

# dplyr 的動詞裡用到的欄名（R CMD check 不把它們當成沒有定義的變數）
utils::globalVariables(c(
  "a", "n_misfit", "age_m", "age_med", "age_reach", "at_edge", "beta", "beta_se",
  "bin", "book",
  "pattern_id",
  "raw_max",
  "raw", "can_fit", "cohort", "cohort_1", "column", "D", "DIF", "dif", "dif_1",
  "dif_2",
  "differ", "dir_1", "dir_2", "domain", "drift", "dropped", "E", "empirical", "extreme", "first",
  "fit_row", "flag", "floor_pct", "form", "gap", "group", "growth", "has_step", "highest_p50", "in_both", "infit",
  "inner", "inside", "low", "item", "item_id", "item_max", "item_score", "items", "judged", "last", "lowest_p50",
  "max_abs_resid", "max_new", "max_ref", "max_score", "mean_1", "mean_diff", "mean_se2", "method", "n",
  "n_1", "n_above", "n_below", "n_items", "n_persons", "n_steps", "n_values", "n_wave", "no_gap", "normed",
  "ok", "outfit", "P50", "p_sup", "pass", "pct", "person", "visit_id", "pvalue",
  "reason",
  "rel",
  "resid", "rir", "row", "score", "scores", "scores_1", "scores_2", "sd", "sd_1", "se", "SE_beta",
  "se_diff", "se_link", "se_link_1", "se_method", "set", "sex", "share", "sigma", "status", "stem",
  "step", "suggest_fix", "tau", "theta", "thetas", "thr", "thr_from", "thr_new", "thr_ref", "thr_se",
  "thr_to", "top", "total", "type", "value", "var_theta", "W", "wave_m", "wave_n", "x", "y", "y_se", "z", "z_both", "z_link",
  "key", "max_raw", "miss", "raw", "reach", "text",
  "median", "P25", "P75", "check", "side", "domain", "wave",
  "free", "group_1", "group_2", "k", "gamma_obs", "gamma_exp", "p"
))
