# kits_algorithm()：一個結果（或一個主題）是哪個算法算的：幾步、每步一條式子、兩層的函數（kits_ 和 algo_）、
# 文獻、在哪裡讀更多（tutorial_algorithms/ 的教學檔、inst/scripts 的教學程式和它的段）。
# 資料驅動：一個主題一個 list（ALGO_TOPICS），一個 print 方法；結果另外加一行它自己的事實。
# 教學的說明和完整的程式在教學檔和 inst/scripts，不在這裡。原始碼只用 ASCII；Unicode 字元用 \u 寫。

# 文獻：一行，作者、年、期刊（或出版社）
ALGO_REFS <- c(
  masters1982 = "Masters, G. N. (1982). A Rasch model for partial credit scoring. Psychometrika, 47, 149-174.",
  tmb2016 = "Kristensen, K., Nielsen, A., Berg, C. W., Skaug, H., & Bell, B. M. (2016). J. Stat. Softw., 70(5).",
  adams2012 = "Adams, R. J., Wu, M. L., & Wilson, M. (2012). Educational and Psychological Measurement, 72, 547-573.",
  warm1989 = "Warm, T. A. (1989). Psychometrika, 54, 427-450.",
  firth1993 = "Firth, D. (1993). Bias reduction of maximum likelihood estimates. Biometrika, 80, 27-38.",
  magis2015 = "Magis, D. (2015). Psychometrika, 80, 200-204 (Jeffreys modal estimation = WLE for the PCM).",
  kolen2014 = "Kolen, M. J., & Brennan, R. L. (2014). Test equating, scaling, and linking (3rd ed.). Springer.",
  marco1977 = "Marco, G. L. (1977). Journal of Educational Measurement, 14, 139-160.",
  weber2019 = "Weber, A. M., Rubio-Codina, M., Walker, S. P., et al. (2019). BMJ Global Health, 4, e001724.",
  jones2009 = "Jones, M. C., & Pewsey, A. (2009). Sinh-arcsinh distributions. Biometrika, 96, 761-780.",
  lenhard2021 = "Lenhard, W., & Lenhard, A. (2021). Educational and Psychological Measurement, 81, 229-261.",
  andersen1973 = "Andersen, E. B. (1973). A goodness of fit test for the Rasch model. Psychometrika, 38, 123-140.",
  bechger2015 = "Bechger, T. M., & Maris, G. (2015). Psychometrika, 80, 317-340.",
  christensen2013 = "Christensen, K. B., Kreiner, S., & Mesbah, M. (Eds.) (2013). Rasch models in health. Wiley.",
  lord1984 = "Lord, F. M., & Wingersky, M. S. (1984). Applied Psychological Measurement, 8, 453-461.",
  thissen1995 = "Thissen, D., Pommerich, M., Billeaud, K., & Williams, V. S. L. (1995). Appl. Psychol. Meas., 19, 39-49.",
  kits = "kits: the rule of this package; see ?\"kits-package\"."
)

# 一步：說明和一條式子
algo_step <- function(text, formula) list(text = text, formula = formula)

# 每個主題：標題、步驟、兩層的函數（|> 的管線）、文獻（ALGO_REFS 的鍵）、教學程式和它的段
ALGO_TOPICS <- list(
  cml = list(
    title = "item calibration: the partial credit model (PCM) by conditional maximum likelihood (CML)",
    steps = list(
      algo_step("Sufficient statistics of each item pattern b (kits_stats())",
                "n_ig = #{p : x_pi \u2265 g},  n_br = #{p \u2208 b : r_p = r}"),
      algo_step("PCM probabilities of the steps (beta_ig = step parameter)",
                "P(X_pi = k | \u03b8_p) = exp(k\u03b8_p \u2212 \u03a3_{g\u2264k} \u03b2_ig) / \u03a3_h exp(h\u03b8_p \u2212 \u03a3_{g\u2264h} \u03b2_ig)"),
      algo_step("Elementary symmetric functions (ESF) of each item pattern",
                "\u03b3_r = [z^r] \u03a0_i (1 + \u03a3_k exp(\u2212\u03a3_{g\u2264k} \u03b2_ig) z^k)"),
      algo_step("Conditional log-likelihood; Newton with the gradient and Hessian of RTMB",
                "log L = \u2212\u03a3_ig \u03b2_ig n_ig \u2212 \u03a3_b \u03a3_r n_br log \u03b3_br,  \u03a3_ig \u03b2_ig = 0"),
      algo_step("Covariance and SE of the betas (eta = the free parameters, beta = beta_0 + A eta)",
                "V(\u03b2\u0302) = A H\u207b\u00b9 A\u1d40,  SE(\u03b2\u0302_ig) = \u221aV_ig,ig"),
      algo_step("Thurstonian thresholds (column thr) and their SE by the delta method",
                "P(X_i \u2265 k | \u03b8 = thr_ik) = RP (0.5),  SE(thr_ik) = \u221a(J\u1d40 V_i J)")),
    components = c("kits_stats(X) |> kits_cml() |> kits_thresholds()           # = kits_item(X)",
                   "kits_stats(X) |> algo_cml() |> algo_thresholds(rp = 0.5)",
                   "algo_pcm(theta, beta)                                      # step 2",
                   "algo_esf(betas)                                            # step 3",
                   "inside: cml_stats(), cml_gamma_matrix(), cml_nll(), cml_optimize()"),
    refs = c("masters1982", "tmb2016", "adams2012"),
    script = "cml.R", parts = c("stats", "cml", "thr", "fixed", "fixed_out")),
  wle = list(
    title = "person scores: the weighted likelihood estimate (WLE), its SE and person fit",
    steps = list(
      algo_step("Log-likelihood of person p from the PCM, over the items answered",
                "log L_p(\u03b8) = r_p \u03b8 \u2212 \u03a3_i log S_i(\u03b8),  S_i(\u03b8) = \u03a3_k exp(k\u03b8 \u2212 \u03a3_{g\u2264k} \u03b2_ig)"),
      algo_step("Test information of the items answered",
                "I_p(\u03b8) = \u03a3_i Var(X_i | \u03b8) = \u03a3_i \u03a3_k (k \u2212 E_i)\u00b2 P_ik"),
      algo_step("Jeffreys modal estimate, for the PCM exactly Warm's WLE; one RTMB Newton solve for all",
                "\u03b8\u0302_p = argmax_\u03b8 [ log L_p(\u03b8) + \u00bd log I_p(\u03b8) ]   (item pattern x raw score)"),
      algo_step("Asymptotic SE of the WLE (Warm, 1989), also for a zero or a perfect score",
                "SE(\u03b8\u0302_p) = 1 / \u221aI_p(\u03b8\u0302_p)"),
      algo_step("Person fit at the WLE (W = variance of x_pi)",
                "infit = \u03a3_i (x \u2212 E)\u00b2 / \u03a3_i W,  outfit = (1/n) \u03a3_i (x \u2212 E)\u00b2 / W;  misfit if max > 2"),
      algo_step("Conversion by kits_convert(), each link in turn",
                "\u03b8* = B + A \u03b8,  SE* = A SE,  se_total = \u221a(SE*\u00b2 + se_link\u00b2)")),
    components = c("kits_subj(X, item) |> kits_convert(link)",
                   "cml <- kits_stats(X) |> algo_cml()",
                   "algo_jme(cml$betas)                                        # theta and se per raw score",
                   "inside: wle_table(), wle_solve(), add_person_fit(), line_apply()"),
    refs = c("warm1989", "firth1993", "magis2015"),
    script = "wle.R", parts = c("wle", "person_fit", "convert", "warm")),
  link = list(
    title = "a link: mean/sigma with one round of screening, and the linking error",
    steps = list(
      algo_step("Points: common steps of two calibrations, or pass step s of the anchors",
                "(from_j, to_j) = (thr_ik, thr_ik)  or  (thr_is, \u03c4_i), \u03c4 = the GSED difficulty"),
      algo_step("mean/sigma line on all points (A = slope, B = intercept)",
                "A = SD(to) / SD(from),  B = mean(to) \u2212 A mean(from)"),
      algo_step("One round of screening, then the line again on the kept points",
                "e_j = to_j \u2212 (B0 + A0 from_j);  drop j if |e_j| > cut SD(e)"),
      algo_step("Linking error: R resamples of whole items, the screening redone in each",
                "SE(B) = SD_r(B_r),  SE(A) = SD_r(A_r),  Cov(A, B) over the resamples"),
      algo_step("Linking error at theta (a chain of links adds up the same way)",
                "SE_link(\u03b8) = \u221a(Var(B) + \u03b8\u00b2 Var(A) + 2\u03b8 Cov(A, B))"),
      algo_step("Domain check: one slope shared by the domains d, one intercept each",
                "A_s = \u221a(\u03a3 (to \u2212 mean_d to)\u00b2 / \u03a3 (from \u2212 mean_d from)\u00b2),  B_d = mean_d to \u2212 A_s mean_d from")),
    components = c("pts <- kits_link_points(from, to)                          # from, to: two kits_item()",
                   "pts |> kits_link_screen(cut = 2) |> kits_link_line() |> kits_link_error(n_boot = 200)",
                   "screen <- algo_screen(pts$from, pts$to, pts$item, cut = 2)",
                   "screen |> algo_mean_sigma()                                # the line: B, A",
                   "screen |> algo_link_error(n_boot = 200)                    # SE(B), SE(A)",
                   "algo_line(link, theta)                                     # new = B + A x old"),
    refs = c("marco1977", "kolen2014", "weber2019"),
    script = "link.R", parts = c("points", "line", "error", "domains", "slope", "rules")),
  norms = list(
    title = "a continuous norm: the sinh-arcsinh (shash) distribution by maximum likelihood",
    steps = list(
      algo_step("Age, transformed and standardised (age_x = sqrt by default)",
                "x = (age_x(age) \u2212 m) / s"),
      algo_step("Four polynomials in x (degrees = their top powers)",
                "\u03bc(x) = \u03a3_j c_j^(\u03bc) x^j;  log \u03c3(x), \u03b5(x), log \u03b4(x) the same way"),
      algo_step("shash density; -log L minimised by nlminb with the gradient and Hessian of RTMB",
                "log f = log \u03b4 \u2212 log \u03c3 + log cosh w \u2212 \u00bd log(2\u03c0(1 + u\u00b2)) \u2212 \u00bd sinh\u00b2 w,  w = \u03b4 asinh u + \u03b5"),
      algo_step("Bounds as in cNORM, and several starts of delta; u = (D - mu) / sigma",
                "log \u03c3, \u03b5 \u2208 [\u221210, 10],  log \u03b4 \u2208 [\u22122, 2];  one fit per \u03b4_0, the best kept"),
      algo_step("Same-age z and percentile rank",
                "z = sinh(\u03b4 asinh((D \u2212 \u03bc) / \u03c3) + \u03b5),  PR = 100 \u03a6(z)"),
      algo_step("Percentile curves (columns P3, P5, ...)",
                "P_q = \u03bc + \u03c3 sinh((asinh(\u03a6\u207b\u00b9(q / 100)) \u2212 \u03b5) / \u03b4)")),
    components = c("kits_norm_data(D, age_m, ages) |> kits_norm_fit() |> kits_norm_table()  # = kits_norms()",
                   "kits_norm_z(D, age_m, norm)                                # z from the norm table",
                   "fit <- algo_shash(sqrt(age_m), D)                          # x = transformed age",
                   "fit |> algo_shash_z(sqrt(age_m), D)                        # z at the exact age",
                   "inside: best_shash_fit(), shash_quantile(), shash_z()"),
    refs = c("jones2009", "lenhard2021", "tmb2016"),
    script = "norms.R", parts = c("data", "fit", "table", "z", "by")),
  norm_screen = list(
    title = "screening the waves of a norm (scores piled at the floor) and the ages of the norm",
    steps = list(
      algo_step("Floor share of each cohort x wave w (D to 3 decimals)",
                "floor_pct_w = round(100 \u00d7 #{p \u2208 w : D_p = min_{q\u2208w} D_q} / n_w, 1)"),
      algo_step("Leave out the waves piled at the floor (default cut = 20)",
                "wave w left out of the norms if floor_pct_w \u2265 cut"),
      algo_step("Age range of the norms, in steps of by (0.5 month)",
                "ages = by \u230aQ_0.01(age | first wave) / by\u230b, ..., by \u2308Q_0.99(age | last wave) / by\u2309")),
    components = c("screen <- D |> kits_norm_screen(wave_m, cohort)",
                   "kits_norm_ages(age_m[screen$normed], wave_m[screen$normed])   # then kits_norms()",
                   "inside: floor_shares()"),
    refs = c("lenhard2021", "kits"),
    script = "norm_screen.R", parts = c("screen", "ages")),
  dif = list(
    title = "DIF: the CML in each group, the overall LR test, Wald z per step, item DIF",
    steps = list(
      algo_step("Groups: group = a vector (e.g. sex), \"raw\" (Andersen's median split) or \"wave\"",
                "\"raw\": low if r_p \u2264 median of r in the item pattern of p (ties go low)"),
      algo_step("CML of the PCM in each group g, and of all groups together",
                "log L_g = \u03a3_{p\u2208g} log P(x_p | r_p; \u03b2_g)"),
      algo_step("Overall likelihood-ratio test (S = steps, c_g = unlinked item sets)",
                "LR = 2 (\u03a3_g log L_g \u2212 log L_all),  df = \u03a3_g (S_g \u2212 c_g) \u2212 (S \u2212 1)"),
      algo_step("Wald z of each step (group betas shifted to the joint mean)",
                "d_s = \u03b2\u0302_s,h \u2212 \u03b2\u0302_s,g,  z_s = d_s / \u221a(SE\u00b2_s,g + SE\u00b2_s,h)"),
      algo_step("Item DIF with two groups: item pairs, r = ref, f = the other group",
                "dif_ig = (\u03b2_ig^(r) \u2212 \u03b2\u0304^(r)) \u2212 (\u03b2_ig^(f) \u2212 \u03b2\u0304^(f)),  dif_i = mean_g dif_ig")),
    components = c("d <- kits_dif(X, group = sex)             # or group = \"raw\", \"wave\"",
                   "d                                         # print: groups, the LR test, the largest step differences",
                   "split(as.data.frame(X), sex) |> lapply(\\(x) kits_stats(x) |> algo_cml())",
                   "kits_dif_combine(dif_cohort1, dif_cohort2)  # two cohorts"),
    refs = c("andersen1973", "bechger2015", "christensen2013"),
    script = "dif.R", parts = c("cml", "dif", "overall")),
  drift = list(
    title = "item drift: a new wave against the reference calibration",
    steps = list(
      algo_step("Free CML of the new wave, as kits_algorithm(\"cml\")",
                "thr_new,ik = the thresholds of kits_item(X_new)"),
      algo_step("Steps that can be judged (min_n = 30)",
                "step k of item i judged if min(n_below,ik, n_above,ik) \u2265 min_n"),
      algo_step("mean/sigma line new -> ref with one round of screening (cut = 2)",
                "thr_ref,ik = B + A thr_new,ik;  drop if |e_ik| > cut SD(e)"),
      algo_step("Drifted items",
                "drift_i = a judged step of item i dropped by the screening"),
      algo_step("Items to fix in the next calibration (drift$fix)",
                paste("suggest_fix_i = (judged and not drift_i), too little information",
                      "or a higher top score than in ref"))),
    components = c("kits_item(X_new) |> kits_linking(ref)        # inside kits_drift(X_new, ref)",
                   "kits_stats(X_new) |> algo_cml() |> algo_thresholds()",
                   "algo_screen(thr_new, thr_ref, item_id, cut = 2) |> algo_mean_sigma()",
                   "kits_item(X_new, fixed_params = drift$fix)   # the next calibration"),
    refs = c("kolen2014", "marco1977", "kits"),
    script = "drift.R", parts = c("points", "line", "items")),
  information = list(
    title = "information, expected score and precision of a form",
    steps = list(
      algo_step("PCM probabilities of each item",
                "P_ik(\u03b8) = exp(k\u03b8 \u2212 \u03a3_{g\u2264k} \u03b2_ig) / \u03a3_h exp(h\u03b8 \u2212 \u03a3_{g\u2264h} \u03b2_ig)"),
      algo_step("Test information (kits_information())",
                "I(\u03b8) = \u03a3_i \u03a3_k (k \u2212 E_i)\u00b2 P_ik,  E_i = \u03a3_k k P_ik"),
      algo_step("Expected raw score (kits_expected_score()); RTMB checks its slope",
                "E(\u03b8) = \u03a3_i E_i(\u03b8),  dE/d\u03b8 = I(\u03b8)"),
      algo_step("Raw-score distribution of a form (Lord-Wingersky)",
                "f_j(x | \u03b8) = \u03a3_k f_{j\u22121}(x \u2212 k | \u03b8) P_jk(\u03b8)"),
      algo_step("Bias and RMSE of the WLE (kits_precision()); usable where RMSE < criterion",
                "bias = \u03a3_r P(r | \u03b8) \u03b8\u0302(r) \u2212 \u03b8,  RMSE = \u221a(bias\u00b2 + SD\u00b2)")),
    components = c("kits_information(item)(theta), kits_expected_score(item)(theta)",
                   "kits_precision(item, forms)",
                   "algo_cml(kits_stats(X))$betas |> algo_lw(theta = theta)    # raw-score distribution",
                   "inside: pcm_function(), tcc(), form_precision()"),
    refs = c("lord1984", "thissen1995", "kolen2014", "warm1989"),
    script = "information.R", parts = c("information", "expected_score"))
)

# 每個主題的教學檔（tutorial_algorithms/）
ALGO_TUTORIALS <- c(cml = "01_calibration", wle = "02_scoring", link = "03_linking", norms = "04_norms",
                    norm_screen = "04_norms", dif = "05_checks", drift = "05_checks", information = "06_scale")

# 結果的類別 -> 主題
ALGO_CLASSES <- c(kits = "cml", kits_item = "cml", kits_stats = "cml", kits_subj = "wle",
                  kits_link = "link", kits_norms = "norms", kits_norm_screen = "norm_screen",
                  kits_dif = "dif", kits_lrtest = "dif", kits_drift = "drift",
                  kits_precision = "information", kits_information = "information",
                  kits_expected_score = "information")

# 結果自己的事實：一行
algo_facts <- function(x) {
  if (inherits(x, "kits")) x <- x$item
  if (inherits(x, "kits_item")) {
    n_pat <- length(unique(x$fit$inputs$design$pattern_id))
    return(sprintf("%d items, %d steps, %d person-waves, %d item pattern(s); scale %s",
                   length(unique(x$thresholds$item_id)), nrow(x$thresholds), x$n_persons, n_pat,
                   scale_label(scale_of(x))))
  }
  if (inherits(x, "kits_stats")) {
    return(sprintf("%d person-waves, %d items, %d steps, %d item pattern(s)", x$n_persons,
                   length(x$items), nrow(x$steps), nrow(x$patterns)))
  }
  if (inherits(x, "kits_subj")) {
    links <- attr(x, "info")$links
    conv <- if (is.null(links)) "not converted" else
      paste("converted", paste(c(links$from[1], links$to), collapse = " -> "))
    return(sprintf("%d person-waves; scale %s (%s)", nrow(x), scale_label(scale_of(x)), conv))
  }
  if (inherits(x, "kits_link")) {
    boot <- if (link_has_draws(x)) sprintf("; linking error over %d resamples", nrow(x$draws)) else ""
    return(sprintf("%s -> %s: %d points, %d dropped (%s); intercept %s, slope %s%s",
                   scale_label(x$from_scale), scale_label(x$to_scale), x$n, length(x$dropped),
                   screen_rule_text(x$settings), fmt(x$intercept, 2), fmt(x$slope, 2), boot))
  }
  if (inherits(x, "kits_norms")) {
    a <- attr(x, "kits_args")
    deg <- if (is.null(a$degrees)) "" else paste0(", degrees ", paste(a$degrees, collapse = ", "))
    return(sprintf("ages %s to %s months, %d person-waves%s; scale %s", format(min(x$age_m)),
                   format(max(x$age_m)), x$n[1], deg, scale_label(scale_of(x))))
  }
  if (inherits(x, "kits_norm_screen")) {
    w <- x$waves
    label <- if (is.null(w$cohort)) paste(w$wave_m) else paste(w$cohort, w$wave_m)
    left <- if (all(w$normed)) "none" else paste(label[!w$normed], collapse = ", ")
    return(sprintf("%d waves, left out (cut %s): %s; %d of %d person-waves in the norms",
                   nrow(w), format(x$cut), left, sum(x$normed), length(x$normed)))
  }
  if (inherits(x, c("kits_dif", "kits_lrtest"))) {   # DIF 的表（屬性 groups、overall），或 LR 檢定的 list
    g <- if (is.data.frame(x)) attr(x, "groups") else x$groups
    ov <- if (is.data.frame(x)) attr(x, "overall") else list(stat = x$LR, df = x$df, p = x$p)
    groups <- if (is.data.frame(g)) paste(g$group, collapse = ", ") else sprintf("%s (ref), %s", g$ref, g$focal)
    test <- if (is.null(ov$stat)) "" else
      sprintf("; overall test %s, df %d, p %s", fmt(ov$stat, 2), as.integer(ov$df), format.pval(ov$p, digits = 3))
    items <- if (is.data.frame(x)) sprintf("%d items; ", nrow(x)) else ""
    return(sprintf("%sgroups %s%s", items, groups, test))
  }
  if (inherits(x, "kits_drift")) {
    it <- as.data.frame(x$items)
    drifted <- it$item_id[which(it$drift)]
    return(sprintf("%d common items, %d judged, drifted: %s; %d to fix; scale of ref %s", nrow(it),
                   sum(!is.na(it$drift)), if (length(drifted) > 0) paste(drifted, collapse = ", ") else "none",
                   sum(it$suggest_fix), scale_label(x$scale)))
  }
  if (inherits(x, "kits_precision")) {
    r <- attr(x, "ranges")
    return(sprintf("%s; RMSE < %s; scale %s", paste(sprintf("form %s (%d items): usable %s to %s", r$form,
                   r$n_items, fmt(r$from, 2), fmt(r$to, 2)), collapse = ", "),
                   format(signif(attr(x, "criterion"), 3)), scale_label(scale_of(x))))
  }
  sprintf("%d item(s); theta on scale %s", length(attr(x, "items")), scale_label(attr(x, "scale")))
}

#' The algorithm behind a result
#'
#' Says, in about 30 lines, which algorithm produced a result: its steps, one formula per step, the
#' functions that run the steps (the `kits_` steps and the `algo_` building blocks, as `|>`
#' pipelines), key references, and where to read more (the tutorial page in
#' `tutorial_algorithms/` and the teaching script in `inst/scripts`, with its parts).
#'
#' @param x A result: [kits_item()] or [kits()], the scores of [kits_subj()] or [kits_convert()],
#'   [kits_linking()], [kits_norms()], [kits_norm_screen()], [kits_dif()], [kits_drift()],
#'   [kits_precision()] or [kits_information()]. Or a topic: `"cml"`, `"wle"`, `"link"`,
#'   `"norms"`, `"norm_screen"`, `"dif"`, `"drift"`, `"information"`.
#' @return Invisibly, a `kits_algorithm` object (a list: `topic`, `title`, `facts`, `steps`,
#'   `components`, `references`, `tutorial`, `script`, `parts`); it prints the algorithm. For a
#'   result, `facts` is one line on the result itself (items, steps, persons, scale).
#' @references
#' Masters, G. N. (1982). A Rasch model for partial credit scoring. Psychometrika, 47, 149-174.
#'
#' Warm, T. A. (1989). Weighted likelihood estimation of ability in item response theory.
#' Psychometrika, 54, 427-450.
#'
#' Kolen, M. J., & Brennan, R. L. (2014). Test equating, scaling, and linking (3rd ed.). Springer.
#' @examples
#' kits_algorithm("cml")
#'
#' betas <- list(a1 = c(-1, 0), a2 = 0.5, a3 = -0.5)
#' item <- kits_sim(rnorm(200), betas, seed = 1) |>
#'   kits_item(scale = "demo")
#' kits_algorithm(item)
#'
#' # The teaching script and its parts
#' algo <- kits_algorithm("link")
#' system.file("scripts", algo$script, package = "kits")
#' algo$parts
#' @export
kits_algorithm <- function(x) {
  if (is.character(x)) {
    abort_if(length(x) != 1 || !x %in% names(ALGO_TOPICS),
             "`x` must be a kits result or one topic: ", paste(names(ALGO_TOPICS), collapse = ", "), ".")
    topic <- x
    facts <- NULL
  } else {
    hit <- names(ALGO_CLASSES)[map_lgl(names(ALGO_CLASSES), \(cl) inherits(x, cl))]
    abort_if(length(hit) == 0,
             "kits_algorithm() covers the main estimates: ", paste(names(ALGO_CLASSES), collapse = ", "),
             "; or give a topic: ", paste(names(ALGO_TOPICS), collapse = ", "), ".")
    topic <- ALGO_CLASSES[[hit[1]]]
    facts <- algo_facts(x)
  }
  a <- ALGO_TOPICS[[topic]]
  out <- list(topic = topic, title = a$title, facts = facts,
              steps = data.frame(step = map_chr(a$steps, "text"), formula = map_chr(a$steps, "formula")),
              components = a$components, references = unname(ALGO_REFS[a$refs]),
              tutorial = sprintf("tutorial_algorithms/%s.html", ALGO_TUTORIALS[[topic]]), script = a$script,
              parts = a$parts)
  print(structure(out, class = "kits_algorithm"))
}

#' @export
print.kits_algorithm <- function(x, ...) {
  cat(sprintf("<kits_algorithm> %s: %s\n", x$topic, x$title))
  if (!is.null(x$facts)) cat("This result: ", x$facts, "\n", sep = "")
  cat("Steps:\n")
  cat(sprintf("  %d. %s\n       %s\n", seq_len(nrow(x$steps)), x$steps$step, x$steps$formula), sep = "")
  cat("Components:\n", paste0("  ", x$components, "\n"), sep = "")
  cat("References:\n", paste0("  ", x$references, "\n"), sep = "")
  cat("Read more:\n  ", x$tutorial, "\n", sep = "")
  cat(sprintf("  system.file(\"scripts\", \"%s\", package = \"kits\"), parts: %s\n", x$script,
              paste(x$parts, collapse = ", ")))
  invisible(x)
}
