# 連結的檢查（Kolen & Brennan, 2014）：篩點把換算值移動多少（summary(link, theta = )），
# 和兩個不用共同題篩點的獨立方法比較世代連結（kits_link_check()）：同齡的兩群孩子配平均和 SD、合併校準（concurrent CML）。

# 篩點把換算值移動多少（內部用，給 summary(link, theta = )）：
# 篩點前的線（intercept0、slope0）和最後的線（intercept、slope）在 theta 的值和差；有重抽時也除以 theta 的連結誤差
# （Kolen & Brennan, 2014, pp. 190, 202-203：比較有和沒有離群題的結果）
screening_shift <- function(link, theta) {
  abort_if(!inherits(link, "kits_link"), "`link` must be a link (kits_linking(), kits_link_line()).")
  need_theta(theta, ", on the scale where the link starts")
  theta <- as.vector(theta)
  first <- list(intercept = link$intercept0, slope = link$slope0)
  out <- data.frame(theta = theta, before = line_apply(first, theta),
                    after = line_apply(link, theta))
  out$shift <- out$after - out$before
  if (link_has_draws(link)) {
    out$se_link <- sqrt(link_moments(list(link), theta)$var)
    out$z <- out$shift / out$se_link
  }
  out <- kits_table(out)
  attr(out, "scale") <- if (is.null(link$to_scale)) NA_character_ else link$to_scale
  out
}


#' Independent checks of a link between two calibrations
#'
#' Estimates a cohort link in two ways that do not use the screened common items, and compares
#' the converted values at chosen theta with the main link (Kolen & Brennan, 2014):
#' same-age groups at the shared waves (a latent normal per cohort and wave, item parameters
#' fixed), and a concurrent CML calibration of the stacked data of both cohorts.
#' Details: the models tutorial (tutorial_models/kits_models.html), section 7, and
#' `vignette("kits-advanced", package = "kits")`, section 10.7.
#'
#' @param link The link to check: a [kits_linking()] between two calibrations (not to anchors, not with domains).
#' @param item_from,item_to The two separate calibrations ([kits_item()]) that the link joins,
#'   on the scales `link$from_scale` and `link$to_scale`.
#' @param X_from,X_to The responses of each calibration (a matrix or data.frame, as for [kits_item()]).
#'   The person ids of the two may overlap; they are kept apart in the concurrent calibration.
#' @param persons Optional person information with `visit_id` and `wave_m` for both cohorts (ids
#'   must not repeat). Without it, or without shared waves, the same-age check is skipped.
#' @param persons_from,persons_to The person information of each cohort on its own (`visit_id`,
#'   `wave_m`), for when the person ids overlap. Give either `persons` or both of these.
#' @param theta Values on the `from` scale at which to compare the lines. By default the median theta of each wave of
#'   `from` (or the quartiles without `persons`).
#' @param tol Tolerance for the difference when the link has no bootstrap draws (units of the target scale; default 0.2).
#' @return A `kits_link_check`: a list with
#'   * `methods`: one row per method: `method`, `intercept`, `slope`, `basis` (what it rests on).
#'   * `at`: one row per method and theta: `theta`, `method`, `value`, `diff`, `se_method`, and
#'     with draws `se_link`, `z` and `z_both`.
#'   * `waves`: one row per shared wave (NULL when none): `n`, latent `mean`, `sd` and their SEs of
#'     each cohort (`_from`, `_to`), and `converged`.
#'   * `concurrent`: the concurrent [kits_item()]; `lines`: the lines of the separate calibrations on its scale.
#'   * `from_scale`, `to_scale`, `settings` (`tol`), `diagnostics` (see [kits_check()]) and `info` (see [kits_info()]).
#' @references
#' Kolen, M. J., & Brennan, R. L. (2014). *Test equating, scaling, and linking* (3rd ed.). Springer.
#' \doi{10.1007/978-1-4939-0317-7}
#' @examples
#' # Two cohorts measured at 36 months; items a4 to a7 are common
#' betas <- lapply(seq(-3, 3, length.out = 10), function(b) b + c(-0.5, 0.5))
#' names(betas) <- paste0("a", 1:10)
#' set.seed(1)
#' wave <- rep(c(24, 36, 36, 48), each = 400)
#' theta <- 2 * log(wave / 24) + rnorm(1600)
#' young <- kits_sim(theta[1:800], betas[1:7], prefix = "y", seed = 1)
#' old <- kits_sim(theta[801:1600], betas[4:10], prefix = "o", seed = 2)
#' persons <- data.frame(visit_id = c(rownames(young), rownames(old)), wave_m = wave)
#' item_young <- kits_item(young, scale = "young")
#' item_old <- kits_item(old, scale = "old")
#' link <- kits_linking(item_old, item_young, n_boot = 100)
#' chk <- kits_link_check(link, item_old, item_young, old, young, persons = persons)
#' chk
#' @export
kits_link_check <- function(link, item_from, item_to, X_from, X_to, persons = NULL, theta = NULL, tol = 0.2,
                            persons_from = NULL, persons_to = NULL) {
  call <- call_text(sys.call())
  abort_if(!inherits(link, "kits_link") || identical(link$settings$kind, "anchors") || !is.null(link$domains),
           "`link` must be a kits_linking() between two calibrations (not to anchors, not with domains).")
  abort_if(!inherits(item_from, "kits_item") || !inherits(item_to, "kits_item"),
           "`item_from` and `item_to` must be the kits_item() results that the link joins.")
  abort_if(!identical(scale_of(item_from), link$from_scale[1]) || !identical(scale_of(item_to), link$to_scale),
           sprintf("The link converts %s -> %s, but the calibrations are on %s and %s.", scale_label(link$from_scale),
                   scale_label(link$to_scale), scale_label(scale_of(item_from)), scale_label(scale_of(item_to))),
           class = "kits_error_scale")
  need_number(tol, "tol", "one positive number (units of the target scale)", above = 0)
  resp_from <- resp_long(X_from, "X_from", gaps = FALSE)
  resp_to <- resp_long(X_to, "X_to", gaps = FALSE)
  cohorts <- cohort_persons(persons, persons_from, persons_to)

  # 同齡的兩群孩子：各自在自己的校準上，估這一波的能力分佈（常態），配平均和 SD
  s_from <- kits_subj(resp_mark(resp_from), item_from, persons = cohorts$from)
  same_age <- NULL
  if (!is.null(cohorts$to) && "wave_m" %in% names(s_from)) {
    s_to <- kits_subj(resp_mark(resp_to), item_to, persons = cohorts$to)
    same_age <- same_age_lines(resp_from, resp_to, item_from, item_to,
                               stats::setNames(s_from$wave_m, s_from$visit_id),
                               stats::setNames(s_to$wave_m, s_to$visit_id))
  }
  if (is.null(theta)) theta <- default_check_theta(s_from)
  need_theta(theta, ", on the scale where the link starts")
  theta <- as.vector(theta)

  # 合併校準：兩邊的 visit_id 分開（加上前綴），所有題一起用 CML 校準
  conc <- concurrent_lines(resp_from, resp_to, item_from, item_to)

  methods <- data.frame(method = "common items (the link)", intercept = link$intercept,
                        slope = link$slope,
                        basis = sprintf("%d common steps of %d items, %d dropped (%s)", link$n, link$n_items,
                                        length(link$dropped),
                                        screen_rule_text(link$settings)))
  methods <- bind_rows(methods, same_age$methods, conc$method) |>
    kits_table()

  at <- pmap(list(methods$method, methods$intercept, methods$slope), \(method, intercept, slope) {
    j <- match(method, same_age$methods$method)
    line <- list(intercept = intercept, slope = slope)
    gap <- line_minus(line, link)   # 兩條線的差也是一條線
    data.frame(theta = theta, method = method, value = line_apply(line, theta),
               diff = line_apply(gap, theta),
               se_method = if (is.na(j)) NA_real_ else same_age_se(same_age$waves, same_age$rows[[j]], theta))
  }) |>
    list_rbind()
  if (link_has_draws(link)) {
    at <- at |>
      mutate(se_link = sqrt(link_moments(list(link), theta)$var),
             z = diff / se_link,
             z_both = diff / sqrt(se_link^2 + if_else(is.na(se_method), 0, se_method^2)))
  }
  at <- kits_table(at)
  attr(at, "scale") <- link$to_scale
  out <- structure(list(methods = methods, at = at, waves = same_age$waves, concurrent = conc$item, lines = conc$lines,
                        from_scale = link$from_scale, to_scale = link$to_scale, settings = list(tol = tol)),
                   class = "kits_link_check")
  out$diagnostics <- link_check_diagnostics(out, link_has_draws(link), !is.null(cohorts$from))
  out$info <- make_info(call, list(tol = tol, theta = theta),
                        conc$item$info$data)   # 兩邊的作答疊起來（visit_id 加上 from、to）
  out
}

# 一個方法和連結差多少，以連結誤差為單位（同齡的線再加上它自己的抽樣誤差）；超過 2 註記
method_z_row <- function(m, d) {
  own <- any(is.finite(d$se_method))
  z_values <- if (own) d$z_both else d$z
  z <- max(abs(z_values), na.rm = TRUE)
  status <- if (z > 2) "note" else "ok"
  units <- if (own) "SEs (the linking error and the sampling error of the groups combined)" else "linking errors"
  detail <- if (z > 2) " (beyond 2: the methods disagree more than these errors explain)" else ""
  diag_row(m, z, status, sprintf("largest |difference| from the link: %.3g, %.2f %s%s", max(abs(d$diff)), z, units, detail))
}

# 沒有重抽時：一個方法和連結差多少，和 tol 比
method_tol_row <- function(m, d, tol) {
  big <- max(abs(d$diff))
  status <- if (big > tol) "note" else "ok"
  message <- sprintf(paste0("largest |difference| from the link: %.3g (tolerance %s; give the link n_boot > 0 to compare with ",
                            "its linking error)"),
                     big, format(tol))
  diag_row(m, big, status, message)
}

# 兩個世代的人的資料：persons（兩個世代的 visit_id 不重複時）或 persons_from 和 persons_to（重複時）
cohort_persons <- function(persons, persons_from, persons_to) {
  persons <- old_ids(persons, "persons")
  persons_from <- old_ids(persons_from, "persons_from")
  persons_to <- old_ids(persons_to, "persons_to")
  by_cohort <- !is.null(persons_from) || !is.null(persons_to)
  abort_if(by_cohort && !is.null(persons),
           "Give either `persons` (person ids that do not overlap) or `persons_from` and `persons_to`, not both.")
  abort_if(by_cohort && (is.null(persons_from) || is.null(persons_to)),
           "Give both `persons_from` and `persons_to`, one table per cohort.")
  if (by_cohort) {
    return(list(from = persons_from, to = persons_to))
  }
  repeated <- is.data.frame(persons) && "visit_id" %in% names(persons) &&
    anyDuplicated(persons$visit_id) > 0
  abort_if(repeated,
           "`persons` has repeated person ids. When the person ids of the two cohorts overlap, give one table per ",
           "cohort: `persons_from` and `persons_to`.")
  list(from = persons, to = persons)
}

# 預設的 theta：from 每個波次的中位數；沒有波次時用四分位數
default_check_theta <- function(s) {
  if ("wave_m" %in% names(s)) {
    medians <- s |>
      summarise(median = stats::median(theta), .by = wave_m) |>
      arrange(wave_m)
    return(unique(medians$median))
  }
  unname(stats::quantile(s$theta, c(0.25, 0.5, 0.75)))
}

# 一群孩子（一個世代的一個波次）的能力分佈 N(mu, sigma^2)：題目參數固定在自己的校準，邊際最大概似（RTMB）。
# PCM 的原始分是充分統計量，所以每個人的邊際概似只和「題目組合 x 原始分」有關：
#   L_p = integral exp(r_p theta - sum_{i in b_p} log S_i(theta)) phi(theta; mu, sigma) d theta（作答型態的常數項不影響 mu、sigma）
# log S_i(theta) 是 PCM 的核心 pcm_kernel()（R/pcm_core.R）算的。
# 積分用固定的標準常態格點 z（-6 到 6，121 點），theta = mu + sigma z；log sum exp 用 RTMB 的 logspace_add()。
# 只寫 -log L：梯度、Hessian 由 RTMB 給（fit_ad()），sd = exp(log sd) 的 SE 由 sdreport() 給。
latent_normal <- function(resp, parms) {
  betas <- betas_by_item(parms)
  persons <- person_sets(resp)
  combos <- summarise(persons, n = n(), .by = c(set, raw))   # 每個「題目組合 x 原始分」一列，n 是人數
  items_of <- map(combos$set, set_items)
  r <- combos$raw
  n <- combos$n
  z <- seq(-6, 6, length.out = 121)
  log_w <- stats::dnorm(z, log = TRUE) - log(sum(stats::dnorm(z)))
  nll <- function(par) {
    th <- par$mu + exp(par$log_sd) * z
    log_s <- map(betas, \(b) pcm_kernel(th, b)$log_S)            # log S_i(theta)，一題一個向量
    ll <- 0
    for (j in seq_along(r)) {
      v <- r[j] * th + log_w
      for (i in items_of[[j]]) v <- v - log_s[[i]]
      lse <- v[1]
      for (q in seq_along(z)[-1]) lse <- RTMB::logspace_add(lse, v[q])   # log sum_q exp(v_q)
      ll <- ll + n[j] * lse
    }
    RTMB::ADREPORT(exp(par$log_sd))
    -ll
  }
  wle <- add_wle(persons, parms)$theta   # 每個人的 WLE（起始值用；kits_subj() 的程式）
  start <- list(mu = mean(wle), log_sd = log(max(stats::sd(wle), 0.1)))
  obj <- RTMB::MakeADFun(nll, start, silent = TRUE)
  fit <- fit_ad(obj)
  rep <- RTMB::sdreport(obj, par.fixed = fit$par, hessian.fixed = fit$hessian)
  par <- unname(fit$par)
  c(n = nrow(persons), mean = par[1], sd = exp(par[2]), se_mean = sqrt(rep$cov.fixed[1, 1]),
    se_sd = summary(rep, "report")[1, "Std. Error"],
    convergence = if (fit$max_gradient <= 1e-4 && isTRUE(fit$hessian_pd)) 0 else 1)   # 和 kits_item() 同一個標準
}

# 兩個世代都有的波次：每個世代在這一波的能力分佈（常態，邊際最大概似）；每個波次一條線（平均和 SD 對齊），
# 再加一條所有共同波次一起的線（每個波次的權重相同）
same_age_lines <- function(resp_from, resp_to, item_from, item_to, wave_from, wave_to) {
  shared <- sort(intersect(unique(wave_from), unique(wave_to)))
  if (length(shared) == 0) {
    return(NULL)
  }
  in_wave <- function(resp, waves, w) filter(resp, visit_id %in% names(waves)[waves == w])
  waves <- map(shared, \(w) {
    f <- latent_normal(in_wave(resp_from, wave_from, w), param_frame(item_from))
    t <- latent_normal(in_wave(resp_to, wave_to, w), param_frame(item_to))
    data.frame(wave_m = w, n_from = as.integer(f[["n"]]), mean_from = f[["mean"]], sd_from = f[["sd"]],
               se_mean_from = f[["se_mean"]], se_sd_from = f[["se_sd"]],
               n_to = as.integer(t[["n"]]), mean_to = t[["mean"]], sd_to = t[["sd"]],
               se_mean_to = t[["se_mean"]], se_sd_to = t[["se_sd"]],
               converged = f[["convergence"]] == 0 && t[["convergence"]] == 0)
  }) |>
    list_rbind()
  # 平均和 SD 配起來。SD 的比 sqrt(mean(sigma_to^2) / mean(sigma_from^2)) 當成 s_to、s_from = 1
  # （和 same_age_se() 的 slope 同一個數字）
  line_of <- function(rows) {
    ratio <- sqrt(mean(waves$sd_to[rows]^2) / mean(waves$sd_from[rows]^2))
    mean_sigma_moments(mean(waves$mean_from[rows]), 1, mean(waves$mean_to[rows]), ratio)
  }
  rows <- if (length(shared) > 1) c(list(seq_along(shared)), as.list(seq_along(shared))) else list(1)
  methods <- map(rows, \(r) {
    l <- line_of(r)
    data.frame(method = sprintf("same-age groups, wave%s %s", if (length(r) > 1) "s" else "", paste(waves$wave_m[r],
        collapse = " + ")),
               intercept = l[["intercept"]], slope = l[["slope"]],
               basis = sprintf("%s person-waves of from, %s of to (normal ability distributions matched)",
                               format(sum(waves$n_from[r])), format(sum(waves$n_to[r]))))
  }) |>
    list_rbind()
  list(methods = methods, waves = kits_table(waves), rows = rows)
}

# 同齡的線在 theta 的抽樣誤差（delta 法；各波次、兩個世代的平均和 SD 當成彼此獨立）。
# value = mean(mu_to) + slope (theta - mean(mu_from))，slope = sqrt(S_to / S_from)，S = 各波次 sigma^2 的平均。
# 只寫 value：它對各平均、SD 的偏微分由 RTMB 給（jacobian，一列一個 theta），SE = sqrt(sum J^2 se^2)
same_age_se <- function(waves, rows, theta) {
  w <- waves[rows, , drop = FALSE]
  W <- length(rows)
  x <- c(w$mean_to, w$mean_from, w$sd_to, w$sd_from)
  se <- c(w$se_mean_to, w$se_mean_from, w$se_sd_to, w$se_sd_from)
  line_at <- function(x) {
    part <- function(j) x[(j - 1) * W + seq_len(W)]
    ratio <- sqrt(sum(part(3)^2) / sum(part(4)^2))
    line_apply(mean_sigma_moments(sum(part(2)) / W, 1, sum(part(1)) / W, ratio), theta)
  }
  J <- matrix(RTMB::MakeTape(line_at, x)$jacobian(x), length(theta))
  sqrt(drop(J^2 %*% se^2))
}

# 合併校準，和每個單獨校準到合併量尺的線（所有門檻的 mean/sigma），推出 from -> to 的線
concurrent_lines <- function(resp_from, resp_to, item_from, item_to) {
  resp_from$visit_id <- paste0("from\u001f", resp_from$visit_id)
  resp_to$visit_id <- paste0("to\u001f", resp_to$visit_id)
  conc <- kits_item(resp_mark(bind_rows(resp_from, resp_to)), scale = "concurrent")
  ct <- select(as.data.frame(conc$thresholds), item_id, item_score, thr)
  one <- function(item) {
    m <- select(as.data.frame(item$thresholds), item_id, item_score, thr) |>
      inner_join(ct, by = c("item_id", "item_score"), suffix = c("_sep", "_conc")) |>
      arrange(sort_rank(item_id), item_score)
    c(algo_mean_sigma(m$thr_sep, m$thr_conc), n = nrow(m))
  }
  f <- one(item_from)
  t <- one(item_to)
  implied <- line_compose(f, algo_line(t, invert = TRUE))   # from -> 合併量尺 -> to
  lines <- kits_table(data.frame(calibration = c("from", "to"), scale = c(scale_of(item_from), scale_of(item_to)),
                                 intercept = c(f[["intercept"]], t[["intercept"]]),
                                 slope = c(f[["slope"]], t[["slope"]]),
                                 steps = as.integer(c(f[["n"]], t[["n"]]))))
  method <- data.frame(method = "concurrent CML",
                       intercept = implied$intercept, slope = implied$slope,
                       basis = sprintf("%d items, %d person-waves calibrated together", length(unique(ct$item_id)),
                                       conc$n_persons))
  list(item = conc, lines = lines, method = method)
}

link_check_diagnostics <- function(x, draws, has_persons) {
  at <- x$at
  rows <- list()
  for (m in x$methods$method[-1]) {
    d <- filter(as.data.frame(at), method == m)
    if (all(is.na(d$diff))) {
      rows[[length(rows) + 1]] <- diag_row(m, NA, "note", "could not be computed")
      next
    }
    rows[[length(rows) + 1]] <- if (draws) method_z_row(m, d) else method_tol_row(m, d, x$settings$tol)
  }
  if (is.null(x$waves)) {
    skipped <- if (has_persons) "no wave is shared by the two calibrations: skipped" else "no `persons` with wave_m: skipped"
    rows[[length(rows) + 1]] <- diag_row("same-age groups", 0, "ok", skipped)
  } else {
    bad <- x$waves$wave_m[!x$waves$converged]
    if (length(bad) > 0) {
      rows[[length(rows) + 1]] <- diag_row("same-age groups: convergence", length(bad), "warn",
        sprintf("the ability distribution at wave(s) %s did not converge", paste(bad, collapse = ", ")))
    }
    rows[[length(rows) + 1]] <- diag_row("same-age groups: assumption", nrow(x$waves), "ok",
      paste("assumes the children of both cohorts at the same age come from one population (equivalent groups,",
            "Kolen & Brennan pp. 431-432); a cohort effect shows up as a difference"))
  }
  diag_table(rows)
}


#' @rdname kits_link_check
#' @param x A `kits_link_check`.
#' @param digits Decimals to print.
#' @param ... Not used.
#' @export
print.kits_link_check <- function(x, digits = 2, ...) {
  need_digits(digits)
  cat(sprintf(paste0("<kits_link_check> %s -> %s: the link and independent lines,",
                     " new value = intercept + slope x old value\n"),
              scale_label(x$from_scale), scale_label(x$to_scale)))
  m <- x$methods
  print(format_table(m[c("method", "intercept", "slope")], digits), row.names = FALSE,
        right = FALSE)
  cat(paste0("  ", m$method, ": ", m$basis), sep = "\n")
  if (!is.null(x$waves)) {
    cat("Shared waves (each cohort on its own calibration):\n")
    print(format_table(x$waves, digits), row.names = FALSE, right = TRUE)
  }
  cat(sprintf("At theta (scale %s; value on %s; diff = method minus the link):\n", scale_label(x$from_scale),
              scale_label(x$to_scale)))
  print(format_table(x$at, digits), row.names = FALSE, right = TRUE)
  checks_line(x$diagnostics)
  invisible(x)
}

#' @rdname summary.kits
#' @export
summary.kits_link_check <- kits_summary
