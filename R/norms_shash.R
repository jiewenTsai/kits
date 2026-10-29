# shash（sinh-arcsinh）用四個數描述一個月齡的分布：mu 位置、sigma 尺度、epsilon 偏態、delta 尾巴（內部用）。
# 公式和 cNORM 的 qshash() 相同。

# z 分數 -> 分數
shash_quantile <- function(z, mu, sigma, epsilon, delta) {
  mu + sigma * sinh((asinh(z) - epsilon) / delta)
}

# 分數 -> z 分數（shash_quantile() 的反函數）；分數可以是 0 或負的
shash_z <- function(D, mu, sigma, epsilon, delta) {
  sinh(delta * asinh((D - mu) / sigma) + epsilon)
}


#' Fit a continuous norm (shash by maximum likelihood in RTMB)
#'
#' Fits a sinh-arcsinh (shash) distribution of the scores whose four parameters are polynomials in
#' transformed age (by default the square root), and returns the norm table: the parameters and the
#' percentiles at each age of `ages`. Scores can be 0 or negative.
#' Norms are given only within the age range of the data, rounded out to whole months.
#' Screen the waves first with [kits_norm_screen()]; the three steps are [kits_norm_steps].
#'
#' Details: `kits_algorithm("norms")` and tutorial_algorithms/04_norms.html; the
#' models tutorial (tutorial_models/kits_models.html), section 8.
#'
#' @param d The score of each person-wave (for example the D score), or a [kits_subj()] result
#'   (its `theta` and `age_m` are used; a message says when its scale is not `"D"`).
#' @param age_m Age in months of each person-wave (not with a kits_subj: taken from its column).
#' @param ages The ages at which to output the norms. With a [kits_subj()] it can be left out: then
#'   the rule of [kits_norm_ages()] sets them, with `pad = 0` unless `pad` is given.
#' @param degrees The degrees of the four polynomials for mu, sigma, epsilon, delta. The default is
#'   for about 12 waves; with 2 to 5 waves use for example `c(3, 2, 1, 1)`, by sex `c(4, 3, 2, 1)`.
#'   At least `10 * sum(degrees + 1)` person-waves must be inside the ages (`kits_error_data`).
#' @param age_x The transformation function for age.
#' @param pct The percentiles to output.
#' @param pad Use data up to this many months outside the range of `ages` (with `pad = 0`, data
#'   outside are not used).
#' @param delta_starts Start values for delta: one fit per value, the largest likelihood is kept.
#' @param delta_max When the kept solution's delta exceeds this number, a warning is shown.
#' @param by Optional groups, one norm each: with a [kits_subj()], the names of its group columns
#'   (for example `"sex"`, or `c("domain", "sex")`); with scores as numbers, a data.frame of the
#'   group columns (one row per score). No group may be NA.
#' @param ... Not used: an argument the function does not have stops with an error (class
#'   `kits_error_input`).
#' @return A `kits_norms` (a data.frame; `plot()` draws the percentile curves) with one row per age:
#'   * `age_m`, and `mu`, `sigma`, `epsilon`, `delta`: location, scale, skewness, tail weight.
#'   * `n`, `neg_loglik` (the -log-likelihood of the kept solution).
#'   * the percentiles `P3`, `P5`, ...
#'   * `attr(, "kits_args")`: the settings (`ages_in`, `degrees`, `pad`, `delta_starts`, `age_x`,
#'     `pct`, `delta_max`, as in kitsPreview).
#'   * `attr(, "scale")`, `attr(, "diagnostics")` ([kits_check()]), `attr(, "info")`
#'     ([kits_info()]).
#'   * With `by`: the group columns first; `attr(, "by")` names them, `attr(, "groups")` keeps
#'     each norm.
#'
#'   Pass it to [kits_norm_z()] for same-age z-scores.
#' @examples
#' set.seed(1)
#' age <- runif(600, 3, 36)
#' D <- 20 + 15 * log(age) + rnorm(600, 0, 3)
#' norm <- kits_norms(D, age, ages = seq(3, 36, by = 3), degrees = c(3, 2, 1, 1))
#' library(dplyr)
#' norm |>
#'   select(age_m, P10, P50, P90)
#'
#' # The same norm in three steps
#' norm_data <- kits_norm_data(D, age, ages = seq(3, 36, by = 3))
#' fit <- kits_norm_fit(norm_data, degrees = c(3, 2, 1, 1))
#' all.equal(kits_norm_table(fit)$P50, norm$P50)
#' # The children binned: no single child in the plot
#' p <- plot(norm, scores = data.frame(theta = D, age_m = age))
#' if (interactive()) p
#'
#' # One norm per sex: one table, the group column first
#' sex <- rep(c("boy", "girl"), 300)
#' by_sex <- kits_norms(D, age, ages = seq(3, 36, by = 3), degrees = c(3, 2, 1, 1),
#'                      by = data.frame(sex = sex))
#' by_sex |>
#'   select(sex, age_m, P50)
#' kits_norm_z(c(50, 60), age_m = c(12, 24), by_sex, by = data.frame(sex = c("boy", "girl")))
#' @export
kits_norms <- function(d, ...) UseMethod("kits_norms")

#' @rdname kits_norms
#' @export
kits_norms.default <- function(d, age_m, ages, degrees = c(5, 4, 3, 2), age_x = sqrt, pct = NORM_PCT,
                               pad = 1, delta_starts = c(0.5, 1, 1.5, 2), delta_max = 3,
                               by = NULL, ...) {
  no_extra_args(..., fun = "kits_norms")
  call <- call_text(sys.call(), name = "kits_norms")
  age_x_text <- age_x_code(substitute(age_x), age_x)   # 例如 "sqrt"（和 kitsPreview 的 kits_args 相同）
  d_expr <- substitute(d)
  if (missing(ages)) {
    kits_abort("Give the ages of the norms: kits_norms(d, age_m, ages).", class = "kits_error_input")
  }
  input <- data_input(d, d_expr)
  settings <- list(degrees = degrees, age_x = age_x, age_x_text = age_x_text, pct = pct, pad = pad,
                   delta_starts = delta_starts, delta_max = delta_max)
  groups <- by_frame(by, length(d))
  fit_norms_by(d, age_m, ages, settings, NA_character_, call, input, groups)
}

#' @rdname kits_norms
#' @export
kits_norms.kits_subj <- function(d, ages, degrees = c(5, 4, 3, 2), age_x = sqrt, pct = NORM_PCT,
                                 pad = 1, delta_starts = c(0.5, 1, 1.5, 2), delta_max = 3,
                                 by = NULL, ...) {
  subj_extra_args(..., fun = "kits_norms", form = "kits_norms(s, ages, ...)")
  call <- call_text(sys.call(), name = "kits_norms")
  age_x_text <- age_x_code(substitute(age_x), age_x)
  d_expr <- substitute(d)
  s <- d
  input <- data_input(s, d_expr)
  auto <- missing(ages)
  if (auto) {   # kits_norm_ages() 的規則：第一波的 1% 到最後一波的 99%；範圍外的資料不用（pad = 0）
    ages <- kits_norm_ages(s)
    if (missing(pad)) pad <- 0
  }
  norm_scale_message(s)
  settings <- list(degrees = degrees, age_x = age_x, age_x_text = age_x_text, pct = pct, pad = pad,
                   delta_starts = delta_starts, delta_max = delta_max)
  groups <- by_columns(by, s)
  age_m <- subj_column(s, "age_m")
  out <- fit_norms_by(s$theta, age_m, ages, settings, scale_of(s), call, input, groups)
  info <- attr(out, "info")
  if (auto) info$settings$ages_from <- "kits_norm_ages()"
  if (auto) attr(out, "fitted")$set_by <- "kits_norm_ages()"
  info$settings$unscreened <- unscreened_waves(s)
  attr(out, "info") <- info
  out
}

# 沒有先篩的波次：kits_norm_screen() 會拿掉的波次（20% 以上在最低分）；沒有時 NULL
unscreened_waves <- function(s) {
  if (!"wave_m" %in% names(s)) {
    return(NULL)
  }
  screen <- tryCatch(suppressMessages(kits_norm_screen(s)), kits_error = function(e) NULL)
  if (is.null(screen) || all(screen$waves$normed)) {
    return(NULL)
  }
  out <- filter(as.data.frame(screen$waves), !normed)
  if ("cohort" %in% names(out)) paste(out$cohort, out$wave_m) else format(out$wave_m)
}

# ---- 照組配常模（by）------------------------------------------------------------------------------------

# 數字的分數：by 是一個 data.frame（每個分數一列）；沒給時 NULL
by_frame <- function(by, n) {
  if (is.null(by)) {
    return(NULL)
  }
  abort_if(is.character(by), "With scores as numbers, give `by` as a data.frame of the group columns, ",
           "e.g. by = data.frame(sex = sex); names of columns work with a kits_subj().")
  abort_if(!is.list(by) || length(by) == 0 || is.null(names(by)) || any(names(by) == ""),
           "`by` must be a data.frame of named group columns, e.g. data.frame(sex = sex).")
  groups <- as.data.frame(by)
  abort_if(nrow(groups) != n, sprintf("`by` needs one row per score (%d); it has %d.", n, nrow(groups)))
  need_groups(groups)
}

# kits_subj：by 是它的欄名
by_columns <- function(by, s) {
  if (is.null(by)) {
    return(NULL)
  }
  abort_if(!is.character(by) || length(by) == 0 || anyNA(by) || anyDuplicated(by),
           "`by` must be the names of group columns of the scores, e.g. by = \"sex\".")
  miss <- setdiff(by, names(s))
  abort_if(length(miss) > 0, "The scores have no column ", paste0("`", miss, "`", collapse = ", "),
           " for `by`; give it in the `persons` of kits_subj(), or add it to the scores.")
  need_groups(as.data.frame(s)[by])
}

# 組的欄：沒有 NA，文字或數字（因子轉成文字）
need_groups <- function(groups) {
  groups <- map(groups, \(v) if (is.factor(v)) as.character(v) else v) |>
    as.data.frame()
  has_na <- names(groups)[map_lgl(groups, anyNA)]
  abort_if(length(has_na) > 0, "The group column(s) ", paste(has_na, collapse = ", "), " of `by` have NA; ",
           "give every person-wave a group, or leave those out first.", class = "kits_error_data")
  groups
}

# 組（不重複、照順序）；每個人次屬於第幾組
group_keys <- function(groups) {
  keys <- groups |>
    distinct() |>
    arrange(pick(everything()))
  list(keys = keys, id = vctrs::vec_match(groups, keys))
}

# 組的名字，例如 "score = lan, sex = boy"
group_label <- function(key) {
  values <- map_chr(key, \(v) format(v))
  paste(names(key), "=", values, collapse = ", ")
}

# 沒有 by 時一個常模；有 by 時每組一個（和只用這組的孩子的 kits_norms() 相同），疊成一張表，組的欄在前
fit_norms_by <- function(d, age_m, ages, settings, scale, call, input, groups) {
  out <- fit_norms_groups(d, age_m, ages, settings, scale, call, input, groups)
  attr(out, "fitted") <- list(from = min(out$age_m), to = max(out$age_m), rows = nrow(out))
  out
}

# 常模的月齡範圍和列數記在 attr(, "fitted")：filter() 之後印出來時，還能說出配的範圍
fit_norms_groups <- function(d, age_m, ages, settings, scale, call, input, groups) {
  if (is.null(groups)) {
    return(fit_norms(d, age_m, ages, settings, scale = scale, call = call, input = input))
  }
  abort_if(length(d) != length(age_m), "`d` and `age_m` must have the same length.")
  g <- group_keys(groups)
  fit_one <- function(i) {
    rows <- g$id == i
    in_group(group_label(g$keys[i, , drop = FALSE]),
             fit_norms(d[rows], age_m[rows], ages, settings, scale = scale, call = call,
                       input = input))
  }
  norms <- vector("list", nrow(g$keys))
  for (i in seq_along(norms)) norms[[i]] <- fit_one(i)   # 不用 map()：錯誤照原樣（和它的類別）傳出去
  grouped_norms(norms, g$keys, call, input, values_fingerprint(d, age_m, groups))
}

# 一組的計算：錯誤和警告前面加上組的名字（錯誤的類別不變）
in_group <- function(label, expr) {
  tryCatch(
    withCallingHandlers(expr, kits_warning = function(w) {
      kits_warn(label, ": ", conditionMessage(w), class = class(w)[1])
      invokeRestart("muffleWarning")
    }),
    kits_error = function(e) kits_abort(label, ": ", conditionMessage(e), class = class(e)[1]))
}

# 每組的常模疊成一張表：組的欄在前；設定、量尺照第一組（每組相同），診斷每一列前面加上組的名字
grouped_norms <- function(norms, keys, call, input, fingerprint) {
  rows <- map(seq_along(norms), \(i) {
    key <- keys[rep(i, nrow(norms[[i]])), , drop = FALSE]
    rownames(key) <- NULL
    bind_cols(key, as.data.frame(norms[[i]]))
  })
  diagnostics <- map(seq_along(norms), \(i) {
    label <- group_label(keys[i, , drop = FALSE])
    mutate(as.data.frame(attr(norms[[i]], "diagnostics")), check = paste0(label, ": ", check))
  })
  out <- kits_table(list_rbind(rows), "kits_norms")
  first <- norms[[1]]
  settings <- c(attr(first, "info")$settings, list(by = names(keys)))
  attr(out, "kits_args") <- c(attr(first, "kits_args"), list(by = names(keys)))
  attr(out, "scale") <- attr(first, "scale")
  attr(out, "diagnostics") <- kits_table(list_rbind(diagnostics), "kits_diagnostics")
  attr(out, "info") <- make_info(call, settings, fingerprint, input = input)
  attr(out, "by") <- names(keys)
  attr(out, "groups") <- norms
  out
}

# 常模的三步（kits_norms() 的兩個方法共用；kits_norm_data()、kits_norm_fit()、kits_norm_table() 是同樣的三步）
fit_norms <- function(d, age_m, ages, settings, scale, call, input) {
  check_norm_settings(d, age_m, settings$degrees, settings$delta_starts)
  norm_data <- norm_data_step(d, age_m, ages, settings$pad, settings$age_x, settings$age_x_text, scale)   # 1 資料
  fit <- norm_fit_step(norm_data, settings$degrees, settings$delta_starts)                                 # 2 配適
  norm_table_step(fit, settings$pct, settings$delta_max, call, input)                                      # 3 常模表
}

check_norm_settings <- function(d, age_m, degrees, delta_starts) {
  need_d_age(d, age_m)
  check_norm_fit_settings(degrees, delta_starts)
}

check_norm_fit_settings <- function(degrees, delta_starts) {
  abort_if(length(degrees) != 4 || any(is.na(degrees)) || any(degrees < 1 | degrees != round(degrees)),
           "`degrees` needs 4 whole numbers of at least 1: the degrees for mu, sigma, epsilon, delta.")
  abort_if(length(delta_starts) == 0 || any(!is.finite(delta_starts) | delta_starts <= 0),
           "`delta_starts` needs at least one positive number.")
}

# 1 配常模用的資料：分數和月齡都有限、在 ages 的範圍（加減 pad）內、age_x 算得出來。
# 常模表的月齡只留資料的範圍內（往外取整到整月）。結果是 kits_norm_data（d、age_m 兩欄），設定在 attr(, "norm")
norm_data_step <- function(d, age_m, ages, pad, age_x, age_x_text, scale) {
  ages <- sort(unique(ages[is.finite(ages)]))
  abort_if(length(ages) < 2, "`ages` needs at least 2 ages.")
  ages_in <- range(ages)   # 篩資料用的範圍
  data <- data.frame(d = d, age_m = age_m) |>
    filter(is.finite(d), is.finite(age_m)) |>
    filter(age_m >= ages_in[1] - pad, age_m <= ages_in[2] + pad) |>
    filter(is.finite(age_x(age_m)))
  abort_if(nrow(data) < 2, nrow(data), " person-wave(s) within the ages; a norm needs data.", class = "kits_error_data")
  abort_if(stats::sd(data$d) == 0,
           "All scores within the ages are the same (", data$d[1], "), so no norm can be fitted.", class = "kits_error_data")
  in_data <- ages >= floor(min(data$age_m)) & ages <= ceiling(max(data$age_m))
  ages <- ages[in_data & is.finite(age_x(ages))]
  abort_if(length(ages) < 2 || length(unique(data$age_m)) < 2,
           "The data cover fewer than 2 of the `ages`; a norm over age cannot be fitted.", class = "kits_error_data")
  out <- kits_table(data, "kits_norm_data")
  attr(out, "norm") <- list(ages = ages, ages_in = ages_in, pad = pad, age_x = age_x, age_x_text = age_x_text, scale = scale)
  out
}

# 2 shash 的最大概似：每個 delta 的起始值配一次，留下概似最大的解；資料太少時停下來
norm_fit_step <- function(norm_data, degrees, delta_starts) {
  setting <- attr(norm_data, "norm")
  abort_if(nrow(norm_data) < 10 * sum(degrees + 1),
           nrow(norm_data), " person-waves within the ages; too few for degrees ", paste(degrees, collapse = ", "),
           " (at least ", 10 * sum(degrees + 1), "). Use lower degrees or more data.", class = "kits_error_data")
  model <- best_shash_fit(setting$age_x(norm_data$age_m), norm_data$d, degrees, delta_starts)
  structure(list(model = model, data = norm_data, degrees = degrees, delta_starts = delta_starts), class = "kits_norm_fit")
}

# 3 常模表：每個月齡的 mu、sigma、epsilon、delta（和配適時一樣夾住）和百分位，加上設定、量尺、診斷、出處
norm_table_step <- function(fit, pct, delta_max, call, input = NULL) {
  setting <- attr(fit$data, "norm")
  ages <- setting$ages
  par <- shash_parameters(fit$model, setting$age_x(ages))
  if (max(par$delta) > delta_max) {
    kits_warn("delta reaches ", round(max(par$delta), 2), " even at the best start (light tails, or degrees too high). ",
              lower_degrees_advice(fit$degrees), class = "kits_warning_norms")
  }
  out <- data.frame(age_m = ages, mu = par$mu, sigma = par$sigma, epsilon = par$epsilon, delta = par$delta,
                    n = nrow(fit$data), neg_loglik = fit$model$value)
  out <- bind_cols(out, norm_percentiles(out, pct))
  attr(out, "kits_args") <- list(ages_in = setting$ages_in, degrees = fit$degrees, pad = setting$pad,
                                 delta_starts = fit$delta_starts, age_x = setting$age_x_text, pct = pct,
                                 delta_max = delta_max)
  out <- kits_table(out, "kits_norms")
  attr(out, "scale") <- setting$scale
  attr(out, "diagnostics") <- norm_diagnostics(fit$model, par, delta_max, fit$data, setting$age_x)
  settings <- list(degrees = fit$degrees, age_x = setting$age_x_text, pad = setting$pad, delta_starts = fit$delta_starts,
                   delta_max = delta_max, ages = range(ages))
  attr(out, "info") <- make_info(call, settings, values_fingerprint(fit$data$d, fit$data$age_m), input = input)
  attr(out, "fitted") <- list(from = min(ages), to = max(ages), rows = nrow(out))
  out
}

# 每個 delta 的起始值配一次，留下概似最大的解；警告先收起來，只發出留下的那個解的
best_shash_fit <- function(x, y, degrees, delta_starts) {
  fits <- map(delta_starts, \(d0) keep_warnings(algo_shash(x, y, degrees, d0)))
  best <- fits[[which.min(map_dbl(fits, \(f) f$result$value))]]
  messages <- map_chr(best$warnings, conditionMessage)
  walk(best$warnings[!duplicated(messages)], warning)
  best$result
}

# 每個百分位一欄（P3、P5、...）：P_q = mu + sigma sinh((asinh(qnorm(q / 100)) - epsilon) / delta)
norm_percentiles <- function(par, pct) {
  columns <- map(pct, \(q) shash_quantile(stats::qnorm(q / 100), par$mu, par$sigma, par$epsilon, par$delta))
  names(columns) <- paste0("P", pct)
  as.data.frame(columns, check.names = FALSE)
}

# expr 的結果和它發出的警告；警告先不發出，留給呼叫的人決定
keep_warnings <- function(expr) {
  kept <- new.env()
  kept$warnings <- list()
  result <- withCallingHandlers(expr, warning = function(w) {
    kept$warnings <- c(kept$warnings, list(w))
    invokeRestart("muffleWarning")
  })
  list(result = result, warnings = kept$warnings)
}

# 收斂的說明（kits_norm_fit 的 print 和診斷）：nlminb 的代碼，和 kits_item() 相同的標準
norm_convergence_text <- function(model) {
  nlminb <- sprintf("nlminb code %d: %s", model$convergence, model$message)
  if (model$converged) {
    return(sprintf("converged (%s; max |gradient| < 1e-4, Hessian of the free coefficients positive definite)", nlminb))
  }
  sprintf("not converged (%s; max |gradient| %.1e, tolerance 1e-4; Hessian of the free coefficients %s)", nlminb,
          model$grad_max, if (model$hessian_pd) "positive definite" else "not positive definite")
}

# kits_norms() 的診斷：收斂、係數在界上、最大的 delta、配適資料的同齡 z（每個年齡十分位）
# data：配適用的資料（kits_norm_data，d、age_m 兩欄）
norm_diagnostics <- function(model, par, delta_max, data, age_x) {
  conv_status <- if (model$converged) "ok" else "warn"
  grad_status <- if (model$grad_max > 1e-4) "warn" else "ok"
  grad_message <- "largest |d(-log L) / d coefficient| over the free coefficients (tolerance 1e-4)"
  at_bound <- length(model$at_bound) > 0
  bound_status <- if (at_bound) "warn" else "ok"
  bound_message <- if (at_bound) paste0("coefficient(s) ", paste(model$at_bound, collapse = ", "),
                                        "; try lower degrees") else "none"
  heavy <- max(par$delta) > delta_max
  delta_status <- if (heavy) "warn" else "ok"
  delta_detail <- if (heavy) sprintf(" (above delta_max = %g)", delta_max) else ""
  rows <- list(
    diag_row("optimizer convergence", model$convergence, conv_status, norm_convergence_text(model)),
    diag_row("max |gradient|", model$grad_max, grad_status, grad_message),
    diag_row("coefficients at their limits", length(model$at_bound), bound_status, bound_message),
    diag_row("max delta", max(par$delta), delta_status, sprintf("largest tail parameter over the ages%s", delta_detail)))
  # 同齡 z：用模型在每個人的月齡的參數（不經過常模表的內插）
  z <- algo_shash_z(model, age_x(data$age_m), data$d)
  by_decile <- z_by_age_decile(z, data$age_m)
  worst_m <- slice(by_decile, which.max(abs(mean)))
  sd_off <- slice(by_decile, which.max(abs(sd - 1)))
  bad_m <- filter(by_decile, abs(mean) > 0.2)$decile
  bad_sd <- filter(by_decile, !is.na(sd), sd < 0.8 | sd > 1.2)$decile
  mean_off <- length(bad_m) > 0
  mean_status <- if (mean_off) "warn" else "ok"
  mean_message <- if (mean_off) {
    paste0("mean z beyond 0.2 at ages ", paste(bad_m, collapse = ", "), ": the curves miss the data there")
  } else {
    sprintf("largest at ages %s", worst_m$decile)
  }
  rows$zm <- diag_row("same-age z: max |mean| by age decile", abs(worst_m$mean), mean_status, mean_message)
  spread_off <- length(bad_sd) > 0
  sd_status <- if (spread_off) "warn" else "ok"
  sd_message <- if (spread_off) {
    paste0("SD of z outside 0.8 to 1.2 at ages ", paste(bad_sd, collapse = ", "))
  } else {
    sprintf("at ages %s", sd_off$decile)
  }
  rows$zs <- diag_row("same-age z: SD furthest from 1 by age decile", sd_off$sd, sd_status, sd_message)
  diag_table(rows)
}

# 同齡 z 的平均和 SD，每個年齡十分位一列（照年齡排；十分位的界都一樣時只有一組 "all"）
z_by_age_decile <- function(z, age_m) {
  breaks <- unique(stats::quantile(age_m, seq(0, 1, by = 0.1), names = FALSE))
  decile <- if (length(breaks) > 2) cut(age_m, breaks, include.lowest = TRUE) else factor(rep("all", length(age_m)))
  data.frame(z = z, decile = decile) |>
    summarise(mean = mean(z), sd = stats::sd(z), .by = decile) |>
    arrange(decile) |>
    mutate(decile = as.character(decile))
}


#' @rdname kits_norms
#' @param x A `kits_norms`.
#' @param digits Decimals to print (only the printout is rounded; the object keeps full precision).
#' @param n Rows to print: by default all rows of a table of up to 30 rows, else the first 10
#'   (see [print.kits_table()]).
#' @param ... Passed to `print()` for a data.frame.
#' @export
print.kits_norms <- function(x, digits = 2, n = NULL, ...) {
  need_digits(digits)
  need_rows(n)
  cat(sprintf("<kits_norms> scale: %s\n", scale_label(scale_of(x))))
  fitted <- attr(x, "fitted")
  if (!is.null(fitted)) {
    from <- fitted$set_by
    by <- attr(x, "by")
    shown <- if (nrow(x) < fitted$rows) sprintf("; showing %s", plural(nrow(x), "row"))
    cat(sprintf("  norm fitted on %s to %s months%s%s%s\n", format(fitted$from), format(fitted$to),
                if (is.null(from)) "" else sprintf(" (set by %s)", from),
                if (is.null(by)) "" else sprintf("; one norm per %s", paste(by, collapse = " x ")),
                if (is.null(shown)) "" else shown))
  }
  print_rows(x, digits, n, ...)
  unscreened <- attr(x, "info")$settings$unscreened
  if (length(unscreened) > 0) {
    cat(sprintf("Note: wave(s) %s have 20%% or more at the lowest score; ",
                paste(unscreened, collapse = ", ")),
        "kits_norm_screen() leaves them out. ",
        "Screen first: s |> filter(kits_norm_screen(s)$normed) |> kits_norms().\n", sep = "")
  }
  cat("Next: kits_norm_z(scores, x), plot(x), plot(x, scores = scores).\n")
  checks_line(attr(x, "diagnostics"))
  invisible(x)
}


#' Plot a continuous norm: percentile curves, one norm or by group
#'
#' One norm: the percentile curves, with the children's scores optionally drawn behind them
#' (binned by default, so the plot and its html file hold no single child).
#' A norm by group (`kits_norms(..., by = )`): one panel per domain x group, each group in its own
#' colour.
#'
#' @seealso `plot()` of a [kits_item_table()]: the growth curves with the items.
#' @param x The result of [kits_norms()] (it can have a `segment` column).
#' @param scores Optional: the scores of the children (a [kits_subj()] result, or a data.frame with
#'   `theta` and `age_m`; for a norm by group also the group columns), drawn behind the curves.
#' @param pct Which percentiles to draw; they must be columns of `x`.
#' @param persons How to draw the children: `"bins"` (default, counts per cell of a 40 x 40 grid)
#'   or `"points"` (one point per child; for your own checks only, not for plots you share).
#' @param min_cell With `persons = "bins"`, cells with fewer children than this are not drawn
#'   (default 5).
#' @param age_ticks By group: the age ticks.
#' @param y_range By group: the range of the y axis, the same in every panel; by default one common
#'   range from the data of all panels (rounded to multiples of 10).
#' @param ... Not used.
#' @return A plotly plot.
#' @examples
#' norm <- data.frame(age_m = c(6, 12, 18), mu = c(40, 55, 62), sigma = 5, epsilon = 0, delta = 1,
#'                    P3 = c(31, 46, 53), P50 = c(40, 55, 62), P97 = c(49, 64, 71))
#' class(norm) <- c("kits_norms", "kits_table", "data.frame")
#' p <- plot(norm, pct = c(3, 50, 97))
#' if (interactive()) p
#'
#' # One norm per sex: one panel per sex
#' set.seed(1)
#' age <- runif(1200, 3, 24)
#' sex <- rep(c("boy", "girl"), 600)
#' D <- 20 + 15 * log(age) + ifelse(sex == "girl", 1, 0) + rnorm(1200, 0, 3)
#' by_sex <- kits_norms(D, age, ages = 3:24, degrees = c(3, 2, 1, 1), by = data.frame(sex = sex))
#' p <- plot(by_sex)
#' if (interactive()) p
#' @export
plot.kits_norms <- function(x, scores = NULL, pct = c(3, 10, 50, 90, 97),
                            persons = c("bins", "points"),
                            min_cell = 5, age_ticks = AGE_TICKS, y_range = NULL, ...) {
  no_extra_args(..., fun = "plot.kits_norms")
  persons <- match_option(persons, c("bins", "points"), "persons")
  need_number(min_cell, "min_cell", "one number of at least 1", lower = 1)
  abort_if(length(pct) == 0, "`pct` needs at least one percentile.")
  children <- norm_children(scores)
  if (!is.null(attr(x, "by"))) {
    return(norm_groups_plot(x, pct, age_ticks, y_range, children, persons, min_cell))
  }
  cols <- paste0("P", pct)
  need_cols(as.data.frame(x), c("age_m", cols), "x")
  p <- plotly::plot_ly()
  if (!is.null(children)) p <- add_children(p, children, persons, min_cell)
  segments <- norm_segments(as.data.frame(x))
  for (s in names(segments)) {
    d <- segments[[s]]
    for (q in cols) {
      p <- plotly::add_lines(p, x = d$age_m, y = d[[q]], name = q, legendgroup = q, showlegend = s == names(segments)[1],
                             line = list(color = PAL$blue, width = if (q == "P50") LINE_W$bold else LINE_W$thin,
                                         dash = if (q %in% c("P3", "P97")) "dot" else "solid"),
                             hovertext = sprintf("%s, %.1f%s: %.1f", q, d$age_m, LAB$month, d[[q]]), hoverinfo = "text")
    }
  }
  plotly::layout(p, xaxis = list(title = LAB$age), yaxis = list(title = LAB$D)) |>
    kits_style(LAB$title_norm_curves)
}

# 畫在曲線後面的孩子：分數表（theta、age_m，和組的欄）
norm_children <- function(scores) {
  if (is.null(scores)) {
    return(NULL)
  }
  abort_if(!is.data.frame(scores),
           "`scores` must be the scores table of the children (theta and age_m), ", "e.g. from kits_subj().")
  need_cols(scores, c("theta", "age_m"), "scores")
  as.data.frame(scores)
}

# 孩子：分組畫（每格的人數，少於 min_cell 的格不畫），或一個孩子一個點
add_children <- function(p, children, persons, min_cell) {
  if (persons == "points") {
    marker <- list(color = PAL$light, size = 3)
    return(plotly::add_markers(p, x = children$age_m, y = children$theta, marker = marker,
                               name = LAB$persons,
                               hoverinfo = "skip"))
  }
  add_person_bins(p, children$theta, children$age_m, min_cell)
}

# 照組的常模：一列一個領域（domain 是組的欄時），一欄一個其他的組（例如性別）；孩子畫在自己那一組的格子
norm_groups_plot <- function(x, pct, age_ticks, y_range, children, persons, min_cell) {
  by <- attr(x, "by")
  others <- setdiff(by, "domain")
  d <- as.data.frame(x)
  need_cols(d, c("age_m", paste0("P", pct)), "x")
  panel_keys <- function(t) {
    t$domain <- if ("domain" %in% by) as.character(t$domain) else NA_character_
    t$sex <- if (length(others) == 0) NA_character_ else group_text(t[others])
    t
  }
  d <- panel_keys(d)
  if (!is.null(children)) {
    need_cols(children, by, "scores")
    children <- panel_keys(children)
  }
  title <- paste0("Continuous norms",
                  if (length(others) > 0) paste0(" by ", paste(others, collapse = " and ")),
                  ": ", paste0("P", pct, collapse = ", "))
  y_title <- if (is.na(scale_of(x))) LAB$D else scale_of(x)
  overlay <- list(children = children, persons = persons, min_cell = min_cell)
  norm_grid_plot(d, NULL, unique(d$sex), pct, age_ticks, y_range, title, y_title, overlay)
}

# 幾個組的欄接成一個名字（例如 "M03, boy"）；一欄時就是它的值
group_text <- function(cols) {
  cols |>
    map(as.character) |>
    purrr::reduce(\(a, b) paste(a, b, sep = ", "))
}

# 孩子的分數分組畫（隱私）：約 40 x 40 格，每格的人數；少於 min_cell 人的格不畫，圖裡沒有單一個孩子
add_person_bins <- function(p, d, age_m, min_cell, n_bins = 40) {
  ok <- is.finite(d) & is.finite(age_m)
  if (!any(ok)) {
    return(p)
  }
  edges <- function(v) {
    r <- range(v)
    if (r[1] == r[2]) r <- r + c(-0.5, 0.5)
    seq(r[1], r[2], length.out = n_bins + 1)
  }
  ex <- edges(age_m[ok])
  ey <- edges(d[ok])
  ix <- cut(age_m[ok], ex, include.lowest = TRUE, labels = FALSE)
  iy <- cut(d[ok], ey, include.lowest = TRUE, labels = FALSE)
  counts <- matrix(tabulate(ix + (iy - 1) * n_bins, nbins = n_bins^2), n_bins, n_bins)   # 列是年齡、欄是分數
  counts[counts < min_cell] <- NA
  if (all(is.na(counts))) {
    return(p)
  }
  mids <- function(e) (utils::head(e, -1) + e[-1]) / 2
  plotly::add_heatmap(p, x = mids(ex), y = mids(ey), z = t(counts), colorscale = list(c(0, PAL$light), c(1, PAL$gray)),
                      showscale = FALSE, name = LAB$persons, hovertemplate = "n = %{z}<extra></extra>")
}

# shash 的界（和 cNORM 相同）：log sigma 和 epsilon 在 -10..10，log delta 在 -2..2。
# 配適（algo_shash() 的概似和係數的界）和常模表（shash_parameters()）都用這一個常數。
SHASH_BOUNDS <- c(log_sigma = 10, epsilon = 10, log_delta = 2)

# pmin(pmax(v, -bound), bound)，寫成 abs() 的形式：RTMB 不能對比大小（pmin、pmax）自動微分，abs() 可以
ad_clamp <- function(v, bound) (abs(v + bound) - abs(v - bound)) / 2

# 四個參數在每一列的值：X 是四個設計矩陣（多項式），coef 是四組係數（mu、sigma、epsilon、delta），
# clamp_fn(v, bound) 把 v 夾在 -bound..bound（配適時用可微分的 ad_clamp()，算常模表時用 pmin/pmax）
shash_link <- function(X, coef, clamp_fn) {
  list(mu = drop(X$mu %*% coef$mu),
       sigma = exp(clamp_fn(drop(X$sigma %*% coef$sigma), SHASH_BOUNDS[["log_sigma"]])),
       epsilon = clamp_fn(drop(X$epsilon %*% coef$epsilon), SHASH_BOUNDS[["epsilon"]]),
       delta = exp(clamp_fn(drop(X$delta %*% coef$delta), SHASH_BOUNDS[["log_delta"]])))
}

# shash 的最大概似（RTMB；x 是已轉換的月齡）：algo_shash()（algo.R，給 kits_norms()）

# 次數太高時的建議：說出參數的名字和較低的值
lower_degrees_advice <- function(degrees) {
  paste0("Lower the polynomial degrees with the `degrees` argument (now c(", paste(degrees, collapse = ", "), ")), ",
         paste0("for example degrees = c(3, 2, 1, 1): the default c(5, 4, 3, 2) suits many waves (KIT has about 12 over 3 to 72 ",
                "months), "),
         "not 2 to 5 waves. Then check the percentile curves with plot().")
}

# nlminb 停在梯度還不夠小的地方時（大 N 的 -log L 很平）：在不在界上的係數走牛頓步（RTMB 的 Hessian），
# 每一步都要讓 -log L 不變大（容許浮點誤差，和 kits_item() 的 fit_ad() 相同）、留在界內，否則停；在界上的係數不動
shash_newton <- function(result, obj, lower, upper, max_newton = 5) {
  par <- result$par
  value <- obj$fn(par)
  for (step in seq_len(max_newton)) {
    free <- !(is.finite(lower) & (par <= lower + 1e-6 | par >= upper - 1e-6))
    H <- obj$he(par)[free, free, drop = FALSE]
    if (any(!is.finite(H)) || !is_pd(H)) break
    delta <- solve(H, obj$gr(par)[free])
    par_new <- par
    par_new[free] <- pmin(pmax(par[free] - delta, lower[free]), upper[free])
    v_new <- obj$fn(par_new)
    if (!is.finite(v_new) || v_new > value + 1e-9 + 1e-12 * abs(value)) break
    par <- par_new
    value <- v_new
    if (max(abs(delta)) < 1e-10) break
  }
  utils::modifyList(result, list(par = par, objective = value))
}

# 收斂的判斷和 kits_item() 相同：不在界上的（自由的）係數，最大 |梯度| <= 1e-4，而且它們的 Hessian 正定
shash_convergence <- function(par, gr, he, lower, upper) {
  at_bound <- is.finite(lower) & (par <= lower + 1e-6 | par >= upper - 1e-6)
  grad_max <- max(abs(gr(par)[!at_bound]), 0)
  hessian_pd <- is_pd(he(par)[!at_bound, !at_bound, drop = FALSE])
  list(converged = grad_max <= 1e-4 && hessian_pd, grad_max = grad_max, hessian_pd = hessian_pd,
       at_bound = which(at_bound))
}

# algo_shash() 的參數在指定的（已轉換）月齡的值。
# 和配適時一樣夾住：log sigma、epsilon 在 -10..10，log delta 在 -2..2（cNORM 的 predictCoefficients_shash() 不夾，
# 所以多項式超出界的月齡，cNORM 的常模表和它配出來的模型不同；kits 的常模表就是配出來的模型）。
shash_parameters <- function(model, x) {
  xs <- (x - model$age_mean) / model$age_sd
  coef <- list(mu = model$mu_est, sigma = model$sigma_est, epsilon = model$epsilon_est,
               delta = model$delta_est)
  X <- map(coef, \(cf) outer(xs, seq_along(cf) - 1, "^"))
  clamp <- function(v, bound) pmin(pmax(v, -bound), bound)
  as.data.frame(shash_link(X, coef, clamp))
}

# age_x 寫成文字：base 的函數用名字（sqrt、log），其他的寫出函數本身（和 kitsPreview 的 kits_args 相同）
age_x_code <- function(expr, fun) {
  if (is.symbol(expr) && exists(as.character(expr), envir = baseenv()) && identical(fun, get(as.character(expr),
      envir = baseenv()))) {
    return(as.character(expr))
  }
  paste(deparse(fun, control = c("keepInteger", "niceNames")), collapse = "\n")
}
