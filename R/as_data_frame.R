# as.data.frame() 給串列的結果：一列（或每個單位一列）的摘要，好用 dplyr::bind_rows() 疊起來比較。
# 數字直接從物件拿（和 print()、summary() 相同）；checks 是診斷裡最嚴重的狀態（warn、note、ok）。

#' One-row summaries of kits results, ready to stack
#'
#' `as.data.frame()` turns a result that is a list into a plain data.frame with one row
#' (or one row per wave, per method, per item pattern), so that several results can be stacked with
#' `dplyr::bind_rows()` and compared, for example the links of every domain.
#' The numbers are those of the object (the same as `print()` and `summary()`).
#' `checks` is the most serious status of the diagnostics (`"warn"`, `"note"` or `"ok"`; see [kits_check()]).
#'
#' The columns:
#' * A link ([kits_linking()], [kits_link_line()]), one row:
#'   `from`, `to` (the scale labels; several `from` scales are joined with ", "), `kind`
#'   (`"calibration"`, `"anchors"`, or NA for your own points), `intercept`, `slope`,
#'   `intercept_error`, `slope_error`, `cov_intercept_slope` (their covariance over the draws),
#'   `intercept0`, `slope0` (the line before screening),
#'   `n_points`, `n_dropped`, `n_items`, `n_failed` (resamples that could not be computed), `cut`,
#'   `spread`, `unit`, `step` (`"all"` or the step, as text), `n_boot`, `seed`, `key` (the dscore
#'   key of the anchors), `checks`. Values that do not apply are NA.
#' * [kits_drift()], one row: `scale` (of `ref`), `intercept`, `slope`, `n_items`, `n_judged`,
#'   `n_drift` (judged items with a dropped step), `n_fix` (items suggested for fixing), `n_steps`,
#'   `n_steps_judged`, `n_steps_dropped`, `cut`, `min_n`, `checks`.
#' * [kits_link_check()], one row per method: `from`, `to`, `method`, `intercept`, `slope`, `basis`,
#'   `max_abs_diff`
#'   (the largest |diff| over `theta`), `max_abs_z`, `max_abs_z_both` (NA without bootstrap draws), `checks`.
#' * [kits_info()] (the provenance of any result), one row: `call`, `version`, `settings` (as text,
#'   `"name = value; ..."`), `n_rows`, `data_md5` (of the responses, or of the values), `persons_md5`, `input`
#'   (the name of the data object, when it was a data.frame or a kits_subj).
#' * [kits_item()], one row: `scale`, `n_persons`, `n_items`, `n_steps`, `n_patterns`, `n_fixed`
#'   (fixed steps), `n_left_out` (fixed items left out), `neg_loglik`, `convergence`,
#'   `max_gradient`, `checks`.
#' * [kits()], one row: the columns of [kits_item()], then `theta_mean`, `theta_sd`, `n_zero`, `n_perfect`,
#'   `n_misfit` (persons by `flag`).
#' * [kits_norm_screen()], one row per (cohort x) wave: its `waves` table.
#' * [kits_norm_fit()], one row: `n`, `degrees` (as text), `neg_loglik`, `convergence`, `converged`, `grad_max`,
#'   `n_at_bound`.
#' * [kits_stats()], one row per item pattern: `pattern_id`, `n_items`, `items` (as text), `n`,
#'   `max_raw`.
#'
#' @param x A kits result (see above).
#' @param row.names,optional Not used (as in [base::as.data.frame()]).
#' @param ... Not used.
#' @return A data.frame.
#' @examples
#' betas <- list(a1 = -2, a2 = -1, a3 = c(-0.5, 0.5), a4 = 0, a5 = c(0.5, 1.5), a6 = 2, a7 = 3)
#' young <- kits_sim(rnorm(500, 0), betas[1:5], seed = 1)
#' old <- kits_sim(rnorm(500, 2), betas[3:7], prefix = "q", seed = 2)
#' item_young <- kits_item(young, scale = "young cohort")
#' item_old <- kits_item(old, scale = "old cohort")
#' link <- kits_linking(item_old, item_young, n_boot = 100)
#' as.data.frame(link)
#'
#' # Several links in one table
#' link_step1 <- kits_linking(item_old, item_young, step = 1, n_boot = 100)
#' dplyr::bind_rows(as.data.frame(link), as.data.frame(link_step1))
#'
#' # The calibrations and their provenance
#' dplyr::bind_rows(as.data.frame(item_young), as.data.frame(item_old))
#' as.data.frame(kits_info(item_young))
#' @name as.data.frame.kits
NULL

# 診斷裡最嚴重的狀態；沒有診斷時 NA
worst_status <- function(diagnostics) {
  if (is.null(diagnostics) || nrow(diagnostics) == 0) {
    return(NA_character_)
  }
  levels <- c("warn", "note", "ok")
  levels[min(match(diagnostics$status, levels))]
}

# 一個值；NULL 或空的時候 NA（例如自己準備的點的 kind）
value_or_na <- function(x, na = NA) {
  if (is.null(x) || length(x) == 0) na else x
}

#' @rdname as.data.frame.kits
#' @export
as.data.frame.kits_link <- function(x, row.names = NULL, optional = FALSE, ...) {
  st <- x$settings
  data.frame(from = paste(x$from_scale, collapse = ", "),
             to = as.character(x$to_scale),
             kind = value_or_na(st$kind, NA_character_),
             intercept = x$intercept,
             slope = x$slope,
             intercept_error = x$intercept_error,
             slope_error = x$slope_error,
             cov_intercept_slope = if (is.matrix(x$vcov)) x$vcov[1, 2] else NA_real_,
             intercept0 = x$intercept0,
             slope0 = x$slope0,
             n_points = x$n,
             n_dropped = length(x$dropped),
             n_items = value_or_na(x$n_items, NA_integer_),
             n_failed = x$n_failed,
             cut = st$cut,
             spread = value_or_na(st$spread, "sd"),
             unit = value_or_na(st$unit, "step"),
             step = as.character(value_or_na(st$step, NA_character_)),
             n_boot = st$n_boot,
             seed = as.numeric(value_or_na(st$seed, NA_real_)),
             key = value_or_na(st$key, NA_character_),
             checks = worst_status(x$diagnostics))
}

#' @rdname as.data.frame.kits
#' @export
as.data.frame.kits_drift <- function(x, row.names = NULL, optional = FALSE, ...) {
  items <- as.data.frame(x$items)
  points <- as.data.frame(x$points)
  settings <- x$info$settings
  data.frame(scale = as.character(value_or_na(x$scale, NA_character_)),
             intercept = x$intercept,
             slope = x$slope,
             n_items = nrow(items),
             n_judged = sum(!is.na(items$drift)),
             n_drift = sum(items$drift, na.rm = TRUE),
             n_fix = sum(items$suggest_fix),
             n_steps = nrow(points),
             n_steps_judged = sum(points$judged),
             n_steps_dropped = sum(points$dropped, na.rm = TRUE),
             cut = settings$cut,
             min_n = settings$min_n,
             checks = worst_status(x$diagnostics))
}

#' @rdname as.data.frame.kits
#' @export
as.data.frame.kits_link_check <- function(x, row.names = NULL, optional = FALSE, ...) {
  at <- as.data.frame(x$at)
  has_z <- "z" %in% names(at)
  by_method <- at |>
    summarise(max_abs_diff = max(abs(diff)),
              max_abs_z = if (has_z) max(abs(z)) else NA_real_,
              max_abs_z_both = if (has_z) max(abs(z_both)) else NA_real_,
              .by = method)
  as.data.frame(x$methods) |>
    left_join(by_method, by = "method") |>
    mutate(from = paste(x$from_scale, collapse = ", "), to = as.character(x$to_scale), .before = 1) |>
    mutate(checks = worst_status(x$diagnostics))
}

#' @rdname as.data.frame.kits
#' @export
as.data.frame.kits_info <- function(x, row.names = NULL, optional = FALSE, ...) {
  data <- x$data
  data.frame(call = value_or_na(x$call, NA_character_),
             version = value_or_na(x$version, NA_character_),
             settings = settings_text(x$settings),
             n_rows = as.numeric(value_or_na(data$n_rows, NA_real_)),
             data_md5 = value_or_na(data$scores_md5, value_or_na(data$values_md5, NA_character_)),
             persons_md5 = value_or_na(data$persons_md5, NA_character_),
             input = value_or_na(x$input$name, NA_character_))
}

# 設定寫成一行字："cut = 2; n_boot = 100; ..."（向量用逗號接起來）
settings_text <- function(settings) {
  if (length(settings) == 0) {
    return(NA_character_)
  }
  values <- map_chr(settings, \(v) if (is.null(v)) "NULL" else paste(format(v), collapse = ", "))
  paste(names(settings), "=", values, collapse = "; ")
}

#' @rdname as.data.frame.kits
#' @export
as.data.frame.kits_item <- function(x, row.names = NULL, optional = FALSE, ...) {
  thr <- as.data.frame(x$thresholds)
  o <- x$fit$optim
  data.frame(scale = as.character(value_or_na(x$scale, NA_character_)),
             n_persons = x$n_persons,
             n_items = n_distinct(thr$item_id),
             n_steps = nrow(thr),
             n_patterns = n_distinct(x$design$pattern_id),
             n_fixed = sum(thr$beta_se == 0),
             n_left_out = n_distinct(value_or_na(x$left_out, character())),
             neg_loglik = x$fit$est$nll,
             convergence = as.integer(value_or_na(o$convergence, NA_integer_)),
             max_gradient = as.numeric(value_or_na(o$max_gradient, NA_real_)),
             checks = worst_status(x$diagnostics))
}

#' @rdname as.data.frame.kits
#' @export
as.data.frame.kits <- function(x, row.names = NULL, optional = FALSE, ...) {
  subj <- as.data.frame(x$subj)
  as.data.frame(x$item) |>
    mutate(theta_mean = mean(subj$theta),
           theta_sd = stats::sd(subj$theta),
           n_zero = sum(subj$flag == "zero"),
           n_perfect = sum(subj$flag == "perfect"),
           n_misfit = sum(subj$flag == "misfit"),
           checks = worst_status(x$diagnostics))
}

#' @rdname as.data.frame.kits
#' @export
as.data.frame.kits_norm_screen <- function(x, row.names = NULL, optional = FALSE, ...) {
  drop_table_class(x$waves)
}

#' @rdname as.data.frame.kits
#' @export
as.data.frame.kits_norm_fit <- function(x, row.names = NULL, optional = FALSE, ...) {
  m <- x$model
  data.frame(n = nrow(x$data),
             degrees = paste(x$degrees, collapse = ", "),
             neg_loglik = m$value,
             convergence = m$convergence,
             converged = m$converged,
             grad_max = m$grad_max,
             n_at_bound = length(m$at_bound))
}

#' @rdname as.data.frame.kits
#' @export
as.data.frame.kits_stats <- function(x, row.names = NULL, optional = FALSE, ...) {
  x$patterns |>
    mutate(n_items = lengths(items), items = map_chr(items, \(its) paste(its, collapse = ", "))) |>
    select(pattern_id, n_items, items, n, max_raw)
}
