# kits_norms() 的三步，一步一個匯出的函數（讓學生一步一步跑；kits_algorithm(norm) 列出它們）：
#   kits_norm_data()：資料 -> kits_norm_fit()：shash 的最大概似 -> kits_norm_table()：常模表
# 和 kits_norms() 是同一段程式（norm_data_step()、norm_fit_step()、norm_table_step()），所以數字相同。

#' A norm step by step: data, fit, table
#'
#' The three steps of [kits_norms()], each a function; with the same settings they give the same
#' table:
#' 1. `kits_norm_data()`: the person-waves with finite score and age within `ages` plus `pad`;
#' 2. `kits_norm_fit()`: the shash by maximum likelihood, one fit per start of delta, the best kept;
#' 3. `kits_norm_table()`: the shash parameters and percentiles at each age, with the diagnostics.
#'
#' Details: `kits_algorithm("norms")` and tutorial_algorithms/04_norms.html.
#'
#' @param d The score of each person-wave (for example the D score). For a [kits_subj()] result `s`,
#'   give `s$theta` and `s$age_m`.
#' @param age_m Age in months of each person-wave.
#' @param ages The ages at which the norms will be given.
#' @param pad Use data up to this many months outside the range of `ages`.
#' @param age_x The transformation function for age.
#' @param ... In `kits_norm_data()`: not used (an argument the function does not have stops with an
#'   error, class `kits_error_input`). In `print()`: not used.
#' @param norm_data A `kits_norm_data()` result.
#' @param degrees The degrees of the four polynomials for mu, sigma, epsilon, delta (see
#'   [kits_norms()]); at least `10 * sum(degrees + 1)` person-waves are needed (`kits_error_data`).
#' @param delta_starts Start values for delta.
#' @param fit A `kits_norm_fit()` result.
#' @param pct The percentiles to output.
#' @param delta_max When the largest delta exceeds this number, a warning is shown.
#' @param x A `kits_norm_data` or a `kits_norm_fit`.
#' @param digits Decimals to print (only the printout is rounded).
#' @param n Rows to print (see [print.kits_table()]).
#' @return
#' * `kits_norm_data()`: a `kits_norm_data`, a data.frame with `d` and `age_m`, one row per
#'   person-wave used; `attr(, "norm")` keeps `ages`, `ages_in`, `pad`, `age_x`, `age_x_text`,
#'   `scale`.
#' * `kits_norm_fit()`: a `kits_norm_fit`, a list with `model` (the polynomial coefficients, -log L
#'   `value`, convergence), `data` (the `kits_norm_data()`), `degrees` and `delta_starts`.
#' * `kits_norm_table()`: a `kits_norms` (see [kits_norms()]).
#' @seealso [kits_norms()] (all three steps in one call), [kits_norm_z()].
#' @examples
#' set.seed(1)
#' age <- runif(600, 3, 36)
#' D <- 20 + 15 * log(age) + rnorm(600, 0, 3)
#'
#' # 1 the data inside the ages
#' norm_data <- kits_norm_data(D, age, ages = seq(3, 36, by = 3))
#' norm_data
#'
#' # 2 the shash by maximum likelihood
#' fit <- kits_norm_fit(norm_data, degrees = c(3, 2, 1, 1))
#' fit
#'
#' # 3 the norm table
#' norm <- kits_norm_table(fit)
#' norm
#' @name kits_norm_steps
NULL

#' @rdname kits_norm_steps
#' @export
kits_norm_data <- function(d, age_m, ages, pad = 1, age_x = sqrt, ...) {
  no_extra_args(..., fun = "kits_norm_data")
  age_x_text <- age_x_code(substitute(age_x), age_x)
  abort_if(inherits(d, "kits_subj"),
           "kits_norm_data() takes the scores and ages as vectors: kits_norm_data(s$theta, s$age_m, ages). ",
           "kits_norms(s, ages) takes a kits_subj directly.")
  need_d_age(d, age_m)
  norm_data_step(d, age_m, ages, pad, age_x, age_x_text, scale = NA_character_)
}

#' @rdname kits_norm_steps
#' @export
print.kits_norm_data <- function(x, digits = 2, n = NULL, ...) {
  setting <- attr(x, "norm")
  cat(sprintf("<kits_norm_data> %d person-waves, ages %s to %s (pad %s); the table at %d ages from %s to %s\n",
              nrow(x), format(setting$ages_in[1]), format(setting$ages_in[2]), format(setting$pad),
              length(setting$ages), format(min(setting$ages)), format(max(setting$ages))))
  print_rows(x, digits, n)
  cat("Next: kits_norm_fit(x).\n")
  invisible(x)
}


#' @rdname kits_norm_steps
#' @export
kits_norm_fit <- function(norm_data, degrees = c(5, 4, 3, 2), delta_starts = c(0.5, 1, 1.5, 2)) {
  need_norm_data(norm_data)
  check_norm_fit_settings(degrees, delta_starts)
  norm_fit_step(norm_data, degrees, delta_starts)
}

#' @rdname kits_norm_steps
#' @export
print.kits_norm_fit <- function(x, ...) {
  m <- x$model
  cat(sprintf("<kits_norm_fit> shash by maximum likelihood: %d person-waves, degrees %s, -log L = %s\n",
              nrow(x$data), paste(x$degrees, collapse = ", "), format(m$value, digits = 10)))
  cat(sprintf("  %s\n", norm_convergence_text(m)))
  if (length(m$at_bound) > 0) cat(sprintf("  coefficient(s) at their limits: %s\n", paste(m$at_bound, collapse = ", ")))
  cat("Next: kits_norm_table(x).\n")
  invisible(x)
}

# 一定要是 kits_norm_data() 的結果（設定在 attr(, "norm")）
need_norm_data <- function(norm_data) {
  abort_if(!inherits(norm_data, "kits_norm_data") || is.null(attr(norm_data, "norm")),
           "`norm_data` must be a kits_norm_data() result, e.g. kits_norm_fit(kits_norm_data(D, age_m, ages)).")
}


#' @rdname kits_norm_steps
#' @export
kits_norm_table <- function(fit, pct = NORM_PCT, delta_max = 3) {
  abort_if(!inherits(fit, "kits_norm_fit"), "`fit` must be a kits_norm_fit() result.")
  norm_table_step(fit, pct, delta_max, call_text(sys.call()))
}
