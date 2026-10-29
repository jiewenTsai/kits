#' Screening waves before fitting norms: waves with scores piled at the floor
#'
#' Computes, for each wave, the percentage of children at the lowest score of that wave (D rounded
#' to 3 decimals), and leaves a wave out of the norms when it reaches `cut`.
#' Waves piled at the ceiling still go into the norms.
#'
#' Details: `kits_algorithm("norm_screen")` and tutorial_algorithms/04_norms.html.
#'
#' @param d The score of each person-wave (for example the D score), or a [kits_subj()] result
#'   (its `theta`, `wave_m` and `cohort` are used; a message says when its scale is not `"D"`).
#' @param wave_m The wave (months) of each person-wave (with a kits_subj: taken from its column).
#' @param cohort The cohort of each person-wave; optional, the cohorts are then computed separately
#'   (with a kits_subj: taken from its column when present).
#' @param cut A wave whose percentage at the lowest score reaches this number is left out of the
#'   norms; the percentage is first rounded to 1 decimal (as printed in `floor_pct`).
#' @return A `kits_norm_screen`: a list that prints briefly. Its elements:
#'   * `waves`: one row per (cohort x) wave: `cohort` (when given), `wave_m`, `n`, `floor_pct` (the
#'     percentage at the lowest score), `normed` (`FALSE` means left out of the norms).
#'   * `normed`: one `TRUE`/`FALSE` per person-wave, in input order (NA in `d` is `FALSE`).
#'   * `cut`, `diagnostics` (see [kits_check()]) and `info` (see [kits_info()]).
#' @param ... Not used: an argument the function does not have stops with an error (class
#'   `kits_error_input`).
#' @examples
#' D <- c(10, 10, 10, 10, 12, 15,  20, 22, 25, 30, 31, 33)
#' wave_m <- rep(c(6, 12), each = 6)
#' screen <- kits_norm_screen(D, wave_m)
#' screen
#' screen$normed
#'
#' # The two cohorts computed separately
#' kits_norm_screen(D, wave_m, cohort = rep(c("M03", "M36"), times = 6))
#' @export
kits_norm_screen <- function(d, wave_m, cohort = NULL, cut = 20, ...) {
  no_extra_args(..., fun = "kits_norm_screen")
  call <- call_text(sys.call())
  d_expr <- substitute(d)
  input <- data_input(d, d_expr)
  if (inherits(d, "kits_subj")) {
    s <- d
    d <- s$theta
    if (missing(wave_m)) wave_m <- subj_column(s, "wave_m")
    if (is.null(cohort) && "cohort" %in% names(s)) cohort <- s$cohort
  }
  abort_if(length(d) != length(wave_m) || (!is.null(cohort) && length(cohort) != length(d)),
           "`d`, `wave_m` (and `cohort`) must have the same length.")
  abort_if(!is.numeric(d) && !all(is.na(d)), "`d` must be numbers.")
  need_number(cut, "cut", "one number (a percentage)", finite = FALSE)
  D <- round(d, 3)   # 比較分數前先捨入到小數點下 3 位（公式裡的 D）
  if (is.null(cohort)) cohort <- rep(NA, length(D))
  scores <- data.frame(D = D, wave_m = wave_m, cohort = cohort)
  ok <- is.finite(D) & !is.na(wave_m)
  abort_if(!any(ok), "No person-wave has both a finite `d` and a `wave_m`.", class = "kits_error_data")
  out <- scores |>
    filter(is.finite(D), !is.na(wave_m)) |>
    floor_shares() |>
    mutate(normed = floor_pct < cut)

  # 每個人次：所在的（世代 x）波次有沒有放進常模；分數是 NA 的人次是 FALSE
  in_normed_wave <- scores |>
    left_join(select(out, cohort, wave_m, normed), by = c("cohort", "wave_m")) |>
    pull(normed)
  person <- ok & in_normed_wave %in% TRUE
  if (all(is.na(out$cohort))) out$cohort <- NULL
  left <- filter(out, !normed)
  label <- if (is.null(left$cohort)) paste(left$wave_m) else paste(left$cohort, left$wave_m)
  any_left <- nrow(left) > 0
  left_status <- if (any_left) "note" else "ok"
  left_message <- if (any_left) sprintf("%g%% or more at the lowest score: %s", cut, paste(label, collapse = ", ")) else "none"
  diagnostics <- diag_table(list(diag_row("waves left out", nrow(left), left_status, left_message)))
  structure(list(waves = kits_table(out), normed = person, cut = cut, diagnostics = diagnostics,
                 info = make_info(call, list(cut = cut), values_fingerprint(D, wave_m, cohort), input = input)),
            class = "kits_norm_screen")
}


# 每個（世代 x）波次：人次，和拿到這一波最低分的百分比（先捨入到小數點下 1 位，和 kitsPreview 相同）。
# n 是 double，和 kitsPreview 的表相同。
floor_shares <- function(scores) {
  scores |>
    summarise(n = as.numeric(n()), floor_pct = round(100 * mean(D == min(D)), 1), .by = c(cohort, wave_m)) |>
    arrange(sort_rank(cohort), wave_m)
}

#' Print a norm screen
#'
#' @param x The result of [kits_norm_screen()].
#' @param digits Decimals to print (only the printout is rounded; the object keeps full precision).
#' @param ... Not used.
#' @return `x` (invisibly).
#' @examples
#' kits_norm_screen(c(10, 10, 10, 10, 12, 15, 20, 22, 25, 30, 31, 33), rep(c(6, 12), each = 6))
#' @export
print.kits_norm_screen <- function(x, digits = 2, ...) {
  need_digits(digits)
  cat(sprintf("<kits_norm_screen> waves with %g%% or more at the lowest score are left out\n", x$cut))
  print(format_table(x$waves, digits), row.names = FALSE, right = TRUE)
  cat(sprintf("Person-waves in the norms: %d of %d (use $normed)\n", sum(x$normed), length(x$normed)))
  checks_line(x$diagnostics)
  invisible(x)
}


#' Age range of the norms
#'
#' From the 1st percentile of age in the first wave to the 99th percentile of age in the last wave,
#' both ends rounded outward to `by` (by default half a month).
#' Pass only the waves to be normed (screen them first with [kits_norm_screen()]).
#'
#' Details: `kits_algorithm("norm_screen")` and tutorial_algorithms/04_norms.html.
#'
#' @param age_m Age in months of each person-wave, or a [kits_subj()] result (its `age_m` and
#'   `wave_m` are used, for example `kits_norm_ages(s[screen$normed, ])`).
#' @param wave_m The wave (months) of each person-wave (with a kits_subj: taken from its column).
#' @param probs Which percentile to take in the first wave and in the last wave.
#' @param by The age interval.
#' @return A vector of ages, which can be passed directly to `ages` of [kits_norms()].
#' @examples
#' set.seed(1)
#' age <- c(rnorm(200, 12, 0.6), 7, 8, rnorm(200, 24, 0.6))
#' wave_m <- c(rep(12, 202), rep(24, 200))
#' range(kits_norm_ages(age, wave_m))
#'
#' # The same from a kits_subj with age_m and wave_m columns
#' X <- kits_sim(rnorm(402), list(a1 = c(-1, 0), a2 = 0.5, a3 = -0.5), seed = 1)
#' persons <- data.frame(visit_id = rownames(X), age_m = age, wave_m = wave_m)
#' s <- kits_subj(X, kits_item(X), persons = persons)
#' range(kits_norm_ages(s))
#' @export
kits_norm_ages <- function(age_m, wave_m, probs = c(0.01, 0.99), by = 0.5) {
  if (inherits(age_m, "kits_subj")) {
    if (!missing(wave_m)) {
      kits_abort("With a kits_subj, wave_m comes from its column: kits_norm_ages(s), other arguments by name.",
                 class = "kits_error_input")
    }
    s <- age_m
    age_m <- subj_column(s, "age_m")
    wave_m <- subj_column(s, "wave_m")
  }
  abort_if(length(age_m) != length(wave_m), "`age_m` and `wave_m` must have the same length.")
  abort_if((!is.numeric(age_m) && !all(is.na(age_m))) || (!is.numeric(wave_m) && !all(is.na(wave_m))),
           "`age_m` and `wave_m` must be numbers (months).")
  need_number(by, "by", "one positive number", above = 0, finite = FALSE)
  ages <- data.frame(age_m = age_m, wave_m = wave_m) |>
    filter(is.finite(age_m), !is.na(wave_m))
  abort_if(nrow(ages) == 0, "No person-wave has both `age_m` and `wave_m`.", class = "kits_error_data")
  first <- filter(ages, wave_m == min(wave_m))$age_m
  last <- filter(ages, wave_m == max(wave_m))$age_m
  lo <- floor(stats::quantile(first, probs[1], names = FALSE) / by) * by
  hi <- ceiling(stats::quantile(last, probs[2], names = FALSE) / by) * by
  abort_if(hi <= lo, "The ages of the last wave are not above those of the first wave (", lo, " to ", hi, "). ",
           "Check `age_m` and `wave_m`.", class = "kits_error_data")
  round(seq(lo, hi, by = by), 10)
}
