#' Default percentiles of the norms
#'
#' The default `pct` of [kits_norms()] and [kits_norm_table()].
#'
#' @format A numeric vector: 3, 5, 10, 15, 25, 50, 75, 85, 90, 95, 97.
#' @examples
#' NORM_PCT
#' @export
NORM_PCT <- c(3, 5, 10, 15, 25, 50, 75, 85, 90, 95, 97)


#' Same-age z-score
#'
#' Looks up the shash parameters for each age in the norm table (interpolating between rows, each
#' `segment` on its own) and returns the same-age z; ages outside the norm range give NA.
#' A perfect score has no upper limit, so its z is not meaningful: set it to NA yourself (see the
#' examples). With a norm by group, each child is put on the curves of its group.
#'
#' Details: `kits_algorithm("norms")` and tutorial_algorithms/04_norms.html.
#'
#' @param d Scores (for example the D score), or a [kits_subj()] result (its `theta` and `age_m`
#'   are used; its scale must match the norm's, else `kits_error_scale`).
#' @param age_m Age in months (not with a kits_subj).
#' @param norm The result of [kits_norms()].
#'   Several segments can be stacked together, marked by a `segment` column.
#' @param by With scores as numbers and a norm by group: a data.frame of the group columns of the
#'   norm, one row per score (a [kits_subj()] needs these columns itself).
#' @param ... Not used: an argument the function does not have stops with an error (class
#'   `kits_error_input`).
#' @return With scores as numbers: the z-scores (the percentile rank is `100 * pnorm(z)`).
#'   With a [kits_subj()]: the scores with the columns `z` and `pct` (`100 * pnorm(z)`).
#' @examples
#' norm <- data.frame(age_m = c(6, 12), mu = c(40, 60), sigma = 5, epsilon = 0, delta = 1)
#' kits_norm_z(c(40, 50, 66), age_m = c(6, 9, 12), norm = norm)
#'
#' # Ages outside the range get no z
#' kits_norm_z(50, age_m = 20, norm = norm)
#'
#' # Perfect scores: no z
#' raw <- c(3, 5, 8)
#' max_raw <- c(8, 8, 8)
#' z <- kits_norm_z(c(40, 50, 66), age_m = c(6, 9, 12), norm = norm)
#' z[raw == max_raw] <- NA
#' z
#' @export
kits_norm_z <- function(d, ...) UseMethod("kits_norm_z")

#' @rdname kits_norm_z
#' @export
kits_norm_z.default <- function(d, age_m, norm, by = NULL, ...) {
  no_extra_args(..., fun = "kits_norm_z")
  if (is.null(attr(norm, "by")) && is.null(by)) {
    return(norm_z(d, age_m, norm))
  }
  abort_if(is.null(by), "The norm is by group (", paste(attr(norm, "by"), collapse = ", "),
           "): give the group of each ", "score in `by`, e.g. by = data.frame(sex = sex).")
  norm_z_by(d, age_m, norm, by_frame(by, length(d)))
}

#' @rdname kits_norm_z
#' @export
kits_norm_z.kits_subj <- function(d, norm, ...) {
  subj_extra_args(..., fun = "kits_norm_z", form = "kits_norm_z(s, norm)")
  s <- d
  if (missing(norm)) {
    kits_abort("Give the norm: kits_norm_z(s, norm).", class = "kits_error_input")
  }
  check_norm_scale(s, norm)
  by <- attr(norm, "by")
  z <- if (is.null(by)) {
    norm_z(s$theta, subj_column(s, "age_m"), norm)
  } else {
    norm_z_by(s$theta, subj_column(s, "age_m"), norm, by_columns(by, s))
  }
  out <- s
  out[["z"]] <- z                           # [[<- 留下類別、量尺和其他屬性
  out[["pct"]] <- 100 * stats::pnorm(z)
  attr(out, "norm_z") <- TRUE   # 給 `$<-`：以前 s$z <- kits_norm_z(s, norm) 的寫法，說清楚新的寫法
  out
}

#' @export
`$<-.kits_subj` <- function(x, name, value) {
  abort_if(inherits(value, "kits_subj") && isTRUE(attr(value, "norm_z")) && NROW(value) == NROW(x),
           "kits_norm_z(s, norm) returns the scores with the columns z and pct ",
           "(kits 0.22.0): write s <- s |> kits_norm_z(norm), not s$z <- kits_norm_z(s, norm).")
  NextMethod()
}

# 照組的常模：每個孩子用自己那一組的曲線；沒有常模的組是 NA（說一聲）
norm_z_by <- function(d, age_m, norm, groups) {
  by <- attr(norm, "by")
  if (is.null(by)) by <- names(groups)   # 自己疊的常模：by 的欄就是組的欄
  miss <- setdiff(by, names(groups))
  abort_if(length(miss) > 0, "The norm is by ", paste(by, collapse = ", "), "; the groups of the scores have no ",
           paste(miss, collapse = ", "), ".")
  table <- as.data.frame(norm)
  need_cols(table, by, "norm")
  keys <- group_keys(table[by])$keys
  groups <- need_groups(groups[by])
  id <- vctrs::vec_match(groups, vctrs::vec_cast(keys, groups))
  no_norm <- distinct(groups[is.na(id), , drop = FALSE])
  if (nrow(no_norm) > 0) {
    labels <- map_chr(seq_len(nrow(no_norm)), \(i) group_label(no_norm[i, , drop = FALSE]))
    kits_inform("No norm for group(s) ", paste(labels, collapse = "; "), ": their z is NA.",
                class = "kits_message_data")
  }
  z <- rep(NA_real_, length(d))
  for (i in seq_len(nrow(keys))) {
    rows <- which(id == i)
    curve <- table[vctrs::vec_match(table[by], keys[i, , drop = FALSE]) %in% 1,
                   setdiff(names(table), by)]
    if (length(rows) > 0) z[rows] <- norm_z(d[rows], age_m[rows], curve)
  }
  z
}

# 分數和常模在不同的量尺上時停下來（兩邊都有記錄時）
check_norm_scale <- function(s, norm) {
  norm_scale <- scale_of(norm)
  s_scale <- scale_of(s)
  need_same_scale(s_scale, norm_scale,
                  sprintf("The scores are on scale %s, but the norm is on scale %s. Convert them first with kits_convert().",
                          s_scale, norm_scale))
}

# 分數 d 和月齡 age_m：數字，一樣長（kits_norm_z()、kits_norms()、kits_norm_data() 用）
need_d_age <- function(d, age_m) {
  abort_if(length(d) != length(age_m), "`d` and `age_m` must have the same length.")
  abort_if(!is.numeric(d) || !is.numeric(age_m), "`d` and `age_m` must be numbers.")
}

# kits_subj 的方法多給的參數：age_m 來自 kits_subj（說清楚），其他照 no_extra_args()
subj_extra_args <- function(..., fun, form) {
  abort_if("age_m" %in% ...names(),
           sprintf("With a kits_subj, the ages come from its age_m column (the `persons` of kits_subj()): use %s ", form),
           "without `age_m`, or give the scores and ages as vectors (s$theta, age_m).")
  no_extra_args(..., fun = fun)
}

# 同齡 z（給 kits_norm_z() 的兩個方法）
norm_z <- function(d, age_m, norm) {
  need_cols(norm, c("age_m", "mu", "sigma", "epsilon", "delta"), "norm")
  abort_if(!all(map_lgl(norm[c("age_m", "mu", "sigma", "epsilon", "delta")], is.numeric)),
           "`norm$age_m`, `mu`, `sigma`, `epsilon`, `delta` must be numbers.")
  abort_if(any(!is.finite(norm$sigma) | norm$sigma <= 0) || any(!is.finite(norm$delta) | norm$delta <= 0),
           "`norm$sigma` and `norm$delta` must be positive numbers (the scale and the tail weight of the shash).")
  need_d_age(d, age_m)
  abort_if("segment" %in% names(norm) && anyNA(norm$segment), "`norm$segment` has NA.")
  segments <- norm_segments(norm)
  check_segments(segments)
  z <- rep(NA_real_, length(d))
  for (tab in segments) {
    hit <- !is.na(age_m) & age_m >= min(tab$age_m) & age_m <= max(tab$age_m)
    if (any(hit)) z[hit] <- segment_z(d[hit], age_m[hit], tab)
  }
  z
}

# 每一段至少 2 列、月齡不重複；段和段的月齡不重疊
check_segments <- function(segments) {
  for (i in seq_along(segments)) {
    tab <- segments[[i]]
    abort_if(nrow(tab) < 2 || anyDuplicated(tab$age_m),
             "Each norm segment needs at least 2 rows with different ages (segment ", names(segments)[i], ").")
    abort_if(i > 1 && min(tab$age_m) <= max(segments[[i - 1]]$age_m),
             "Norm segments overlap in age (", names(segments)[i - 1], " and ", names(segments)[i], ").")
  }
}

# 一段常模裡的同齡 z：shash 的四個參數在孩子的月齡上內插
segment_z <- function(D, age_m, tab) {
  at_age <- function(v) stats::approx(tab$age_m, tab[[v]], age_m)$y
  shash_z(D, at_age("mu"), at_age("sigma"), at_age("epsilon"), at_age("delta"))
}


# kits_subj 的一欄（例如 age_m，來自 persons）；沒有時說清楚
subj_column <- function(s, col) {
  abort_if(!col %in% names(s),
           "The kits_subj has no `", col, "` column: give `persons` with that column to kits_subj(), ",
           "or give the vectors instead (for example s$theta and age_m).")
  s[[col]]
}

# 常模用的分數不是 D 時，說一聲常模會在哪一條量尺上
norm_scale_message <- function(s) {
  sc <- scale_of(s)
  if (is.na(sc) || sc != "D") {
    kits_inform(sprintf(paste0("The norms will be on scale %s (not D). To have norms in D, ",
                               "convert the scores with kits_convert() first."),
                        scale_label(sc)), class = "kits_message_scale")
  }
}
