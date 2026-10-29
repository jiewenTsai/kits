#' Print a kits table
#'
#' Tables from kits (for example [kits_subj()], [kits_itemfit()], [kits_subj_table()], `coef()` of
#' [kits_item()]) are data.frames that print `digits` decimals (2 by default; the object keeps its
#' full precision); a table that keeps diagnostics ends with one line, "Checks: ...".
#' A table of up to 30 rows prints in full, a longer one its first 10 rows (`print(x, n = 50)` prints
#' more, `as.data.frame(x)` gives all rows).
#' Selecting rows or columns and dplyr verbs keep the class and the recorded scale, and drop the
#' diagnostics and the provenance.
#' `rbind()` and `dplyr::bind_rows()` check that the recorded scales match: different scales stop
#' with an error of class `kits_error_scale`, a table without a recorded scale gives a warning of
#' class `kits_warning_scale`; `summary()` and [kits_check()] work on the stacked table.
#' Details: `vignette("kits-flow", package = "kits")`, section 4.
#'
#' @param x A table from kits.
#' @param digits Decimals to print.
#' @param n Rows to print. By default all rows of a table of up to 30 rows, and the first 10 rows
#'   of a longer table. `Inf` prints every row.
#' @param ... Passed to `print()` for a data.frame, for example `row.names = FALSE`.
#' @return `print()`: `x` (invisibly). `rbind()`: the stacked table.
#' @examples
#' betas <- list(a1 = c(-1, 0), a2 = c(0, 1), a3 = 0.5, a4 = -0.5)
#' X <- kits_sim(seq(-3, 3, length.out = 300), betas, seed = 1)
#' item <- kits_item(X)
#' coef(item)
#' print(coef(item), digits = 4)
#' coef(item)$thr[1]
#'
#' # A long table prints its first rows
#' s <- kits_subj(X, item)
#' s
#' print(s, n = 3)
#'
#' # rbind() keeps the scale when the scales match
#' nrow(rbind(s[1:5, ], s[6:10, ]))
#'
#' # dplyr verbs keep the class and the scale
#' library(dplyr)
#' s |>
#'   filter(flag == "ok") |>
#'   attr("scale")
#' @export
print.kits_table <- function(x, digits = 2, n = NULL, ...) {
  print_rows(x, digits, n, ...)
  checks_line(attr(x, "diagnostics"))
  invisible(x)
}

need_rows <- function(n) {
  need_number(n, "n", "one number of rows (0 or more; Inf for all rows)", lower = 0, finite = FALSE, null = TRUE)
}

# 印表格的列：30 列以內全印，再長只印前 10 列（或前 n 列），最後一行說還有幾列（內部用，給各個表格的 print）
# more：最後一行只說「... k more <more>」（例如 print.kits_stats() 的 "steps"）
print_rows <- function(x, digits, n = NULL, ..., more = NULL) {
  need_digits(digits)
  need_rows(n)
  total <- nrow(x)
  if (is.null(n)) n <- if (total <= 30) total else 10
  print(format_table(utils::head(as.data.frame(x), n), digits), right = TRUE, ...)
  if (total > n && !is.null(more)) {
    cat(sprintf("... %d more %s\n", total - n, more))
  } else if (total > n) {
    cat(sprintf("... %d more rows (more: print(x, n = %d); all rows: as.data.frame(x))\n", total - n, min(total, 50)))
  }
  invisible(x)
}

#' @rdname print.kits_table
#' @export
print.kits_itemfit <- function(x, digits = 2, n = NULL, ...) {
  type <- attr(x, "type")   # 第一行說是哪一種題目適配（舊的物件沒有記錄時不印）
  if (!is.null(type)) cat(itemfit_type_text(type), "\n", sep = "")
  print.kits_table(x, digits = digits, n = n, ...)
}

# kits_itemfit 的第一行
itemfit_type_text <- function(type) {
  how <- c(conditional = "given the raw score within each item pattern",
           wle = "at the WLE of each person-wave")
  sprintf("Item fit (%s): infit and outfit %s", type, how[[type]])
}

#' @rdname print.kits_table
#' @param i,j Rows and columns, as for a data.frame.
#' @export
`[.kits_table` <- function(x, i, j, ...) {
  out <- NextMethod()
  if (is.data.frame(out)) {
    keep <- c("scale", "beta_scale", "source", "rp", "links_without_se", "by",
              "items", "domain", "cohort", "raw_theta", "fitted", "type", "dif_labels")
    for (a in keep) attr(out, a) <- attr(x, a)
    if (!all(attr(out, "by") %in% names(out))) attr(out, "by") <- NULL   # 組的欄拿掉了：不再照組
    attr(out, "diagnostics") <- NULL   # 診斷和來源描述整張表，選出一部分後不再適用
    attr(out, "info") <- NULL
  }
  out
}

# dplyr 的動詞（filter()、mutate()、arrange()、join 等）做完之後：類別和量尺留著，
# 診斷和來源描述整張表，所以拿掉（和 `[.kits_table` 相同）
#' @export
dplyr_reconstruct.kits_table <- function(data, template) {
  out <- NextMethod()
  for (a in c("scale", "beta_scale", "rp")) {   # bind_rows()：第一張表沒有記錄時，用 vctrs 求出的共同量尺
    if (is.null(attr(template, a)) && !is.null(attr(data, a))) attr(out, a) <- attr(data, a)
  }
  # bind_rows() 只留第一張表的紀錄（dplyr 重建時用第一張表當樣板）：列變多、又有紀錄時說一聲
  has_record <- !is.null(attr(template, "links_without_se")) ||
    !is.null(attr(template, "raw_theta"))
  if (has_record && nrow(out) > nrow(template)) {
    kits_warn("dplyr::bind_rows() keeps the records of the first table only ",
              "(links without se_link, raw-score tables); ",
              "stack kits tables with rbind() or purrr::list_rbind(), which keep them all.",
              class = "kits_warning_scale")
  }
  attr(out, "diagnostics") <- NULL
  attr(out, "info") <- NULL
  out
}

#' @rdname print.kits_table
#' @param deparse.level Passed to `rbind()`.
#' @export
rbind.kits_table <- function(..., deparse.level = 1) {
  parts <- compact(list(...))
  is_kits <- map_lgl(parts, \(p) inherits(p, "kits_table"))
  scales <- map_chr(parts, \(p) if (inherits(p, "kits_table")) scale_of(p) else NA_character_)
  recorded <- unique(scales[!is.na(scales)])
  abort_if(length(recorded) > 1, "The tables are on different scales: ", paste(recorded, collapse = " and "),
           ". Convert them to one scale with kits_convert() first.", class = "kits_error_scale")
  unchecked <- !is_kits & map_lgl(parts, is.data.frame)
  if (length(recorded) == 1) unchecked <- unchecked | (is_kits & is.na(scales))
  if (any(unchecked)) {
    kits_warn(sum(unchecked), " table(s) without a recorded scale bound to tables on scale ", scale_label(recorded),
              ": their scale was not checked.", class = "kits_warning_scale")
  }
  # 先拿掉 kits 的類別和屬性，用 data.frame 的 rbind 疊起來（列名、因子的規則和 base 的 rbind 相同），
  # 再加回類別和共同的量尺
  plain <- map(parts, plain_table)
  # 連結誤差的欄（kits_convert() 的 se_link、se_total）只有部分的表有時，其他的表補 NA
  extra <- intersect(c("se_link", "se_total"), unlist(map(plain, \(p) if (is.data.frame(p)) names(p))))
  plain <- map(plain, \(p) add_na_columns(p, extra))
  out <- do.call(rbind, c(plain, list(deparse.level = deparse.level)))
  class(out) <- class(parts[[which(is_kits)[1]]])
  if (length(recorded) == 1) attr(out, "scale") <- recorded
  beta_scales <- unique(map_chr(parts[is_kits], beta_scale_of))
  if (length(recorded) == 1 && length(beta_scales) == 1 && !identical(beta_scales, recorded)) attr(out,
      "beta_scale") <- beta_scales
  rps <- unique(map_dbl(parts[is_kits], rp_of))
  abort_if(length(rps) > 1,
           "The tables hold thresholds at different response probabilities (RP): ", paste(100 * rps, collapse = " and "),
           ".")
  if (length(rps) == 1 && rps != 0.5) attr(out, "rp") <- rps
  attr(out, "links_without_se") <- combine_no_se_records(parts)   # 有一張表的連結誤差沒有帶下來：疊起來也沒有
  attr(out, "raw_theta") <- combine_raw_theta(parts)                # 每張分數表的原始分 -> theta 的來源
  attr(out, "diagnostics") <- stacked_diagnostics(out, parts[is_kits])
  out
}

#' @rdname print.kits_table
#' @export
rbind.kits_subj <- rbind.kits_table

# 疊起來的表的診斷：kits_subj 從疊好的表重算；其他的表把各表的診斷疊起來，檢查的名字前面加 "table i: "
# （各表的診斷還是描述它自己的列）。沒有任何診斷時是 NULL。
stacked_diagnostics <- function(out, parts) {
  if (inherits(out, "kits_subj") && all(c("theta", "se", "flag") %in% names(out)) && nrow(out) > 0) {
    first <- diag_row("stacked tables", length(parts), "ok",
                      sprintf("diagnostics recomputed for the %d person-waves of %d stacked tables", nrow(out), length(parts)))
    return(kits_table(rbind(first, as.data.frame(subj_diagnostics(out))), "kits_diagnostics"))
  }
  each <- map(seq_along(parts), \(i) {
    d <- attr(parts[[i]], "diagnostics")
    if (is.null(d)) {
      return(NULL)
    }
    d <- as.data.frame(d)
    d$check <- paste0("table ", i, ": ", d$check)
    d
  })
  each <- compact(each)
  if (length(each) == 0) {
    return(NULL)
  }
  kits_table(list_rbind(each), "kits_diagnostics")
}

# 只留 data.frame 的屬性（欄名、列名、類別）；不是 data.frame 的照原樣
plain_table <- function(p) {
  if (!is.data.frame(p)) {
    return(p)
  }
  p <- as.data.frame(p)
  attributes(p) <- attributes(p)[c("names", "row.names", "class")]
  p
}

# 缺的欄補 NA（給 rbind() 的 se_link、se_total）
add_na_columns <- function(p, columns) {
  if (!is.data.frame(p)) {
    return(p)
  }
  for (v in setdiff(columns, names(p))) p[[v]] <- rep(NA_real_, nrow(p))
  p
}
