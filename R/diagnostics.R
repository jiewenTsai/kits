# 診斷表：每個主要結果記下檢查的結果（list 放在 x$diagnostics，data.frame 放在 attr(x, "diagnostics")）。
# 欄位：check、value、status（"ok"、"note"、"warn"）、message。

diag_row <- function(check, value, status, message) {
  data.frame(check = check, value = as.numeric(value), status = status, message = message, stringsAsFactors = FALSE)
}

# 一個是或否的檢查：bad 時是 status 和 bad_msg，否則 "ok" 和 ok_msg（只算用到的那個訊息）
diag_flag <- function(check, value, bad, status, bad_msg, ok_msg) {
  diag_row(check, value, if (bad) status else "ok", if (bad) bad_msg else ok_msg)
}

diag_table <- function(rows) {
  kits_table(list_rbind(unname(rows)), "kits_diagnostics")
}

# 物件的診斷。kits_subj 疊起來或選過之後沒有診斷：從表裡重算（0 分、滿分、不適配的比例，每個波次的信度）
diagnostics_of <- function(x) {
  d <- if (is.data.frame(x)) attr(x, "diagnostics") else if (is.list(x) && !is.null(names(x))) x$diagnostics
  if (is.null(d) && inherits(x, "kits_subj") && all(c("theta", "se", "flag") %in% names(x)) && nrow(x) > 0) {
    d <- recomputed_subj_diagnostics(x)
  }
  d
}

recomputed_subj_diagnostics <- function(x) {
  d <- subj_diagnostics(x)
  first <- diag_row("diagnostics recomputed", nrow(x), "ok",
                    sprintf("for these %d person-waves (the table was stacked or subset)", nrow(x)))
  kits_table(rbind(first, as.data.frame(d)), "kits_diagnostics")
}

plural <- function(n, word) sprintf("%d %s%s", n, word, if (n == 1) "" else "s")

# 印在每個結果最後的一行
checks_line <- function(diag) {
  if (is.null(diag)) {
    return(invisible())
  }
  n_warn <- sum(diag$status == "warn")
  n_note <- sum(diag$status == "note")
  if (n_warn + n_note == 0) {
    cat("Checks: OK\n")
  } else {
    cat(sprintf("Checks: %s (see summary() or kits_check())\n",
                paste(c(if (n_warn > 0) plural(n_warn, "warning"), if (n_note > 0) plural(n_note, "note")), collapse = ", ")))
  }
  invisible()
}

# 印診斷表：一個檢查一行；整張表照 object、status、check、value、message，select() 之後照選的欄；value 印 3 位有效數字，
# message 放最後、不補空白
#' @export
print.kits_diagnostics <- function(x, ...) {
  d <- as.data.frame(x)
  if (ncol(d) == 0 || nrow(d) == 0) {
    return(invisible(x))
  }
  std <- c("object", "status", "check", "value", "message")   # 整張表：照這個順序
  if (all(std[-1] %in% names(d))) d <- d[c(intersect(std, names(d)), setdiff(names(d), std))]
  cols <- map(d, \(v) if (is.numeric(v)) diag_number(v) else ifelse(is.na(v), "", as.character(v)))
  last <- if (names(cols)[length(cols)] == "message") length(cols) else 0
  padded <- imap(cols, \(v, nm) {
    v <- c(nm, v)
    if (match(nm, names(cols)) == last) return(v)
    formatC(v, width = if (is.numeric(d[[nm]])) max(nchar(v)) else -max(nchar(v)))
  })
  lines <- pmap_chr(unname(padded), paste, sep = "  ")   # 一列一個檢查，欄和欄之間兩個空格
  cat(lines, sep = "\n")
  invisible(x)
}

# 診斷表的數字：3 位有效數字，NA 空白
diag_number <- function(v) ifelse(is.na(v), "", formatC(signif(v, 3), format = "g", digits = 3))


#' Diagnostics of one or more results
#'
#' Every main result of kits checks itself when it is made and keeps a diagnostics table
#' (`x$diagnostics` for lists, `attr(x, "diagnostics")` for tables), one row per check; its
#' printout ends with "Checks: OK" or the number of warnings and notes. `kits_check()` collects
#' the tables of several results; an object without diagnostics gets a row saying so (a stacked
#' [kits_subj()] gets its diagnostics recomputed).
#'
#' The checks:
#' * [kits_item()]: convergence, gradient, Hessian, extreme estimates, category counts, item
#'   pattern sizes (note below 400 person-waves) and item fit (infit or outfit outside 0.6 to 1.4).
#' * [kits_subj()]: shares of zero, perfect and misfit scores, and the reliability of each wave.
#' * [kits_subj_table()]: person-waves who answered items outside the standard form of their wave.
#' * [kits_itemfit()]: whether the item parameters were calibrated on these data.
#' * Links ([kits_linking()], [kits_link_steps]): points kept and dropped, the correlation of the
#'   kept points, the linking error, slope vs 1, the common-item share and the domain balance.
#' * [kits_link_check()]: a check method that differs from the link by more than 2 errors
#'   (or `tol`).
#' * [kits_item_table()] with `data`: items separated by age.
#' * [kits_norms()]: convergence, gradient, coefficients at their limits, the largest delta, and
#'   the same-age z by age decile.
#' * [kits_norm_screen()]: waves left out. [kits_drift()]: drifted items.
#' * [kits_dif()]: the LR test, the largest step difference, the overall DIF test, the largest |dif|,
#'   items not compared, steps not observed in a group.
#' * [kits_vertical()]: scale shrinkage, falling means, and cohorts that differ at a shared wave.
#' * [kits_precision()]: the usable range of each form.
#'
#' Details: `vignette("kits-flow", package = "kits")`, section 6.
#'
#' @section Quality control for a new wave:
#' A short checklist, after Kolen and Brennan (2014; Table 8.8, p. 336, and p. 333):
#' 1. Drift check before fixing: [kits_drift()] against the original calibration (pp. 288-289).
#' 2. Fix to the original calibration: `kits_item(X_new, fixed_params = drift$fix)` (pp. 293-294).
#' 3. Double link: a second route to the target scale ([kits_link_check()]; pp. 295-296).
#' 4. Item pattern n: about 400 person-waves per form (p. 304; the check of [kits_item()]).
#' 5. Common-item share: at least 20% of the items, a content miniature ([kits_linking()]).
#' 6. Consistency with the other cohort at the same age (p. 333; [kits_vertical()]).
#' 7. md5 provenance: compare the fingerprints of [kits_info()] (p. 333).
#' @references
#' Kolen, M. J., & Brennan, R. L. (2014). *Test equating, scaling, and linking* (3rd ed.). Springer.
#' \doi{10.1007/978-1-4939-0317-7}
#'
#' @param ... Results of kits, named or not (the names label the rows).
#' @param all `FALSE` (the default) keeps only notes and warnings; `TRUE` keeps every check.
#' @return A table (a data.frame), warnings first: `object`, and the diagnostics columns `check`,
#'   `value`, `status` (`"ok"`, `"note"` or `"warn"`) and `message`.
#' @examples
#' betas <- list(a1 = c(-1, 0), a2 = c(0, 1), a3 = 0.5, a4 = -0.5)
#' X <- kits_sim(seq(-3, 3, length.out = 300), betas, seed = 1)
#' item <- kits_item(X)
#' s <- kits_subj(X, item)
#' kits_check(item, s)
#' kits_check(calibration = item, all = TRUE)
#' summary(item)
#' @export
kits_check <- function(..., all = FALSE) {
  objs <- list(...)
  abort_if(length(objs) == 0, "Give one or more results of kits, e.g. kits_check(item, s).")
  exprs <- as.list(substitute(list(...)))[-1]
  nms <- names(objs)
  labels <- map_chr(seq_along(objs), \(i) {
    if (!is.null(nms) && nms[i] != "") nms[i] else object_label(exprs[[i]], i)
  })
  out <- list_rbind(map2(objs, labels, check_rows))
  if (!isTRUE(all)) out <- filter(out, status != "ok")
  out <- arrange(out, match(status, c("warn", "note", "ok")))
  kits_table(out, c("kits_check", "kits_diagnostics"))
}

# kits_check() 的列的名字：物件的名字；不是名字時（例如 do.call() 給的值）寫 "object i"（不帶資料）
object_label <- function(expr, i) {
  if (is.symbol(expr)) as.character(expr) else call_text_or_index(expr, i)
}

call_text_or_index <- function(expr, i) {
  if (is.call(expr)) {
    txt <- call_text(expr)
    if (!grepl("<data>", txt, fixed = TRUE)) return(txt)
  }
  sprintf("object %d", i)
}

# 一個物件的診斷列。沒有診斷時不停下來，給一列說明（這樣其他物件的檢查還在）：
# 疊起來或選過的 kits 表（rbind()、bind_rows()、dplyr 的動詞）給註記；kits_subj 會重算；不是 kits 的結果給警告
check_rows <- function(obj, label) {
  d <- diagnostics_of(obj)
  if (is.null(d)) d <- missing_diagnostics(obj)
  data.frame(object = label, as.data.frame(d))
}

# 沒有診斷的物件的說明列
missing_diagnostics <- function(obj) {
  if (inherits(obj, "kits_table")) {
    return(diag_table(list(diag_row("diagnostics", NA, "note", paste(
      "none kept: the table was stacked (rbind(), bind_rows()) or subset (x[i, ], dplyr verbs), and the diagnostics",
      "described the whole table; check the tables before stacking, e.g. kits_check(sA, sB)")))))
  }
  what <- if (is.object(obj)) class(obj)[1] else if (is.matrix(obj)) "matrix" else typeof(obj)
  diag_table(list(diag_row("diagnostics", NA, "warn",
                           sprintf("not a result of a kits function (a %s), so it has no checks; leave it out", what))))
}

#' @export
print.kits_check <- function(x, ...) {
  if (nrow(x) == 0) {
    cat("Checks: OK (no notes or warnings)\n")
    return(invisible(x))
  }
  NextMethod()
}


#' Summary of a kits result: a short print and the diagnostics
#'
#' Prints a result with its whole diagnostics table. For a link, `theta` adds how much the
#' screening moves the converted values (the line before and after screening at `theta`); for
#' scores converted by [kits_convert()], the summary also shows the links applied, the line they
#' make together and, when every link has draws, the linking error at the mean theta.
#'
#' Details: tutorial_algorithms/05_checks.html.
#'
#' @param object A result of kits: [kits_item()], [kits_subj()], [kits_itemfit()],
#'   [kits_subj_table()], [kits()], a link ([kits_linking()], [kits_link_line()]), [kits_drift()],
#'   [kits_norm_screen()], [kits_norms()], [kits_dif()], [kits_vertical()], [kits_precision()] or
#'   [kits_link_check()].
#' @param digits Decimals to print.
#' @param theta For a link: optional values on the scale the link starts from, for example the
#'   median theta of each wave.
#' @param ... Not used.
#' @return A `kits_summary` (printed): `object` and `diagnostics`, and
#'   * for a link with `theta`: `screening`, a table with `theta`, `before` (intercept0 + slope0
#'     theta), `after` (intercept + slope theta), `shift`, and with draws `se_link` and
#'     `z` = shift / se_link, on the scale where the link ends (`attr(, "scale")`);
#'   * for converted scores: `conversion`, a list with `chain` (one row per link: `step`, `from`,
#'     `to`, `intercept`, `slope`), the composed `intercept` and `slope`, `from_scale`, `to_scale`,
#'     `theta_mean`, `value_mean` and `se_link` (NA when a link has no draws); one list per
#'     calibration for stacked scores.
#' @examples
#' betas <- list(a1 = c(-1, 0), a2 = c(0, 1), a3 = 0.5, a4 = -0.5)
#' X <- kits_sim(seq(-3, 3, length.out = 300), betas, seed = 1)
#' summary(kits_item(X))
#' summary(kits_subj(X, kits_item(X)))
#'
#' # Converted scores: the links applied, as one line, and the linking error at the mean
#' item <- kits_item(X, scale = "x")
#' anchors <- data.frame(item_id = c("a1", "a2", "a3", "a4"), tau = c(42, 50, 54, 46))
#' d_link <- kits_linking(item, anchors, n_boot = 50)
#' s <- kits_subj(X, item) |>
#'   kits_convert(d_link)
#' summary(s)
#' summary(s)$conversion$se_link
#'
#' # A link: how much the screening moves the converted values
#' points <- data.frame(id = paste0("p", 1:8), item = paste0("p", 1:8),
#'                      from = c(-2, -1, 0, 1, 2, 3, 0.5, 1.5),
#'                      to = c(-1, 0, 1, 2, 3, 4, 1.5, 5.5))
#' link <- points |>
#'   kits_link_screen() |>
#'   kits_link_line() |>
#'   kits_link_error(n_boot = 100)
#' summary(link, theta = c(-1, 0, 1, 2))
#' summary(link, theta = c(-1, 0, 1, 2))$screening
#' @name summary.kits
NULL

kits_summary <- function(object, digits = 2, ...) {
  need_digits(digits)
  structure(list(object = object, diagnostics = diagnostics_of(object), digits = digits), class = "kits_summary")
}

#' @rdname summary.kits
#' @export
summary.kits_item <- kits_summary
#' @rdname summary.kits
#' @export
summary.kits_table <- function(object, digits = 2, ...) {
  out <- kits_summary(object, digits)
  out$conversion <- conversion_summary(object)   # 換算過的分數：一串連結合成的線（R/line.R）
  out
}
#' @rdname summary.kits
#' @export
summary.kits <- kits_summary
#' @rdname summary.kits
#' @export
summary.kits_link <- function(object, digits = 2, theta = NULL, ...) {
  out <- kits_summary(object, digits)
  if (!is.null(theta)) out$screening <- screening_shift(object, theta)
  out
}
#' @rdname summary.kits
#' @export
summary.kits_drift <- kits_summary
#' @rdname summary.kits
#' @export
summary.kits_norm_screen <- kits_summary

#' @rdname summary.kits
#' @param x A `kits_summary`.
#' @export
print.kits_summary <- function(x, ...) {
  obj <- x$object
  if (is.data.frame(obj)) {
    cat(table_header(obj, x$digits), sep = "\n")
  } else {
    print(obj, digits = x$digits)
  }
  if (!is.null(x$screening)) {
    cat("Screening: the line before (all points) and after, at theta (shift = after - before):\n")
    print(format_table(x$screening, x$digits), row.names = FALSE, right = TRUE)
  }
  if (!is.null(x$conversion)) print_conversion(x$conversion, x$digits)
  if (is.null(x$diagnostics)) {
    cat("No diagnostics recorded.\n")
  } else {
    cat("Diagnostics:\n")
    print(x$diagnostics)
  }
  invisible(x)
}

# 表格的簡短說明（summary() 用；kits_subj 的 print 也用）
table_header <- function(x, digits = 2) {
  if (inherits(x, "kits_subj")) {
    flags <- table(factor(x$flag, levels = c("ok", "misfit", "zero", "perfect")))
    th <- x$theta
    return(c(sprintf("<kits_subj> %d person-waves; scale: %s", nrow(x), scale_label(scale_of(x))),
             if (length(th) > 0) sprintf("  theta mean %s, SD %s, range %s to %s; flags: %s", fmt(mean(th), digits),
                                         fmt(stats::sd(th), digits), fmt(min(th), digits), fmt(max(th), digits),
                                         paste(names(flags), flags, sep = " ", collapse = ", "))))
  }
  if (inherits(x, "kits_norms")) {
    return(sprintf("<kits_norms> ages %s to %s months (%d rows); scale: %s", format(min(x$age_m)), format(max(x$age_m)),
                   nrow(x), scale_label(scale_of(x))))
  }
  if (inherits(x, "kits_itemfit")) {
    return(sprintf("<kits_itemfit> %d items; infit %s to %s, outfit %s to %s", nrow(x), fmt(min(x$infit, na.rm = TRUE), digits),
                   fmt(max(x$infit, na.rm = TRUE), digits), fmt(min(x$outfit, na.rm = TRUE), digits),
                   fmt(max(x$outfit, na.rm = TRUE), digits)))
  }
  if (inherits(x, "kits_dif")) {
    g <- attr(x, "groups")
    if (is.null(g)) return(sprintf("<kits_dif> %d groups; no item DIF", NROW(attr(x, "lrtest")$groups)))
    return(sprintf("<kits_dif> %d items; %s minus %s; dif %s to %s", nrow(x), g$ref, g$focal, fmt(min(x$dif), digits),
                   fmt(max(x$dif), digits)))
  }
  if (inherits(x, "kits_subj_table")) {
    return(sprintf("<kits_subj_table> standard form of %d wave(s), %d rows; scale: %s", length(unique(x$wave_m)), nrow(x),
                   scale_label(scale_of(x))))
  }
  if (inherits(x, "kits_vertical")) {
    return(sprintf("<kits_vertical> %d cohort x wave group(s); scale: %s", nrow(x), scale_label(scale_of(x))))
  }
  if (inherits(x, "kits_precision")) {
    return(sprintf("<kits_precision> %d form(s); RMSE criterion %s; scale: %s", length(unique(x$form)),
                   format(signif(attr(x, "criterion"), 3)), scale_label(scale_of(x))))
  }
  sprintf("<%s> %d rows", class(x)[1], nrow(x))
}
