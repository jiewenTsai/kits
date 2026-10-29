# 量尺上每個位置的精準度（Kolen & Brennan, 2014, pp. 405, 413, 463-465）：
# 每個題本，在 theta 的格點上，原始分的分佈 P(r | theta)（Lord-Wingersky 的遞迴，多點計分版，p. 232, Eq. 6.41），
# 原始分換成 WLE（和 kits_subj_table() 同一個 WLE），再算 E[theta_hat | theta]、偏誤、SD、RMSE。

# 題本：kits_subj_table（每個波次的標準題本）、有名字的 list（每個元素一組題號）、或一組題號（內部用）
as_forms <- function(forms, parms) {
  if (inherits(forms, "kits_subj_table")) {
    standard <- filter(attr(forms, "forms"), type == "standard")
    out <- strsplit(standard$items, ", ", fixed = TRUE)
    names(out) <- paste("wave", wave_label(standard$wave_m))
  } else if (is.character(forms)) {
    out <- list(form = forms)
  } else if (is.list(forms) && length(forms) > 0 && all(map_lgl(forms, is.character))) {
    out <- forms
    if (is.null(names(out))) names(out) <- paste("form", seq_along(out))
    names(out)[names(out) == ""] <- paste("form", which(names(out) == ""))
  } else {
    kits_abort("`forms` must be a kits_subj_table() result, ",
               "a named list of item patterns (character vectors), ",
               "or one item pattern.", class = "kits_error_input")
  }
  abort_if(anyDuplicated(names(out)), "The forms need different names.")
  for (f in names(out)) {
    its <- unique(out[[f]])
    if (length(its) == 0) kits_abort("Form ", f, " has no items.", class = "kits_error_input")
    need_params(its, parms$item_id, paste0("Form ", f, " has item(s) without parameters in `item`: "))
    out[[f]] <- sort(its)
  }
  out
}

# 一串連結（一條，或 list）；檢查都是連結
as_links <- function(link) {
  if (is.null(link)) {
    return(list())
  }
  if (inherits(link, "kits_link")) link <- list(link)
  abort_if(!is.list(link) || !all(map_lgl(link, \(l) inherits(l, "kits_link"))),
           "`link` must be a link (kits_linking(), kits_link_line()) or a list of links, ",
           "in the order they are applied.")
  link
}

# 最長的一段連續格點（含 RMSE 最小的點）：RMSE < criterion
usable_run <- function(ok, best) {
  if (!ok[best]) {
    return(c(NA_integer_, NA_integer_))
  }
  lo <- best
  hi <- best
  while (lo > 1 && ok[lo - 1]) lo <- lo - 1
  while (hi < length(ok) && ok[hi + 1]) hi <- hi + 1
  c(lo, hi)
}


#' Precision of each form across the scale
#'
#' For each form (a set of items), how well does the WLE recover theta at each point of the scale?
#' At each theta of a grid it gives the expected WLE, its bias, SD and RMSE, from the raw-score
#' distribution of the Lord-Wingersky recursion (Thissen et al., 1995; Kolen & Brennan, 2014),
#' and the usable range of each form, where RMSE < `criterion`.
#' Details: `kits_algorithm("information")` and tutorial_algorithms/06_scale.html.
#'
#' @param item The item parameters: a [kits_item()] result (its step parameters `beta`, on its calibration scale),
#'   or a [kits()], or a data.frame with `item_id`, `item_score`, `beta`.
#' @param forms The forms: a [kits_subj_table()] result (the standard form of each wave, from `attr(, "forms")`),
#'   a named list of item patterns (character vectors of item ids), or one item pattern.
#' @param theta The grid, on the calibration scale of `item` (logits).
#' @param link Optional: a link or a list of links (applied in order) to report the results on another scale, for example
#'   `list(co, d_link)` for D.
#' @param criterion The largest RMSE of the usable range, on the scale of the result.
#'   `NULL` (the default) means 1 logit of the calibration scale, that is the product of the slopes
#'   of `link`.
#' @return A table (a data.frame) with one row per form and theta: `form`, `theta`, `expected` (E\[theta-hat | theta\]),
#'   `bias`, `sd`, `rmse`; on the scale in `attr(, "scale")`.
#'   `attr(, "ranges")`: one row per form, `form`, `n_items`, `max_raw`, `from`, `to` (the usable range; NA when no
#'   point of the grid reaches the criterion), `min_rmse` and `at` (where it is smallest).
#'   `attr(, "criterion")`, diagnostics and [kits_info()] too. `plot()` draws the RMSE curves.
#' @references
#' Kolen, M. J., & Brennan, R. L. (2014). *Test equating, scaling, and linking* (3rd ed.). Springer.
#' \doi{10.1007/978-1-4939-0317-7}
#'
#' Thissen, D., Pommerich, M., Billeaud, K., & Williams, V. S. L. (1995). Item response theory for scores on tests
#' including polytomous items with ordered responses. *Applied Psychological Measurement, 19*(1), 39-49.
#' \doi{10.1177/014662169501900105}
#' @examples
#' betas <- lapply(seq(-3, 3, length.out = 12), function(b) b + c(-0.5, 0.5))
#' names(betas) <- sprintf("a%02d", 1:12)
#' X <- kits_sim(rnorm(800, 0, 1.5), betas, seed = 1)
#' item <- kits_item(X)
#' pr <- kits_precision(item, list(easy = sprintf("a%02d", 1:8), hard = sprintf("a%02d", 5:12)),
#'                      criterion = 0.6)
#' pr
#' p <- plot(pr)
#' if (interactive()) p
#' @export
kits_precision <- function(item, forms, theta = seq(-6, 6, by = 0.1), link = NULL, criterion = NULL) {
  call <- call_text(sys.call())
  parms <- param_frame(item)
  need_items(data.frame(visit_id = "p", item_id = parms$item_id, item_score = 0L), parms)
  abort_if(!is.numeric(theta) || length(theta) < 2 || any(!is.finite(theta)),
           "`theta` must be a grid of at least 2 finite numbers.")
  theta <- sort(unique(theta))
  links <- as_links(link)
  need_number(criterion, "criterion", "one positive number (the largest RMSE of the usable range)", above = 0, null = TRUE)
  sets <- as_forms(forms, parms)
  scale <- beta_scale_of(item)
  chain <- algo_line(links)
  if (is.null(criterion)) criterion <- line_stretch(chain, 1)   # 1 logit：6 到 8 題的題本 RMSE 最小約 0.6 到 0.7，達不到 0.5
  convert <- function(v) convert_by_links(v, links, from = if (is.na(scale)) NULL else scale)

  per_form <- imap(sets, \(its, f) form_precision(f, its, parms, theta, convert, chain, criterion))
  out <- per_form |>
    map("curve") |>
    list_rbind() |>
    kits_table("kits_precision")
  ranges <- per_form |>
    map("range") |>
    list_rbind()
  end_scale <- if (length(links) == 0) scale else links[[length(links)]]$to_scale
  attr(out, "scale") <- if (is.null(end_scale)) NA_character_ else end_scale
  attr(out, "criterion") <- criterion
  attr(out, "ranges") <- kits_table(select(ranges, -at_edge))
  attr(out, "diagnostics") <- precision_diagnostics(ranges, criterion)
  attr(out, "info") <- make_info(call, list(scale = scale, criterion = criterion, grid = range(theta)),
                                 item_data = source_of(item),
                                 links = if (length(links) > 0) {
                                   data.frame(intercept = map_dbl(links, "intercept"),
                                              slope = map_dbl(links, "slope"))
                                 })
  out
}

# 一個題本的精準度：每個 theta 的 E[WLE | theta]、偏誤、SD、RMSE（curve），和可用的範圍（range）
form_precision <- function(f, its, parms, theta, convert, chain, criterion) {
  form <- data.frame(pattern_id = "form", item_id = its)
  tab <- wle_table(filter(parms, item_id %in% its), design = form)
  P <- algo_lw(theta, betas_by_item(parms)[its])
  est <- tab$theta
  expected <- drop(P %*% est)
  sd <- sqrt(pmax(drop(P %*% est^2) - expected^2, 0))
  bias <- expected - theta
  rmse <- sqrt(bias^2 + sd^2)
  curve <- data.frame(form = f, theta = convert(theta), expected = convert(expected),
                      bias = line_stretch(chain, bias), sd = line_stretch(chain, sd),
                      rmse = line_stretch(chain, rmse))
  best <- which.min(curve$rmse)
  run <- usable_run(curve$rmse < criterion, best)
  range <- data.frame(form = f, n_items = length(its), max_raw = ncol(P) - 1,
                      from = curve$theta[run[1]], to = curve$theta[run[2]], min_rmse = curve$rmse[best], at = curve$theta[best],
                      at_edge = !is.na(run[1]) && (run[1] == 1 || run[2] == nrow(curve)))
  list(curve = curve, range = range)
}

precision_diagnostics <- function(ranges, criterion) {
  rows <- map(seq_len(nrow(ranges)), \(i) {
    r <- ranges[i, ]
    if (is.na(r$from)) {
      return(diag_row(paste("usable range,", r$form), NA, "note",
                      sprintf("RMSE never below %.3g on the grid (smallest %.3g): the form is too short for this criterion",
                              criterion, r$min_rmse)))
    }
    status <- if (r$at_edge) "note" else "ok"
    detail <- if (r$at_edge) "; the range reaches the end of the grid, so widen `theta`" else ""
    diag_row(paste("usable range,", r$form), r$to - r$from, status,
             sprintf("RMSE < %.3g from %.2f to %.2f%s", criterion, r$from, r$to, detail))
  })
  diag_table(rows)
}

#' @rdname kits_precision
#' @param x A `kits_precision`.
#' @param digits Decimals to print.
#' @param ... Not used.
#' @export
print.kits_precision <- function(x, digits = 2, ...) {
  need_digits(digits)
  cat(sprintf("<kits_precision> %d form(s), %d theta values each; scale: %s\n", length(unique(x$form)),
              nrow(x) / length(unique(x$form)), scale_label(scale_of(x))))
  cat(sprintf("Usable range (RMSE < %s):\n", format(signif(attr(x, "criterion"), 3))))
  print(format_table(attr(x, "ranges"), digits), row.names = FALSE, right = TRUE)
  cat("The curves: as.data.frame(x); plot(x).\n")
  checks_line(attr(x, "diagnostics"))
  invisible(x)
}

#' @rdname kits_precision
#' @export
plot.kits_precision <- function(x, ...) {
  cols <- c(PAL$blue, PAL$orange, "#009E73", "#CC79A7", "#E69F00", "#56B4E9", PAL$black, PAL$gray)
  forms <- unique(x$form)
  crit <- attr(x, "criterion")
  p <- plotly::plot_ly()
  for (i in seq_along(forms)) {
    d <- filter(as.data.frame(x), form == forms[i])
    hover <- sprintf("%s<br>theta %.2f: RMSE %.2f (bias %.2f, SD %.2f)", forms[i], d$theta, d$rmse, d$bias, d$sd)
    p <- plotly::add_lines(p, x = d$theta, y = d$rmse, name = forms[i], hovertext = hover, hoverinfo = "text",
                           line = list(color = cols[(i - 1) %% length(cols) + 1], width = LINE_W$main))
  }
  p <- plotly::add_lines(p, x = range(x$theta), y = rep(crit, 2), name = sprintf(LAB$criterion, crit), hoverinfo = "skip",
                         line = list(color = PAL$gray, width = LINE_W$thin, dash = "dash"))
  scale <- scale_of(x)
  x_title <- if (is.na(scale)) "theta" else if (scale == "D") "D" else sprintf("theta (%s)", scale)
  plotly::layout(p, xaxis = list(title = x_title), yaxis = list(title = LAB$rmse, rangemode = "tozero")) |>
    kits_style(LAB$title_precision)
}
