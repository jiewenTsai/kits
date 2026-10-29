# mean/sigma 的直線：algo_mean_sigma()（algo.R）；新值 = intercept + slope x 舊值。slope = 兩邊 SD 的比，
# intercept 讓兩邊的平均對齊（Kolen & Brennan 的 B、A）。不篩題的結果是連結的 intercept0、slope0。

# 共用斜率、各組截距的 mean/sigma（內部用，給 kits_linking() 的領域）。
# 每組的 from、to 先減掉自己組的平均；合併後的 SD 比是共用斜率；每組的截距 = 該組 to 的平均 - 斜率 x 該組
# from 的平均。
mean_sigma_shared_slope <- function(from, to, group) {
  d <- data.frame(from = from, to = to, group = group) |>
    mutate(dx = from - mean(from), dy = to - mean(to), .by = group)   # 減掉自己組的平均
  slope <- sqrt(sum(d$dy^2) / sum(d$dx^2))
  per_group <- d |>
    summarise(intercept = mean(to) - slope * mean(from), .by = group) |>
    arrange(sort_rank(group))
  list(intercept = stats::setNames(per_group$intercept, per_group$group), slope = slope)
}

# kits_item -> 它的門檻表；data.frame 照用。
# domains = TRUE 時也收有名字的 list（一個領域一個 kits_item 或門檻表），疊起來並加 domain 欄。
as_steps <- function(x, arg, domains = FALSE) {
  if (inherits(x, "kits")) x <- x$item
  if (inherits(x, "kits_item")) return(x$thresholds)
  if (domains && is.list(x) && !is.data.frame(x)) {
    abort_if(is.null(names(x)) || any(names(x) == "") || anyDuplicated(names(x)),
             sprintf("`%s` is a list: name each element by its domain, e.g. list(lan = item_lan, cog = item_cog).", arg))
    parts <- imap(x, \(item, d) {
      s <- as.data.frame(as_steps(item, arg))
      need_cols(s, c("item_id", "item_score", "thr"), arg)
      data.frame(item_id = s$item_id, item_score = s$item_score, thr = s$thr,
                 thr_se = if ("thr_se" %in% names(s)) s$thr_se else NA_real_, domain = d)
    })
    return(list_rbind(unname(parts)))
  }
  abort_if(!is.data.frame(x),
           sprintf("`%s` must be a kits_item() result or a table with item_id, item_score, thr.", arg))
  abort_if(rp_of(x) != 0.5,
           sprintf("`%s` holds RP%g thresholds (coef(item, rp = %g)). Links use the RP50 thresholds, ", arg, 100 * rp_of(x),
                   rp_of(x)), "as the GSED tau is an RP50 value: give the kits_item() or coef(item).")
  x
}

# 檢查門檻表：要有這些欄，(item_id, item_score) 不能重複
need_steps <- function(steps, arg) {
  need_cols(steps, c("item_id", "item_score", "thr"), arg)
  abort_if(!is.numeric(steps$thr) || any(!is.finite(steps$thr)),
           sprintf("`%s$thr` must be finite numbers (no NA or Inf).", arg), class = "kits_error_items")
  dup <- duplicated(steps[c("item_id", "item_score")])
  abort_if(any(dup),
           sprintf("`%s` has repeated item steps (for example %s:%s). Use the table of one cohort only.",
                   arg, steps$item_id[dup][1], steps$item_score[dup][1]), class = "kits_error_items")
}


# 種子：一個整數，或 NULL（不能重現）
need_seed <- function(seed) {
  need_number(seed, "seed", "one whole number (or NULL for a run that cannot be repeated)",
              lower = -.Machine$integer.max, upper = .Machine$integer.max, null = TRUE)
}

# 點的值：數字，沒有 NA 或無限大
need_link_values <- function(from, to) {
  abort_if(!is.numeric(from) || !is.numeric(to), "`from` and `to` must be numbers.")
  abort_if(anyNA(from) || anyNA(to) || any(!is.finite(c(from, to))),
           "`from` and `to` must not contain NA or infinite values.")
}

need_cut <- function(cut) {
  need_number(cut, "cut", "one positive number (in residual SDs)", above = 0, finite = FALSE)
}

need_n_boot <- function(n_boot) {
  need_number(n_boot, "n_boot", "0 (no linking error) or a whole number of resamples of at least 2", lower = 0,
              whole = TRUE, not = 1)
}

# 連結的兩端（量尺的標籤）：to 一個，from 一個或多個（領域）；兩端不能相同
link_ends <- function(from_scale, to_scale) {
  from_scale <- as.character(from_scale)
  to_scale <- as.character(to_scale)
  abort_if(length(to_scale) != 1 || length(from_scale) == 0,
           "`to_scale` must be one label, and `from_scale` one or more labels (or NA).")
  abort_if(!is.na(to_scale) && to_scale %in% from_scale,
           sprintf("Both ends of the link are on scale %s; a link joins two different scales.", to_scale),
           class = "kits_error_scale")
  list(from_scale = from_scale, to_scale = to_scale)
}

# 選項（像 match.arg，但錯誤有類別）：預設給整個向量時用第一個
match_option <- function(x, choices, arg) {
  if (identical(x, choices)) {
    return(choices[1])
  }
  abort_if(!is.character(x) || length(x) != 1 || !x %in% choices,
           sprintf("`%s` must be one of %s.", arg, paste0("\"", choices, "\"", collapse = ", ")))
  x
}

# 篩點的規則（內部用）：spread、unit，spread = "se" 時每個點的 SE（to 的 SE 沒給時是 0）
screen_rule <- function(spread = "sd", unit = "step", cut = 2, se_from = NULL, se_to = NULL, n = NULL) {
  if (spread != "se") {
    return(list(spread = spread, unit = unit))
  }
  abort_if(is.null(se_from),
           "`spread = \"se\"` needs the SE of each point: give `se_from` (and `se_to` when the target has SEs).")
  if (is.null(se_to)) se_to <- rep(0, length(se_from))
  if (!is.numeric(se_from) || !is.numeric(se_to) || length(se_from) != n || length(se_to) != n ||
      any(!is.finite(c(se_from, se_to))) || any(c(se_from, se_to) < 0)) {
    kits_abort("`se_from` and `se_to` must be finite numbers of at least 0, one per point.", class = "kits_error_input")
  }
  list(spread = spread, unit = unit, se_from = se_from, se_to = se_to)
}

# 篩點的尺度（每個點一個值）：|e| > cut x 尺度的點被標記。spread = "se" 時，from 的 SE 乘斜率換到目標量尺。
screen_scale <- function(resid, rule, slope0) {
  s <- switch(rule$spread,
    sd = stats::sd(resid),
    mad = stats::mad(resid),
    se = sqrt(rule$se_to^2 + slope0^2 * rule$se_from^2))
  abort_if(rule$spread == "se" && any(s <= 0),
           "`spread = \"se\"`: some points have SE 0 on both sides (fixed items?); their residuals cannot be standardized.")
  rep_len(s, length(resid))
}

# 篩點留下不到一半的點：規則太嚴。spread = "se" 在大樣本時最常見（門檻的 SE 很小）
half_kept_message <- function(spread, n_kept, n) {
  kept <- sprintf("the screening keeps %d of %d points, fewer than half", n_kept, n)
  if (identical(spread, "se")) {
    return(paste0(kept, ": spread = \"se\" divides by the SE of the thresholds, which is tiny in large samples ",
                  "(in KIT, thr_se is about 0.015 while the residuals are about 0.2), so it flags most points; ",
                  "use spread = \"sd\" (the default) or \"mad\""))
  }
  paste0(kept, ": the rule is too strict for these points; check them, or use a larger `cut`")
}

# 連結的診斷：留下的點、連結誤差、斜率和 1 比、共同題的比例、領域
link_diagnostics <- function(x) {
  rows <- c(screening_rows(x), list(link_error_row(x), slope_row(x), common_share_row(x)),
            domain_rows(x))
  diag_table(compact(rows))
}

# 篩點：留下的點、丟掉的比例、留下的點的相關；留下不到一半時警告
screening_rows <- function(x) {
  pts <- x$points
  kept <- !pts$dropped
  n_kept <- sum(kept)
  share <- 100 * mean(pts$dropped)
  level <- 1 + (share > 10) + (share > 25)   # ok、note（10% 以上）、warn（25% 以上）
  r <- stats::cor(pts$from[kept], pts$to[kept])
  list(
    diag_flag("points kept", n_kept, n_kept < 5, "warn",
              "fewer than 5 points: the line rests on very few points", sprintf("%d of %d", n_kept, nrow(pts))),
    diag_row("share dropped", share, c("ok", "note", "warn")[level],
             sprintf("%.1f%% of the points dropped by the screening%s", share,
                     c("", " (above 10%)", " (above 25%)")[level])),
    diag_flag("correlation of kept points", r, !is.finite(r) || r < 0.8, "warn",
              "below 0.8: the two scales do not agree well on these points", "from and to agree"),
    if (n_kept < nrow(pts) / 2) {
      diag_row("screening keeps fewer than half", n_kept / nrow(pts), "warn",
               half_kept_message(x$settings$spread, n_kept, nrow(pts)))
    })
}

# 連結誤差（n_boot > 0 時）
link_error_row <- function(x) {
  n_boot <- x$settings$n_boot
  if (is.na(x$intercept_error)) {
    return(if (isTRUE(n_boot > 0)) {
      diag_row("linking error", NA, "note", "could not be computed (all points belong to one item)")
    })
  }
  many_failed <- x$n_failed / n_boot > 0.1
  diag_row("linking error", x$intercept_error, if (many_failed) "note" else "ok",
           sprintf("SD of the intercept over %d resamples (slope: %.3g)%s", n_boot, x$slope_error,
                   if (many_failed) sprintf("; %d resamples failed", x$n_failed) else ""))
}

# 斜率和 1 比（兩個校準之間，有重抽時）：PCM 的連結只需要平移
slope_row <- function(x) {
  if (is.na(x$slope_error) || x$slope_error <= 0 || identical(x$settings$kind, "anchors")) {
    return(NULL)
  }
  z <- (x$slope - 1) / x$slope_error
  diag_flag("slope vs 1", z, abs(z) > 2, "note",
            sprintf(paste0("(slope - 1) / slope_error = %.2f: the slope %.3f is not 1, so the scale stretches ",
                           "or shrinks between the two calibrations (the items discriminate differently); ",
                           "worth reporting"), z, x$slope),
            "(slope - 1) / slope_error: the slope is within 2 linking errors of 1, as the PCM expects")
}

# 留下的共同題佔各校準題數的比例（Kolen & Brennan p. 288 的 20%）
common_share_row <- function(x) {
  if (is.null(x$items_total)) {
    return(NULL)
  }
  kept_items <- length(unique(x$points$item[!x$points$dropped]))
  share_items <- kept_items / x$items_total
  low <- any(share_items < 0.2)
  diag_row("common items / items", 100 * min(share_items), if (low) "note" else "ok",
           sprintf("%d kept common items: %.0f%% of the %d items of from, %.0f%% of the %d items of to%s",
                   kept_items, 100 * share_items[["from"]], x$items_total[["from"]], 100 * share_items[["to"]],
                   x$items_total[["to"]],
                   if (low) "; below 20% (Kolen & Brennan pp. 288, 336, a rule for tests of 40+ items)" else ""))
}

# 有領域時：留下的點在各領域的比例（p. 239）、各領域截距離共同的 a 多遠、沒有留下點的領域
domain_rows <- function(x) {
  dom <- x$domains
  if (is.null(dom)) {
    return(NULL)
  }
  top <- which.max(dom$n)
  unbalanced <- dom$n[top] / sum(dom$n) > 2 / 3
  gaps <- abs(dom$intercept - x$intercept)   # 沒有留下點的領域沒有截距（NA）
  dev <- if (all(is.na(gaps))) NA_real_ else max(gaps, na.rm = TRUE)
  big <- isTRUE(!is.na(x$intercept_error) && dev > 2 * x$intercept_error)
  dev_msg <- if (big) {
    "a domain intercept is more than 2 linking errors from the common intercept: one common line may not fit"
  } else if (is.na(dev)) {
    "no domain has kept points"
  } else if (is.na(x$intercept_error)) {
    "largest distance of a domain intercept from the common intercept (domains with kept points)"
  } else {
    "within 2 linking errors"
  }
  lost <- dom$domain[dom$n == 0]
  list(
    if (nrow(dom) > 1) {
      diag_row("domain balance of kept points", dom$n[top] / sum(dom$n), if (unbalanced) "note" else "ok",
               paste0(domain_balance_text(dom), if (unbalanced) {
                 sprintf(paste0("; %s has more than 2/3: the common items are not a miniature of the test ",
                                "(Kolen & Brennan p. 239)"), dom$domain[top])
               }))
    },
    diag_row("max |domain intercept - intercept|", dev, if (big) "note" else "ok", dev_msg),
    if (length(lost) > 0) diag_row("domains without kept anchors", length(lost), "warn", paste(lost, collapse = ", ")))
}


#' Print or plot a link
#'
#' `print()` shows the line, the linking error (when `n_boot > 0`), and the dropped points;
#' `coef()` gives the intercept and the slope.
#' `plot()` draws every point (`from` against `to`) with the final line, the line before
#' screening and the screening band (dashed); dropped points are marked.
#'
#' @param x The result of [kits_linking()] or [kits_link_line()].
#' @param digits Decimals to print (only the printout is rounded; the object keeps full precision).
#' @param ... Not used.
#' @return `print()`: `x` (invisibly). `plot()`: a plotly plot. `coef()`: `c(intercept = , slope =
#'   )`.
#'   `$`, `[[`: the element.
#' @examples
#' points <- data.frame(id = paste0("p", 1:6), item = paste0("p", 1:6),
#'                      from = c(-2, -1, 0, 1, 2, 3), to = c(-1, 0, 1, 2, 3.2, 6))
#' link <- kits_link_line(kits_link_screen(points))
#' link
#' coef(link)
#' p <- plot(link)
#' if (interactive()) p
#' @export
print.kits_link <- function(x, digits = 2, ...) {
  need_digits(digits)
  err <- function(e) if (length(e) != 1 || is.na(e)) "" else sprintf("  (linking error %s)", fmt(e, digits))
  cat("<kits_link> new value = intercept + slope x old value\n")
  if (!is.null(x$from_scale) && !(all(is.na(x$from_scale)) && all(is.na(x$to_scale)))) {
    cat(sprintf("  scales: %s -> %s\n", scale_label(x$from_scale), scale_label(x$to_scale)))
  }
  cat(sprintf("  intercept = %s%s\n", fmt(x$intercept, digits), err(x$intercept_error)))
  cat(sprintf("  slope     = %s%s\n", fmt(x$slope, digits), err(x$slope_error)))
  if (!is.null(link_points_text(x))) cat(sprintf("  points from: %s\n", link_points_text(x)))
  cat(sprintf("  points: %d%s; dropped: %d%s\n", x$n, if (is.null(x$n_items)) "" else sprintf(" (items: %d)", x$n_items),
              length(x$dropped), if (length(x$dropped) > 0) paste0(" (", paste(x$dropped, collapse = ", "), ")") else ""))
  cat(sprintf("  before screening: intercept %s, slope %s; screening: %s\n", fmt(x$intercept0, digits),
              fmt(x$slope0, digits), screen_rule_text(x$settings)))
  if (!is.null(x$domains) && nrow(x$domains) > 1) cat(sprintf("  kept points by domain: %s\n", domain_balance_text(x$domains)))
  if (!is.null(x$settings$no_anchors)) {
    cat(sprintf("  same line, no anchors: %s (KIT scale scores in D units, not GSED D)\n",
                paste(x$settings$no_anchors, collapse = ", ")))
  }
  if (!is.null(x$domains)) {
    cat("  domain check: one slope shared by the domains, fitted with domain intercepts (kept points);\n",
        "  it is not the slope of the line above:\n", sep = "")
    print(format_table(x$domains, digits), row.names = FALSE, right = TRUE)
  }
  cat("Next: kits_convert(scores, x), plot(x), summary(x, theta = ), confint(x).\n")
  checks_line(x$diagnostics)
  invisible(x)
}

#' @rdname print.kits_link
#' @param object A link.
#' @export
coef.kits_link <- function(object, ...) {
  c(intercept = object$intercept, slope = object$slope)
}

#' @rdname print.kits_link
#' @export
plot.kits_link <- function(x, ...) {
  pts <- as.data.frame(x$points)
  pts$hover <- sprintf("%s<br>from %.2f, to %.2f<br>residual %.2f", pts$id, pts$from, pts$to, pts$resid_final)
  line_x <- range(pts$from)
  first <- x$intercept0 + x$slope0 * line_x   # 篩點之前的線
  spread <- if (is.null(x$settings$spread)) "sd" else x$settings$spread
  p <- plotly::plot_ly() |>
    plotly::add_lines(x = line_x, y = line_apply(x, line_x),
                      line = list(color = PAL$black, width = LINE_W$main),
                      name = sprintf(LAB$final_line, x$intercept, x$slope), hoverinfo = "skip") |>
    plotly::add_lines(x = line_x, y = first,
                      line = list(color = PAL$gray, width = LINE_W$thin, dash = "dot"),
                      name = sprintf(LAB$first_line, x$intercept0, x$slope0), hoverinfo = "skip")
  if (spread != "se") {   # 篩選的邊界：第一條線 +/- cut x 尺度（"se" 每個點的尺度不同，不畫）
    band <- x$settings$cut * if (spread == "mad") stats::mad(pts$resid) else stats::sd(pts$resid)
    p <- plotly::add_lines(p, x = c(line_x, NA, line_x), y = c(first + band, NA, first - band),
                           line = list(color = PAL$orange, width = LINE_W$thin, dash = "dash"),
                           name = sprintf(if (spread == "mad") LAB$band_mad else LAB$band, x$settings$cut), hoverinfo = "skip")
  }
  for (drop in c(FALSE, TRUE)) {
    d <- pts[pts$dropped == drop, ]
    if (nrow(d) == 0) next
    p <- plotly::add_markers(p, x = d$from, y = d$to, hovertext = d$hover, hoverinfo = "text",
                             name = if (drop) LAB$dropped else LAB$kept,
                             marker = list(color = if (drop) PAL$orange else PAL$blue, size = MARKER_SIZE + 1,
                                           symbol = if (drop) "x" else "circle"))
  }
  plotly::layout(p, xaxis = list(title = LAB$link_from), yaxis = list(title = LAB$link_to)) |>
    kits_style(LAB$title_link)
}

# 整題重抽 n_boot 次，每次重跑篩點和直線：algo_link_error()（algo.R）


#' Link two scales: between two calibrations, or to D through anchor items
#'
#' Fits a straight mean/sigma line, new value = intercept + slope x old value, with one round of
#' screening, and (with `n_boot > 0`) a linking error that resamples whole items.
#' Between two calibrations each common step is one point; to an external scale each anchor item
#' is one point (its pass step on this scale against its `tau`). Several domains can share one
#' line. The four steps are also functions of their own ([kits_link_steps]).
#' Details: `kits_algorithm("link")` and tutorial_algorithms/03_linking.html.
#'
#' @references
#' Kolen, M. J., & Brennan, R. L. (2014). *Test equating, scaling, and linking* (3rd ed.). Springer.
#' \doi{10.1007/978-1-4939-0317-7}
#'
#' @param from The calibration to be converted: a [kits_item()] result (for example the M36 cohort), a [kits()],
#'   its thresholds table (optionally with a `domain` column), or a named list of them, one per domain.
#' @param to The target: a calibration ([kits_item()], [kits()] or a thresholds table), or a
#'   data.frame of anchor items with `item_id`, `tau` and optionally `step` (the pass step).
#' @param step Which step of each item gives the points.
#'   Between two calibrations: `NULL` (every common step) or one whole number.
#'   To anchor items: the pass step, "pass" is a score of at least `step`; `NULL` means 1, the milestone step.
#' @param to_scale The label of the target scale.
#'   Between two calibrations it is the scale of `to` (a different label stops with `kits_error_scale`);
#'   to anchor items it is `"D"` unless given.
#' @param cut The cutoff for dropping, in units of `spread` (a positive number). At least 3 points
#'   must be kept.
#' @param spread The spread of the screening: `"sd"` (default, the SD of the residuals), `"mad"`
#'   (1.4826 x their MAD) or `"se"` (the SE of each point, from `thr_se`; for small samples only).
#' @param unit What is dropped: `"step"` (default, the flagged points) or `"item"` (every point of
#'   an item with a flagged point).
#' @param n_boot Number of item resamples for the bootstrap linking error; 0 (default) skips it.
#'   200 is enough to see its size; use more (e.g. 2000) for reporting.
#' @param seed Random seed for the resampling, or `NULL` for a run that cannot be repeated.
#'   Your current random state is left unchanged.
#' @return A `kits_link`: a list that prints briefly ([print.kits_link()]). Its elements:
#'   * `intercept`, `slope`: the line after screening (the B and A of Kolen & Brennan, 2014).
#'   * `intercept_error`, `slope_error`, `n_failed`: the linking error (NA when `n_boot = 0`).
#'   * `draws`, `vcov` (only when `n_boot > 0`): the resampled lines and their 2 x 2 covariance.
#'   * `intercept0`, `slope0`: the line before screening.
#'   * `n`: the number of points.
#'   * `dropped`: the dropped points.
#'   * `points`: one row per point: `id`, `item`, `from`, `to`, `resid`, `resid_final`, `z`,
#'   `flagged`, `dropped` (and `se_from`, `se_to`, `step`, `domain` when they apply).
#'   * `settings`: `cut`, `n_boot`, `seed`, `spread`, `unit`, `kind`, `step` (and `key`).
#'   * `from_scale`, `to_scale`: the scales the line converts from and to (NA when not given).
#'   * `diagnostics`: the checks of the link (see [kits_check()]); `info`: see [kits_info()].
#'   * `n_items` and, between two calibrations, `items_total`.
#'   * `domains` (with domains): `domain`, `n`, `intercept`, and the shared `slope`.
#' @examples
#' # Each cohort answers part of the items; a3, a4, a5 are common
#' betas <- list(a1 = -2, a2 = -1, a3 = c(-0.5, 0.5), a4 = 0, a5 = c(0.5, 1.5), a6 = 2, a7 = 3)
#' young <- kits_sim(rnorm(500, 0), betas[1:5], seed = 1)
#' old <- kits_sim(rnorm(500, 2), betas[3:7], prefix = "q", seed = 2)
#' item_old <- kits_item(old, scale = "old cohort")
#' item_young <- kits_item(young, scale = "young cohort")
#' link <- kits_linking(item_old, item_young, n_boot = 200)
#' link
#'
#' # The older cohort's theta on the younger cohort's scale
#' s_old <- kits_subj(old, item_old)
#' head(kits_convert(s_old, link))
#'
#' # To D through anchor items, two domains sharing one line
#' set.seed(1)
#' steps <- data.frame(item_id = paste0("a", 1:12), item_score = 1, thr = rnorm(12),
#'                     domain = rep(c("lan", "cog"), each = 6))
#' anchors <- data.frame(item_id = steps$item_id, tau = 64 + 7 * steps$thr + rnorm(12, 0, 0.5))
#' d_link <- kits_linking(steps, anchors, n_boot = 200)
#' d_link
#' kits_convert(c(-1, 0, 1), d_link, from = "my scale")
#' @export
kits_linking <- function(from, to, step = NULL, cut = 2, n_boot = 0, seed = 1, to_scale = NULL,
                         spread = c("sd", "mad", "se"), unit = c("step", "item")) {
  call <- call_text(sys.call())
  spread <- match_option(spread, c("sd", "mad", "se"), "spread")
  unit <- match_option(unit, c("step", "item"), "unit")
  need_seed(seed)
  need_cut(cut)
  need_n_boot(n_boot)
  link <- points_step(from, to, step, to_scale) |>   # 1 每個共同步驟（或錨題）一個點
    screen_step(cut, spread, unit) |>                # 2 篩一次
    line_step() |>                                   # 3 留下的點的直線
    error_step(n_boot, seed)                         # 4 連結誤差（n_boot > 0 時）
  link$info <- make_info(call, link$settings, link$info$data)
  link
}

# to 是 data.frame 但不是錨題表也不是門檻表時，說清楚兩種表各要哪些欄
need_link_target <- function(to) {
  if (!is.data.frame(to) || inherits(to, "kits_item") || any(c("thr", "tau") %in% names(to))) {
    return(invisible())
  }
  kits_abort("`to` has neither a `tau` nor a `thr` column (it has: ", paste(names(to), collapse = ", "), "). ",
             "An anchor table (a link to D) needs item_id and tau (the GSED difficulty); ",
             "a thresholds table (a link between two calibrations) needs item_id, item_score and thr, as coef(item). ",
             "Rename the columns, e.g. dplyr::rename(anchors, item_id = item, tau = difficulty).", class = "kits_error_input")
}

# to 是錨題表（item_id、tau，沒有門檻 thr）時，是接到外部量尺的連結
is_anchors <- function(to) {
  is.data.frame(to) && !inherits(to, "kits_item") && "tau" %in% names(to) && !"thr" %in% names(to)
}

# 兩個校準之間的連結：目標量尺是 to 記錄的量尺；給了不一樣的 to_scale 時停下來
target_scale <- function(to, to_scale) {
  recorded <- scale_of(to)
  if (is.null(to_scale)) {
    return(recorded)
  }
  abort_if(!is.na(recorded) && recorded != to_scale,
           sprintf("`to` is on scale %s, not %s. Leave `to_scale` out, or relabel the calibration.", recorded, to_scale),
           class = "kits_error_scale")
  to_scale
}

# 兩個校準共同的題目步驟：一個步驟一個點（step 給了時只用那個步驟）
common_points <- function(steps, to, step) {
  to <- as_steps(to, "to")
  need_steps(to, "to")
  common <- inner_join(select(as.data.frame(steps), any_of(c("item_id", "item_score", "thr", "thr_se", "domain"))),
                       select(as.data.frame(to), any_of(c("item_id", "item_score", "thr", "thr_se"))),
                       by = c("item_id", "item_score"), suffix = c("_from", "_to")) |>
    arrange(sort_rank(item_id), item_score)
  if (!is.null(step)) common <- filter(common, item_score == step)
  abort_if(nrow(common) < 3,
           "The two tables have ", nrow(common), " common item step(s)", if (!is.null(step)) sprintf(" at step %d", step),
           "; at least 3 are needed.", class = "kits_error_items")
  out <- transmute(common, id = paste0(item_id, ":", item_score), item = item_id, from = thr_from, to = thr_to,
                   across(any_of("domain")))
  if (all(c("thr_se_from", "thr_se_to") %in% names(common))) {
    out <- mutate(out, se_from = common$thr_se_from, se_to = common$thr_se_to)
  }
  out
}

# 錨題：每題通過步驟（至少 step 分）的門檻是 from，外部量尺上的難度 tau 是 to
anchor_points <- function(steps, anchors, step) {
  need_cols(anchors, c("item_id", "tau"), "to")
  if (anyDuplicated(anchors$item_id)) {
    dup <- anchors$item_id[duplicated(anchors$item_id)][1]
    keys <- if ("key" %in% names(anchors)) unique(as.character(anchors$key)) else character()
    abort_if(length(keys) > 1,
             "The anchors (`to`) have repeated item_id (for example ", dup, ") because they hold more than one key (",
             paste(keys, collapse = ", "), "). The keys put D on different metrics: choose one key, ",
             "e.g. kits_anchors(\"gsed2406\").")
    kits_abort("The anchors (`to`) have repeated item_id: ", dup, ". Give each anchor item once ",
               "(with anchors of several keys, choose one key).", class = "kits_error_input")
  }
  abort_if(!is.numeric(anchors$tau) || any(!is.finite(anchors$tau)), "`to$tau` must be finite numbers.")
  pass_step <- step   # 錨題表的 step 欄是 NA 時用這個
  anchors <- anchors |>
    as.data.frame() |>
    mutate(step = if ("step" %in% names(anchors)) ifelse(is.na(step), pass_step, step) else pass_step)
  pass <- select(as.data.frame(steps), any_of(c("item_id", "item_score", "thr", "thr_se", "domain"))) |>
    inner_join(select(anchors, item_id, step, tau), by = "item_id") |>
    arrange(sort_rank(item_id))
  missing_step <- pass |>
    summarise(has_step = any(item_score == step), .by = item_id) |>
    filter(!has_step) |>
    pull(item_id)
  abort_if(length(missing_step) > 0,
           "Anchor item(s) without the pass step in `from`: ", paste(missing_step, collapse = ", "), ".",
           class = "kits_error_items")
  an <- filter(pass, item_score == step)
  abort_if(nrow(an) < 3,
           nrow(an), " anchor item(s) found in `from`; at least 3 are needed.", class = "kits_error_items")
  out <- transmute(an, id = item_id, item = item_id, from = thr, to = tau, step, across(any_of("domain")))
  if ("thr_se" %in% names(an)) out$se_from <- an$thr_se   # 錨題（GSED 的 tau）沒有 SE：只用 from 的 SE
  out
}

# 有領域時：點記下領域；用同樣留下的點算共用斜率和各領域截距
add_domains <- function(link, domain) {
  link$points$domain <- domain
  keep <- !link$points$dropped
  shared <- mean_sigma_shared_slope(link$points$from[keep], link$points$to[keep], domain[keep])
  link$domains <- as.data.frame(link$points) |>
    summarise(n = sum(!dropped), .by = domain) |>
    arrange(sort_rank(domain)) |>
    mutate(intercept = unname(shared$intercept[domain]), slope = shared$slope) |>
    kits_table()
  lost <- link$domains$domain[link$domains$n == 0]
  if (length(lost) > 0) {
    kits_warn("All points of domain(s) ", paste(lost, collapse = ", "), " were dropped by the screening.",
              class = "kits_warning_scale")
  }
  link$diagnostics <- link_diagnostics(link)
  link
}

# 留下的點在各領域的數目和比例，例如 "motor 22 (88%), lan 2 (8%), cog 1 (4%)"（多的在前）
domain_balance_text <- function(domains) {
  d <- as.data.frame(domains) |>
    mutate(share = 100 * n / sum(n)) |>
    arrange(desc(n), sort_rank(domain))
  paste(sprintf("%s %d (%.0f%%)", d$domain, as.integer(d$n), d$share), collapse = ", ")
}

# 篩點的規則寫成一行字（給 print、kits_algorithm()）
screen_rule_text <- function(st) {
  spread <- if (is.null(st$spread)) "sd" else st$spread
  what <- switch(spread, sd = "|e| > %s SD(e)", mad = "|e| > %s x 1.4826 MAD(e)",
                 se = "|e| / sqrt(se_to^2 + slope0^2 se_from^2) > %s")
  paste0(sprintf(what, format(st$cut)), if (isTRUE(st$unit == "item")) ", whole items")
}

# 連結用了哪些點（給 print）；自己準備的點（data.frame）是 NULL
link_points_text <- function(x) {
  st <- x$settings
  if (is.null(st$kind)) {
    return(NULL)
  }
  if (st$kind == "calibration") {
    return(if (identical(st$step, "all")) "every common item step" else sprintf("step %s of each common item", st$step))
  }
  own <- !is.null(x$points$step) && any(x$points$step != st$step)
  paste0(sprintf("the pass step of each anchor item (score >= %s%s)", st$step, if (own) "; some anchors set their own" else ""),
         if (!is.null(st$key)) sprintf(", tau of dscore key %s", st$key) else "")
}
