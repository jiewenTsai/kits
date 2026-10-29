#' @include style.R
NULL

# plot() 的方法（plotly，和其他的圖同一個樣式）：校準的類別機率、性別 DIF、題目適配、通過月齡、偏移、篩波次。
# 換一種圖用 type =（第一個是預設），選題用 items =（題號）。

# 要畫的題：items 是題號；沒給時全部
pick_items <- function(all, items) {
  if (is.null(items)) {
    return(all)
  }
  abort_if(!is.character(items) || length(items) == 0, "`items` must be item ids, e.g. items = c(\"a1\", \"a2\").")
  miss <- setdiff(items, all)
  abort_if(length(miss) > 0, "Item(s) not in the result: ", paste(utils::head(miss, 10), collapse = ", "), ".")
  intersect(all, items)
}

# 類別的顏色：0 分灰色，再來藍、橘、綠……
CATEGORY_COLS <- c(PAL$gray, PAL$blue, PAL$orange, "#009E73", "#CC79A7", "#E69F00", "#56B4E9",
                   PAL$black)

#' @rdname kits_item
#' @param type What `plot()` draws: `"category"` (the default; the probability of each score
#'   against theta, one panel per item), `"expected"` (the expected score of each item), `"wright"`
#'   (a Wright map; needs `scores`, on the scale of the thresholds), `"fit"` (`x$itemfit`) or
#'   `"observed"` (needs `data`: the share of children who pass each item against their theta, one
#'   line per group, in bins of 0.5 logits with at least 10 children; the ticks at the top are the
#'   thetas of those who pass, at the bottom of those who fail; grey: the model; bottom right:
#'   infit and outfit of [kits_itemfit()] on `data`).
#' @param scores For `type = "wright"`: the scores of the children ([kits_subj()], or a data.frame
#'   with `theta`). For `type = "observed"`: optional, a data.frame with `visit_id` and `theta` on
#'   the scale of the betas; by default `kits_subj(data, x)`.
#' @param data For `type = "observed"`: the responses (as in [kits_item()]).
#' @param group For `type = "observed"`: the groups, one line each, as in [kits_dif()] (for
#'   example `"wave"`, `"sex"` or one value per row of `data`); by default one line.
#' @param step For `type = "observed"`: pass is `item_score >= step`, as in [kits_item_table()]
#'   (one number, or numbers named by item).
#' @param persons For `type = "observed"`: optional data.frame with `visit_id` and the column named
#'   by `group`, as in [kits_dif()].
#' @param items Which items to draw (item ids); by default all of them, or the first 9 when there
#'   are more (with a message).
#' @param theta The grid of theta (on the calibration scale); by default the range of the thresholds
#'   plus 3 logits on each side.
#' @examples
#'
#' # The probability of each score of each item
#' p <- plot(item)
#' p <- plot(item, items = c("a1", "a2"))
#' # The expected score of each item
#' p <- plot(item, type = "expected")
#' # A Wright map: the children beside the thresholds
#' p <- plot(item, type = "wright", scores = kits_subj(X, item))
#' # Item fit: infit and outfit of each item (x$itemfit)
#' p <- plot(item, type = "fit")
#' # Observed share passing against theta, by a group, beside the model (grey)
#' p <- plot(item, type = "observed", data = X, group = rep(c("a", "b"), length.out = nrow(X)))
#' if (interactive()) p
#' @export
plot.kits_item <- function(x, type = c("category", "expected", "wright", "fit", "observed"),
                           items = NULL, theta = NULL, scores = NULL, data = NULL, group = NULL,
                           step = 1, persons = NULL, ...) {
  type <- match_option(type, c("category", "expected", "wright", "fit", "observed"), "type")
  if (type == "observed") {
    return(item_observed_plot(x, items, data, group, step, persons, scores))
  }
  if (type == "wright") {
    return(item_wright_plot(x, items, scores))
  }
  if (type == "fit") {
    return(item_fit_plot(x, items))
  }
  parms <- param_frame(x)
  all <- sort(unique(parms$item_id))
  shown <- pick_items(all, items)
  if (is.null(items) && length(shown) > 9) {
    kits_inform("Drawing the first 9 of ", length(shown), " items; choose others with `items =`.",
                class = "kits_message_data")
    shown <- shown[1:9]
  }
  if (is.null(theta)) theta <- seq(min(parms$beta) - 3, max(parms$beta) + 3, by = 0.05)
  abort_if(!is.numeric(theta) || length(theta) < 2 || any(!is.finite(theta)),
           "`theta` must be a grid of at least 2 finite numbers.")
  beta_of <- betas_by_item(parms)
  x_title <- sprintf("theta (%s)", scale_label(beta_scale_of(x)))
  if (type == "expected") {
    return(item_expected_plot(shown, beta_of, theta, x_title))
  }
  panels <- map(seq_along(shown),
                \(i) item_curves_panel(shown[i], beta_of[[shown[i]]], theta, i == 1))
  rows <- ceiling(length(panels) / 3)
  plotly::subplot(panels, nrows = rows, shareX = TRUE, shareY = TRUE, margin = 0.04) |>
    plotly::layout(xaxis = list(title = x_title), yaxis = list(title = LAB$prob)) |>
    kits_style(LAB$title_item_curves)
}

# Wright map：孩子的分數和門檻要在同一條量尺上（兩邊都有記錄時檢查）
item_wright_plot <- function(x, items, scores) {
  abort_if(is.null(scores) || !is.data.frame(scores) || !"theta" %in% names(scores),
           "A Wright map needs the children: plot(item, type = \"wright\", scores = scores).")
  thr_scale <- scale_of(x)
  need_same_scale(thr_scale, scale_of(scores),
                  sprintf("The scores are on scale %s, the thresholds on %s; ", scale_of(scores), thr_scale),
                  "convert one of them first.")
  thresholds <- as.data.frame(x$thresholds)
  shown <- pick_items(unique(thresholds$item_id), items)
  wright_plot(filter(thresholds, item_id %in% shown), scores$theta)
}

# 題目適配的圖：x$itemfit 畫成 plot(kits_itemfit()) 的圖；算不出來（NULL）或舊的物件時停下來
item_fit_plot <- function(x, items) {
  abort_if(is.null(x$itemfit), "This calibration has no item fit (x$itemfit is NULL; see summary(x)). ",
           "Compute it on data with kits_itemfit(data, x), then plot() it.")
  plot.kits_itemfit(x$itemfit, items = items)
}

# 一題的類別機率曲線：每個分數一條線
item_curves_panel <- function(item, beta, theta, legend) {
  category_lines(item, theta, pcm_probs(theta, beta), legend) |>
    plotly::layout(annotations = panel_title(item), yaxis = list(range = c(0, 1)))
}

# 一題的類別機率曲線（plot.kits_item()、plot.kits_pcm_prob()）：P 的第 k 欄是 k - 1 分，一欄一條線
# （CATEGORY_COLS：灰、藍、橘……）。legend 是多格圖的這一格要不要列圖例；NULL：只有一格，不分組
category_lines <- function(item, theta, P, legend = NULL) {
  p <- plotly::plot_ly()
  for (k in seq_len(ncol(P))) {
    score <- k - 1
    hover <- sprintf("%s, score %d<br>theta %.2f: P = %.3f", item, score, theta, P[, k])
    line <- list(color = CATEGORY_COLS[(k - 1) %% length(CATEGORY_COLS) + 1], width = LINE_W$main)
    group <- if (!is.null(legend)) paste("score", score)
    p <- plotly::add_lines(p, x = theta, y = P[, k], name = sprintf("Score %d", score),
                           legendgroup = group, showlegend = legend, line = line,
                           hovertext = hover, hoverinfo = "text")
  }
  p
}

# 每題的期望分數：一題一條線，在同一格
item_expected_plot <- function(shown, beta_of, theta, x_title) {
  p <- plotly::plot_ly()
  for (i in seq_along(shown)) {
    expected <- algo_pcm(theta, beta_of[[shown[i]]])$E
    hover <- sprintf("%s<br>theta %.2f: expected score %.2f", shown[i], theta, expected)
    line <- list(color = COL_MANY[(i - 1) %% length(COL_MANY) + 1], width = LINE_W$main)
    p <- plotly::add_lines(p, x = theta, y = expected, name = shown[i], line = line,
                           hovertext = hover, hoverinfo = "text")
  }
  plotly::layout(p, xaxis = list(title = x_title), yaxis = list(title = LAB$expected_item)) |>
    kits_style(LAB$title_item_expected)
}

# 實際的通過率對 theta（D-score 書的等化組經驗曲線）：一題一格、一組一條線，灰線是模型的 P(X >= step)。
# theta 分成 0.5 logit 一段，一段至少 10 人才畫；上方的短線是通過的人的 theta，下方是沒通過的。
OBS_BIN <- 0.5
OBS_MIN_N <- 10

item_observed_plot <- function(x, items, data, group, step, persons, scores) {
  abort_if(is.null(data), "type = \"observed\" needs the responses: ",
           "plot(item, type = \"observed\", data = X).")
  need_pass_step(step)
  d <- observed_data(x, data, group, persons, scores)
  beta_of <- betas_by_item(param_frame(x))
  shown <- pick_items(intersect(names(beta_of), unique(d$item_id)), items)
  if (is.null(items) && length(shown) > 9) {
    kits_inform("Drawing the first 9 of ", length(shown), " items; choose others with `items =`.",
                class = "kits_message_data")
    shown <- shown[1:9]
  }
  unknown <- setdiff(names(step), names(beta_of))
  abort_if(length(unknown) > 0, "`step` names item(s) not in the parameters: ",
           some_ids(unknown, 10), ".")
  too_high <- shown[step_of(shown, step) > lengths(beta_of[shown])]
  abort_if(length(too_high) > 0, "`step` is above the highest step of item(s) ",
           some_ids(too_high, 10), ".")
  d <- filter(d, item_id %in% shown) |>
    mutate(pass = item_score >= step_of(item_id, step))
  # 適配算不出來時（例如題太少）不寫適配
  fit <- tryCatch(as.data.frame(kits_itemfit(data, x)),
                  error = \(e) data.frame(item_id = character()))
  theta <- seq(floor(min(d$theta)), ceiling(max(d$theta)), by = 0.05)
  levels <- attr(d, "levels")
  # 每組的圖例放在這組第一次出現的格子
  has <- \(it, g) any(d$item_id == it & d$group == g)
  first <- map_int(levels, \(g) match(TRUE, map_lgl(shown, \(it) has(it, g))))
  cols <- group_colours(length(levels))
  panels <- map(seq_along(shown), \(i) {
    observed_panel(filter(d, item_id == shown[i]), shown[i], beta_of[[shown[i]]], theta, levels,
                   cols, step_of(shown[i], step), fit[fit$item_id == shown[i], , drop = FALSE],
                   first == i, i == 1)
  })
  plotly::subplot(panels, nrows = ceiling(length(panels) / 3), shareX = TRUE, shareY = TRUE,
                  margin = 0.04) |>
    plotly::layout(xaxis = list(title = sprintf("theta (%s)", scale_label(beta_scale_of(x)))),
                   yaxis = list(title = LAB$pass_pct),
                   legend = list(title = list(text = attr(d, "label")))) |>
    kits_style(sprintf(LAB$title_item_observed, attr(d, "label")))
}

# 長格式的作答，每列接上組和 theta（沒有組或 theta 的列不用）
observed_data <- function(x, data, group, persons, scores) {
  resp <- resp_long(data, gaps = FALSE)
  if (is.null(scores)) {
    scores <- subj_theta(resp, param_frame(x))   # 只要 theta：kits_subj() 的 WLE，不算個人適配
  } else {
    abort_if(!is.data.frame(scores) || !all(c("visit_id", "theta") %in% names(scores)),
             "`scores` must be a data.frame with `visit_id` and `theta`, ",
             "e.g. kits_subj(data, item).")
    need_same_scale(beta_scale_of(x), scale_of(scores),
                    sprintf("The scores are on scale %s, the betas on %s; ", scale_of(scores),
                            beta_scale_of(x)), "use kits_subj(data, item).")
  }
  sp <- if (is.null(group)) {
    all <- LAB$all_children
    list(group = rep(all, nrow(resp)), levels = all, label = "all children")
  } else {
    lrt_split(group, data, resp, persons)
  }
  d <- as.data.frame(resp)[c("visit_id", "item_id", "item_score")]
  d$group <- sp$group
  d$theta <- scores$theta[match(d$visit_id, scores$visit_id)]
  lrt_n_na(d)   # 沒有組的人次說一聲
  d <- filter(d, !is.na(group), !is.na(theta), !is.na(item_score))
  abort_if(nrow(d) == 0, "No responses with both a group and a theta.")
  structure(d, levels = sp$levels[sp$levels %in% d$group], label = sp$label)
}

# 一題的一格：每組的通過率（點和線）、上下的短線、灰色的模型曲線、右下角的適配
observed_panel <- function(d, item, beta, theta, levels, cols, step, fit, legend, legend_model) {
  P <- pcm_probs(theta, beta)
  model <- 100 * rowSums(P[, -seq_len(step), drop = FALSE])
  p <- plotly::plot_ly() |>
    plotly::add_lines(x = theta, y = model, name = LAB$model, legendgroup = "model",
                      showlegend = legend_model,
                      line = list(color = PAL$gray, width = LINE_W$main, dash = "dash"),
                      hoverinfo = "skip")
  for (k in which(levels %in% d$group)) {   # 這題沒有的組不畫
    p <- observed_group(p, d[d$group == levels[k], ], item, levels[k], cols[k], legend[k])
  }
  note <- if (nrow(fit) == 1) sprintf("infit %.2f<br>outfit %.2f", fit$infit, fit$outfit)
  title <- if (step > 1) sprintf("%s (score >= %d)", item, step) else item
  ann <- c(panel_title(title), if (!is.null(note)) {
    list(list(x = 0.97, y = 0.1, xref = "paper", yref = "paper", text = note, showarrow = FALSE,
              xanchor = "right", yanchor = "bottom", align = "right", font = list(size = 10)))
  })
  y_axis <- list(range = c(-9, 109), tickvals = seq(0, 100, 20))
  plotly::layout(p, annotations = ann, yaxis = y_axis)
}

# 組的顏色：6 組以內照 COL_MANY；更多組（例如波次）從藍經灰到橘，照組的順序
group_colours <- function(n) {
  if (n <= length(COL_MANY)) COL_MANY[seq_len(n)] else
    grDevices::colorRampPalette(c(PAL$blue, PAL$gray, PAL$orange))(n)
}

# 一組：一段至少 OBS_MIN_N 人的通過率；上下的短線（theta 取到 0.05，重複的只畫一次）
observed_group <- function(p, g, item, name, colour, legend) {
  bins <- g |>
    mutate(bin = (floor(theta / OBS_BIN) + 0.5) * OBS_BIN) |>
    summarise(n = n(), n_pass = sum(pass), theta = mean(theta), .by = bin) |>
    filter(n >= OBS_MIN_N) |>
    arrange(bin)
  hover <- sprintf("%s, %s<br>theta %.2f: %d of %d pass (%.0f%%)", item, name, bins$theta,
                   bins$n_pass, bins$n, 100 * bins$n_pass / bins$n)
  ticks <- distinct(g, pass, theta = round(theta / 0.05) * 0.05)
  tick_marker <- list(symbol = "line-ns-open", size = 8, color = colour, opacity = 0.6,
                      line = list(width = LINE_W$thin))
  p |>
    plotly::add_trace(x = bins$theta, y = 100 * bins$n_pass / bins$n, type = "scatter",
                      mode = "lines+markers", name = name, legendgroup = name, showlegend = legend,
                      line = list(color = colour, width = LINE_W$main),
                      marker = list(color = colour, size = MARKER_SIZE - 2),
                      hovertext = hover, hoverinfo = "text") |>
    plotly::add_markers(x = ticks$theta, y = ifelse(ticks$pass, 104, -4), marker = tick_marker,
                        legendgroup = name, showlegend = FALSE, hoverinfo = "skip")
}


#' @rdname kits_dif
#' @param type What `plot()` draws: `"items"` (the default with two groups; each item's `dif`
#'   with its 95% CI, beyond +/- `cut` in colour) or `"steps"` (the graphical model check, the
#'   default with more groups: with two groups the beta of each step in the second group against
#'   the first, orange beyond Wald |z| 1.96; with more groups one panel per group against the joint
#'   calibration).
#' @param items Which items to draw (item ids); by default all of them.
#' @param cut For `type = "items"`: the size of DIF highlighted (logits; the `cut` of
#'   [kits_dif_combine()]).
#' @examples
#'
#' # Each item with its 95% CI; beyond +/- cut in colour; then each step in both groups
#' p <- plot(dif)
#' p <- plot(dif, type = "steps")
#' if (interactive()) p
#' @export
plot.kits_dif <- function(x, type = c("items", "steps"), items = NULL, cut = 0.3, ...) {
  if (identical(type, c("items", "steps")) && is.null(attr(x, "groups"))) type <- "steps"
  type <- match_option(type, c("items", "steps"), "type")
  if (type == "steps") {
    return(dif_steps_plot(x, items))
  }
  abort_if(is.null(attr(x, "groups")), "type = \"items\" needs the item DIF of two groups; use type = \"steps\".")
  need_dif_cut(cut)
  shown <- pick_items(x$item_id, items)
  se <- attr(x, "se")
  ref <- attr(x, "groups")$ref
  focal <- attr(x, "groups")$focal
  labels <- dif_labels(focal, ref)
  d <- as.data.frame(x) |>
    filter(item_id %in% shown) |>
    mutate(se = if (is.null(se)) NA_real_ else unname(se[item_id]),
           side = names(labels)[match(dif_direction(dif, cut, focal, ref), labels)],
           side = coalesce(side, "none")) |>
    arrange(dif)
  sides <- list(up = list(col = PAL$orange, name = sprintf("Favors %s (> %g)", focal, cut)),
                down = list(col = PAL$blue, name = sprintf("Favors %s (< -%g)", ref, cut)),
                none = list(col = PAL$gray, name = "Within the cut"))
  p <- plotly::plot_ly()
  for (kind in names(sides)) {
    s <- filter(d, side == kind)
    if (nrow(s) == 0) next
    hover <- sprintf("%s: dif %.2f (95%% CI %.2f to %.2f)", s$item_id, s$dif, s$dif - 1.96 * s$se,
                     s$dif + 1.96 * s$se)
    bars <- list(type = "data", array = 1.96 * s$se, color = sides[[kind]]$col)
    p <- plotly::add_markers(p, x = s$dif, y = s$item_id, name = sides[[kind]]$name, error_x = bars,
                             marker = list(color = sides[[kind]]$col, size = MARKER_SIZE + 1),
                             hovertext = hover, hoverinfo = "text")
  }
  y_axis <- list(title = LAB$item, type = "category", categoryorder = "array",
                 categoryarray = d$item_id)
  p |>
    plotly::layout(xaxis = list(title = LAB$dif_axis, zeroline = FALSE), yaxis = y_axis,
                   shapes = cut_lines(c(-cut, 0, cut), "x")) |>
    kits_style(sprintf(LAB$title_dif, ref, focal, cut))
}

# 圖形的模型檢查：兩組是組對組（|Wald z| > 1.96 橘），多組是每組對全部一起的校準，一組一格
dif_steps_plot <- function(x, items) {
  full <- as.data.frame(attr(x, "full"))
  shown <- pick_items(unique(full$item_id), items)
  steps <- as.data.frame(attr(x, "steps"))
  steps <- steps[steps$item_id %in% shown, ]
  groups <- attr(x, "lrtest")$groups$group
  step_id <- function(d) sprintf("%s:%s", d$item_id, d$item_score)   # 沒有列時是空的（不是 ":"）
  if (length(groups) == 2) {
    d <- as.data.frame(attr(x, "diff"))
    d <- d[d$item_id %in% shown, ]
    s1 <- steps[steps$group == groups[1], ]
    s2 <- steps[steps$group == groups[2], ]
    j1 <- match(step_id(d), step_id(s1))
    j2 <- match(step_id(d), step_id(s2))
    pts <- data.frame(step = step_id(d), x = s1$beta[j1], x_se = s1$se[j1], y = s2$beta[j2],
                      y_se = s2$se[j2], flag = abs(d$z) > 1.96, z = d$z)
    p <- dif_panel(pts, paste("Beta in group", groups[1]), paste("Beta in group", groups[2]),
                   c("|Wald z| <= 1.96", "|Wald z| > 1.96"))
    return(kits_style(p, sprintf(paste0("Graphical model check: beta in group %s against group ",
                                        "%s<br>bars = +/- 1.96 SE; gray line = identity"),
                                 groups[2], groups[1])))
  }
  pts <- map(groups, \(g) {
    s <- steps[steps$group == g, ]
    j <- match(step_id(s), step_id(full))
    data.frame(step = step_id(s), x = full$beta[j], x_se = full$se[j], y = s$beta, y_se = s$se) |>
      mutate(z = (y - x) / y_se, flag = abs(z) > 1.96)
  })
  # 圖例：每一種點只在第一個有它的格子裡列出來
  first_of <- function(k) which(map_lgl(pts, \(d) any(d$flag %in% (k == 2))))[1]
  legend_at <- c(first_of(1), first_of(2))
  panels <- map(seq_along(groups), \(i) {
    dif_panel(pts[[i]], "Beta, joint calibration", paste("Beta, group", groups[i]),
              c("Within 1.96 SE of the joint beta", "Beyond 1.96 SE"),
              legend_kind = (i == legend_at) %in% TRUE, identity_legend = i == 1)
  })
  n_col <- min(3, length(groups))
  plotly::subplot(panels, nrows = ceiling(length(groups) / n_col), titleX = TRUE, titleY = TRUE,
                  margin = c(0.04, 0.02, 0.06, 0.06)) |>
    kits_style(paste0("Graphical model check: each group against the joint calibration<br>",
                      "bars = +/- 1.96 SE; gray line = identity"))
}

# 一格：identity 線（灰）和每個步驟（藍；|z| > 1.96 橘），x、y 都有 +/- 1.96 SE 的線。
# legend_kind：兩種點（藍、橘）各自要不要列在圖例；identity_legend：identity 線要不要列
dif_panel <- function(pts, x_title, y_title, names, legend_kind = c(TRUE, TRUE),
                      identity_legend = TRUE) {
  ends <- c(pts$x - 1.96 * pts$x_se, pts$y - 1.96 * pts$y_se, pts$x + 1.96 * pts$x_se,
            pts$y + 1.96 * pts$y_se)
  ends <- if (any(is.finite(ends))) range(ends, na.rm = TRUE) else c(-1, 1)   # 這格沒有步驟
  p <- plotly::plot_ly() |>
    plotly::add_lines(x = ends, y = ends, name = "Identity", showlegend = identity_legend,
                      legendgroup = "id", line = list(color = PAL$gray, width = LINE_W$main),
                      hoverinfo = "skip")
  for (k in 1:2) {
    s <- pts[pts$flag == (k == 2) & !is.na(pts$flag), ]
    if (nrow(s) == 0) next
    hover <- sprintf("%s: x %.2f (SE %.2f), y %.2f (SE %.2f), z %.1f", s$step, s$x, s$x_se, s$y,
                     s$y_se, s$z)
    col <- c(PAL$blue, PAL$orange)[k]
    p <- plotly::add_markers(p, x = s$x, y = s$y, name = names[k], legendgroup = names[k],
                             showlegend = legend_kind[k],
                             error_x = list(type = "data", array = 1.96 * s$x_se, color = col),
                             error_y = list(type = "data", array = 1.96 * s$y_se, color = col),
                             marker = list(color = col, symbol = c("circle", "x")[k],
                                           size = MARKER_SIZE),
                             hovertext = hover, hoverinfo = "text")
  }
  plotly::layout(p, xaxis = list(title = x_title), yaxis = list(title = y_title))
}

# 虛線：在 x（或 y）的這些值上，貫穿整張圖
cut_lines <- function(at, axis) {
  map(at, \(v) {
    line <- list(color = PAL$gray, width = LINE_W$thin, dash = if (v == 0) "solid" else "dash")
    if (axis == "x") {
      return(list(type = "line", x0 = v, x1 = v, y0 = 0, y1 = 1, xref = "x", yref = "paper",
                  line = line))
    }
    list(type = "line", x0 = 0, x1 = 1, y0 = v, y1 = v, xref = "paper", yref = "y", line = line)
  })
}


#' @rdname kits_itemfit
#' @param x A `kits_itemfit`.
#' @param type What `plot()` draws: `"both"` (the default; infit and outfit), `"infit"` or
#'   `"outfit"`.
#' @param items Which items to draw (item ids); by default all of them.
#' @param range The acceptable range of the mean squares, shaded (default 0.6 to 1.4); the points
#'   are coloured by the flag against it (circles infit, squares outfit).
#' @examples
#'
#' # Infit and outfit of each item, the acceptable range shaded
#' p <- plot(fit)
#' p <- plot(fit, type = "outfit", range = c(0.5, 1.5))
#' if (interactive()) p
#' @export
plot.kits_itemfit <- function(x, type = c("both", "infit", "outfit"), items = NULL,
                              range = c(0.6, 1.4), ...) {
  type <- match_option(type, c("both", "infit", "outfit"), "type")
  abort_if(!is.numeric(range) || length(range) != 2 || any(!is.finite(range)) || range[1] >= range[2],
           "`range` must be two increasing numbers, e.g. c(0.6, 1.4).")
  shown <- pick_items(x$item_id, items)
  d <- filter(as.data.frame(x), item_id %in% shown)
  kinds <- if (type == "both") c("infit", "outfit") else type
  flag <- itemfit_flag(d$infit, d$outfit, range)   # 點的顏色：這張圖的範圍算的 flag
  shown_flags <- intersect(c("underfit", "overfit", "ok"), flag)
  cols <- c(underfit = PAL$orange, overfit = PAL$blue, ok = PAL$gray)
  symbols <- c(infit = "circle", outfit = "square")
  p <- plotly::plot_ly()
  for (st in kinds) {
    for (fl in shown_flags) {   # 一種統計量、一種 flag 一條：圖例是例如 infit: underfit
      e <- d[flag %in% fl, ]
      hover <- sprintf("%s: %s %.2f (n = %d; %s)", e$item_id, st, e[[st]], as.integer(e$n), fl)
      marker <- list(color = cols[[fl]], symbol = symbols[[st]], size = MARKER_SIZE + 1)
      p <- plotly::add_markers(p, x = e$item_id, y = e[[st]], name = paste0(st, ": ", fl),
                               marker = marker, hovertext = hover, hoverinfo = "text")
    }
  }
  band <- list(type = "rect", x0 = 0, x1 = 1, y0 = range[1], y1 = range[2], xref = "paper",
               yref = "y",
               fillcolor = "rgba(127,127,127,0.15)", line = list(width = 0), layer = "below")
  x_axis <- list(title = LAB$item, type = "category", categoryorder = "array",
                 categoryarray = d$item_id)
  p |>
    plotly::layout(xaxis = x_axis, yaxis = list(title = LAB$mean_square), shapes = list(band)) |>
    kits_style(sprintf(LAB$title_itemfit, itemfit_label(x), range[1], range[2]))
}


#' @rdname kits_drift
#' @param type What `plot()` draws: `"scatter"` (the default; the new threshold of each step against
#'   the old one, with the line of the judged steps) or `"difference"` (the distance of each step
#'   from that line).
#' @param items Which items to draw (item ids); by default all of them.
#' @examples
#'
#' # New against old thresholds; the steps of the drifted item in orange
#' p <- plot(drift)
#' p <- plot(drift, type = "difference")
#' if (interactive()) p
#' @export
plot.kits_drift <- function(x, type = c("scatter", "difference"), items = NULL, ...) {
  type <- match_option(type, c("scatter", "difference"), "type")
  shown <- pick_items(x$items$item_id, items)
  status <- as.data.frame(x$items) |>
    transmute(item_id, status = case_when(is.na(drift) ~ "not judged", drift ~ "drift",
                                          suggest_fix ~ "fix", .default = "other"))
  d <- as.data.frame(x$points) |>
    filter(item_id %in% shown) |>
    left_join(status, by = "item_id") |>
    mutate(step = paste0(item_id, ":", item_score))
  looks <- list(fix = list(col = PAL$blue, symbol = "circle", name = "Suggested to fix"),
                drift = list(col = PAL$orange, symbol = "x", name = "Drift (step dropped)"),
                other = list(col = PAL$black, symbol = "square", name = "Not suggested"),
                `not judged` = list(col = PAL$gray, symbol = "circle-open", name = "Not judged"))
  p <- plotly::plot_ly()
  if (type == "scatter") {
    ends <- range(c(d$thr_new, d$thr_ref), na.rm = TRUE)
    p <- plotly::add_lines(p, x = ends, y = line_apply(x, ends),
                           name = sprintf(LAB$final_line, x$intercept, x$slope),
                           line = list(color = PAL$black, width = LINE_W$main), hoverinfo = "skip")
  }
  for (st in names(looks)) {
    s <- filter(d, status == st)
    if (nrow(s) == 0) next
    hover <- sprintf("%s: new %.2f, old %.2f, old - line %.2f", s$step, s$thr_new, s$thr_ref,
                     s$resid)
    marker <- list(color = looks[[st]]$col, symbol = looks[[st]]$symbol, size = MARKER_SIZE + 1)
    xs <- if (type == "scatter") s$thr_new else s$step
    ys <- if (type == "scatter") s$thr_ref else s$resid
    p <- plotly::add_markers(p, x = xs, y = ys, name = looks[[st]]$name, marker = marker,
                             hovertext = hover, hoverinfo = "text")
  }
  if (type == "scatter") {
    return(plotly::layout(p, xaxis = list(title = LAB$drift_new),
                          yaxis = list(title = LAB$drift_old)) |>
             kits_style(LAB$title_drift))
  }
  x_axis <- list(title = LAB$step, type = "category", categoryorder = "array",
                 categoryarray = d$step)
  plotly::layout(p, xaxis = x_axis, yaxis = list(title = LAB$drift_resid),
                 shapes = cut_lines(0, "y")) |>
    kits_style(LAB$title_drift_diff)
}


#' @rdname kits_norm_screen
#' @param x A `kits_norm_screen`.
#' @examples
#'
#' # The share at the lowest score of each wave; the waves left out in orange
#' p <- plot(screen)
#' if (interactive()) p
#' @export
plot.kits_norm_screen <- function(x, ...) {
  d <- as.data.frame(x$waves)
  has_cohort <- "cohort" %in% names(d) && !all(is.na(d$cohort))
  d$wave <- if (has_cohort) paste(d$cohort, d$wave_m) else as.character(d$wave_m)
  looks <- list(kept = list(col = PAL$blue, name = "In the norms"),
                out = list(col = PAL$orange, name = "Left out"))
  p <- plotly::plot_ly()
  for (st in names(looks)) {
    s <- filter(d, normed == (st == "kept"))
    if (nrow(s) == 0) next
    hover <- sprintf("wave %s: %.1f%% at the lowest score, n = %d", s$wave, s$floor_pct,
                     as.integer(s$n))
    p <- plotly::add_bars(p, x = s$wave, y = s$floor_pct, name = looks[[st]]$name,
                          marker = list(color = looks[[st]]$col), hovertext = hover,
                          hoverinfo = "text")
  }
  x_axis <- list(title = LAB$wave_axis, type = "category", categoryorder = "array",
                 categoryarray = d$wave)
  plotly::layout(p, xaxis = x_axis, yaxis = list(title = LAB$floor_pct, rangemode = "tozero"),
                 shapes = cut_lines(x$cut, "y")) |>
    kits_style(sprintf(LAB$title_screen, x$cut))
}
