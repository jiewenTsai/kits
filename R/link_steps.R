# kits_linking() 的四步，一步一個匯出的函數（讓學生一步一步跑；kits_algorithm(link) 列出它們）：
#   kits_link_points()：點 -> kits_link_screen()：篩一次 -> kits_link_line()：直線 -> kits_link_error()：連結誤差
# 和 kits_linking() 是同一段程式（points_step()、screen_step()、line_step()、error_step()），所以數字相同。
# 自己準備的點（data.frame）從第 2 步開始。

# ---- 1 點 ------------------------------------------------------------------------------------------------------

# 兩個校準共同的步驟（一個步驟一個點），或錨題（一題一個點）。連結的設定（kind、step、key、兩端的量尺、題數）在 attr(, "link")
points_step <- function(from, to, step, to_scale) {
  need_number(step, "step", "NULL or one whole number (1 = the milestone step)", lower = 1, whole = TRUE, finite = FALSE,
              null = TRUE)
  if (!is.null(to_scale)) need_scale_arg(to_scale, "to_scale")
  is_domains <- is.list(from) && !is.data.frame(from) && !inherits(from, "kits_item") && !inherits(from, "kits")
  from_scale <- if (is_domains) map_chr(from, scale_of) else scale_of(from)
  need_link_target(to)
  kind <- if (is_anchors(to)) "anchors" else "calibration"
  if (kind == "anchors") {
    if (is.null(to_scale)) to_scale <- "D"
    if (is.null(step)) step <- 1
  } else {
    to_scale <- target_scale(to, to_scale)
  }
  check_link_scales(from_scale, to_scale, kind)
  steps <- as_steps(from, "from", domains = TRUE)
  need_steps(steps, "from")
  abort_if("domain" %in% names(steps) && anyNA(steps$domain), "`from$domain` has NA; give every item a domain.")
  pts <- if (kind == "anchors") anchor_points(steps, to, step) else common_points(steps, to, step)
  out <- kits_table(pts, "kits_link_points")
  attr(out, "link") <- list(kind = kind, step = if (is.null(step)) "all" else step,
                            no_anchors = domains_without_anchors(steps, pts, kind),
                            key = if (kind == "anchors") anchor_key(to),
                            from_scale = from_scale, to_scale = to_scale, n_items = length(unique(pts$item)),
                            items_total = if (kind == "calibration") items_total(steps, to))
  out
}

# 接到 D 時沒有錨題的領域（例如 KIT 的社會、情緒）：用同一條線（有錨題的領域估的），分數是「D 單位的 KIT 量尺分數」，
# 不是 GSED 的 D。沒有時 NULL
domains_without_anchors <- function(steps, pts, kind) {
  if (kind != "anchors" || !"domain" %in% names(steps)) {
    return(NULL)
  }
  none <- setdiff(unique(steps$domain), unique(pts$domain))
  if (length(none) == 0) NULL else none
}

# 兩端在同一條量尺上時停下來
check_link_scales <- function(from_scale, to_scale, kind) {
  if (is.na(to_scale) || !to_scale %in% from_scale) {
    return(invisible())
  }
  abort_if(kind == "calibration",
           sprintf("`from` and `to` are both on scale %s. A link between two calibrations joins two scales: ", to_scale),
           "give the calibrations different labels, e.g. kits_item(X36, scale = \"lan M36\").", class = "kits_error_scale")
  kits_abort(sprintf("`from` is already on scale %s; a link to %s joins two different scales.", to_scale, to_scale),
             class = "kits_error_scale")
}

# 兩個校準各有幾題（共同題的比例用，Kolen & Brennan p. 288 的 20%）
items_total <- function(steps, to) {
  c(from = length(unique(steps$item_id)), to = length(unique(as_steps(to, "to")$item_id)))
}

# ---- 2 篩一次 ---------------------------------------------------------------------------------------------------

# 所有的點的直線（intercept0、slope0），殘差 e，|e| > cut x 尺度的點標記、丟掉（unit = "item" 時整題丟掉）
screen_step <- function(points, cut, spread, unit) {
  link <- attr(points, "link")
  pts <- as.data.frame(points)
  need_point_se(pts, spread, link$kind)
  warn_negative(pts$from, pts$to)
  screened <- algo_screen(pts$from, pts$to, pts$item, cut, spread, unit, if (spread == "se") pts$se_from, pts$se_to)
  out <- pts |>
    mutate(resid = screened$resid, z = screened$z, flagged = screened$flagged, dropped = screened$dropped) |>
    kits_table("kits_link_screen")
  attr(out, "link") <- c(link, list(cut = cut, spread = spread, unit = unit, rule = screened$rule,
                                    intercept0 = screened$intercept0, slope0 = screened$slope0))
  out
}

# spread = "se" 要每個點的 SE：兩個校準之間兩邊都要，錨題（GSED 的 tau 沒有 SE）只要 from 的
need_point_se <- function(pts, spread, kind) {
  if (spread != "se") {
    return(invisible())
  }
  no_from <- is.null(pts$se_from) || anyNA(pts$se_from)
  no_to <- identical(kind, "calibration") && (is.null(pts$se_to) || anyNA(pts$se_to))
  abort_if(is.null(kind) && no_from,
           "`spread = \"se\"` needs the SE of each point: give `se_from` (and `se_to` when the target has SEs).")
  if (no_from || no_to) {
    whose <- if (identical(kind, "calibration")) "both calibrations" else "`from`"
    kits_abort("`spread = \"se\"` needs the SE of the thresholds (`thr_se`) of ", whose,
               "; give kits_item() results or tables with thr_se.", class = "kits_error_input")
  }
}

# from 和 to 負相關：mean/sigma 的斜率一定是正的，點可能沒有對好
warn_negative <- function(from, to) {
  usable <- length(from) >= 3 && is.finite(stats::sd(from)) && is.finite(stats::sd(to)) && stats::sd(from) > 0 &&
    stats::sd(to) > 0
  if (usable && stats::cor(from, to) < 0) {
    kits_warn("`from` and `to` are negatively related (r = ", round(stats::cor(from, to), 2),
              "). mean/sigma always gives a positive slope; check that the points are matched.", class = "kits_warning_scale")
  }
}

# ---- 3 直線 ------------------------------------------------------------------------------------------------------

# 留下的點的 mean/sigma：slope = SD(to) / SD(from)，intercept = mean(to) - slope mean(from)。結果是 kits_link
# （還沒有連結誤差）
line_step <- function(screened, call = NULL) {
  st <- attr(screened, "link")
  pts <- as.data.frame(screened)
  kept <- !pts$dropped
  final <- algo_mean_sigma(pts$from[kept], pts$to[kept])
  points <- data.frame(id = pts$id, item = pts$item, from = pts$from, to = pts$to, resid = pts$resid,
                       resid_final = line_resid(final, pts$from, pts$to), z = pts$z,
                       flagged = pts$flagged,
                       dropped = pts$dropped)
  if (st$spread == "se") {   # 每個點的 SE 留在 points 裡（教學程式 link.R 也用）
    points$se_from <- st$rule$se_from
    points$se_to <- st$rule$se_to
  }
  if (identical(st$kind, "anchors")) points$step <- as.integer(pts$step)
  from_scale <- if (is.null(st$from_scale)) NA_character_ else st$from_scale
  to_scale <- if (is.null(st$to_scale)) NA_character_ else st$to_scale
  link <- structure(list(intercept = final$intercept, slope = final$slope,
                         intercept_error = NA_real_,
                         slope_error = NA_real_, n_failed = NA_integer_, intercept0 = st$intercept0,
                         slope0 = st$slope0, n = nrow(pts), dropped = pts$id[pts$dropped],
                         points = kits_table(points), settings = line_settings(st),
                         from_scale = from_scale, to_scale = to_scale),
                    class = "kits_link")
  link$info <- make_info(call, link$settings, values_fingerprint(pts$from, pts$to, pts$id))
  link$n_items <- st$n_items
  link$items_total <- st$items_total
  if ("domain" %in% names(pts)) link <- add_domains(link, pts$domain)
  link$diagnostics <- link_diagnostics(link)
  link_in_order(link)
}

# 連結的設定：篩點的規則，重抽（還沒做：n_boot = 0），點從哪裡來（kits_link_points() 的點）
line_settings <- function(st) {
  settings <- list(cut = st$cut, n_boot = 0, seed = 1, spread = st$spread, unit = st$unit)
  if (!is.null(st$kind)) settings <- c(settings, list(kind = st$kind, step = st$step))
  if (!is.null(st$key)) settings$key <- st$key   # kits_anchors() 的 key
  if (!is.null(st$no_anchors)) settings$no_anchors <- st$no_anchors   # 沒有錨題、用同一條線的領域
  settings
}

# kits_link 的元素照固定的順序（和 0.17.0 相同；0.23.0 起 a、b 叫 intercept、slope，ab_cov 是 vcov）。
# 沒有元素以 a 或 b 開頭，所以舊的 link$a、link$b 不會被部分比對抓到別的元素
LINK_ORDER <- c("intercept", "slope", "intercept_error", "slope_error", "n_failed", "intercept0",
                "slope0",
                "n", "dropped", "points", "settings", "vcov", "draws", "from_scale", "to_scale",
                "diagnostics", "info", "n_items", "items_total", "domains")

link_in_order <- function(link) {
  structure(unclass(link)[intersect(LINK_ORDER, names(link))], class = "kits_link")
}

# ---- 4 連結誤差 ---------------------------------------------------------------------------------------------------

# 整題重抽 n_boot 次，每次重做第 2、3 步；截距、斜率在重抽間的 SD 是連結誤差。n_boot = 0 時只記下設定
error_step <- function(link, n_boot, seed) {
  link[c("intercept_error", "slope_error", "n_failed")] <- list(NA_real_, NA_real_, NA_integer_)
  link[c("vcov", "draws")] <- NULL
  link$settings$n_boot <- n_boot
  link$settings["seed"] <- list(seed)   # seed 可以是 NULL
  pts <- link$points
  if (n_boot > 0 && length(unique(pts$item)) < 2) {
    kits_warn("All points belong to one item, so the items cannot be resampled; no linking error.", class = "kits_warning_scale")
  } else if (n_boot > 0) {
    st <- link$settings
    boot <- algo_screen(pts$from, pts$to, pts$item, st$cut, st$spread, st$unit, pts$se_from, pts$se_to) |>
      algo_link_error(n_boot, seed)
    errors <- c("intercept_error", "slope_error", "n_failed")
    link[errors] <- boot[errors]
    link$vcov <- boot$vcov
    link$draws <- boot$draws
  }
  link$info$settings <- link$settings
  link$diagnostics <- link_diagnostics(link)
  link_in_order(link)
}


#' A link step by step: points, screening, line, linking error
#'
#' The four steps of [kits_linking()] as functions of their own:
#' `kits_link_points()` (the points), `kits_link_screen()` (one round of screening),
#' `kits_link_line()` (the mean/sigma line) and `kits_link_error()` (the item-resampling linking
#' error). Piped in this order they give the same link as `kits_linking(..., n_boot = )`.
#' Details: `kits_algorithm("link")` and tutorial_algorithms/03_linking.html.
#'
#' @inheritParams kits_linking
#' @param points A [kits_link_points()] result, or a data.frame with `id`, `item`, `from`, `to`
#'   (and `se_from`, `se_to` for `spread = "se"`), for example points from elsewhere.
#' @param from_scale,to_scale In `kits_link_screen()`: the labels of the two scales of your own
#'   points (optional; only with a data.frame).
#'   In `kits_link_points()`: `to_scale` as in [kits_linking()].
#' @param screened A `kits_link_screen()` result.
#' @param link A link: a `kits_link_line()` result (or any [kits_linking()] result, whose linking
#'   error is then computed again).
#' @param n_boot Number of item resamples; 0 removes the linking error.
#'   200 is enough to see its size; use more (e.g. 2000) for reporting.
#' @return
#' * `kits_link_points()`: a `kits_link_points`, one row per point, `id`, `item`, `from`, `to`
#'   (and `step`, `domain`, `se_from`, `se_to` when available).
#' * `kits_link_screen()`: a `kits_link_screen`, the points with `resid`, `z`, `flagged` and
#'   `dropped`.
#' * `kits_link_line()`: a `kits_link` (see [kits_linking()]) with NA linking errors.
#' * `kits_link_error()`: the link with `intercept_error`, `slope_error`, `n_failed`, `draws` and
#'   `vcov`.
#' @seealso [kits_linking()] (all four steps in one call), [kits_convert()].
#' @examples
#' # Two calibrations with common items a3, a4, a5
#' betas <- list(a1 = -2, a2 = -1, a3 = c(-0.5, 0.5), a4 = 0, a5 = c(0.5, 1.5), a6 = 2, a7 = 3)
#' young <- kits_sim(rnorm(500, 0), betas[1:5], seed = 1)
#' old <- kits_sim(rnorm(500, 2), betas[3:7], prefix = "q", seed = 2)
#' item_young <- kits_item(young, scale = "young cohort")
#' item_old <- kits_item(old, scale = "old cohort")
#'
#' points <- kits_link_points(item_old, item_young)   # 1 the points
#' screened <- kits_link_screen(points)               # 2 one round of screening
#' link <- kits_link_line(screened)                   # 3 the line
#' link <- kits_link_error(link, n_boot = 200, seed = 1)  # 4 the linking error
#' link
#'
#' # Your own points: point 8 is shifted and dropped
#' own <- data.frame(id = paste0("p", 1:8), item = paste0("p", 1:8),
#'                   from = c(-2, -1, 0, 1, 2, 3, 0.5, 1.5),
#'                   to = c(-1, 0, 1, 2, 3, 4, 1.5, 5.5))
#' kits_link_screen(own, from_scale = "x", to_scale = "y", spread = "mad")$dropped
#' @name kits_link_steps
NULL

#' @rdname kits_link_steps
#' @export
kits_link_points <- function(from, to, step = NULL, to_scale = NULL) {
  points_step(from, to, step, to_scale)
}

#' @rdname kits_link_steps
#' @export
kits_link_screen <- function(points, cut = 2, spread = c("sd", "mad", "se"),
                             unit = c("step", "item"), from_scale = NULL, to_scale = NULL) {
  spread <- match_option(spread, c("sd", "mad", "se"), "spread")
  unit <- match_option(unit, c("step", "item"), "unit")
  need_cols(points, c("id", "item", "from", "to"), "points")
  need_link_values(points$from, points$to)
  need_cut(cut)
  points <- own_points_scales(points, from_scale, to_scale)
  screen_step(points, cut, spread, unit)
}

# 自己的點（data.frame）的兩端：from_scale、to_scale，沒給時 NA。kits_link_points() 的點已經記下兩端
own_points_scales <- function(points, from_scale, to_scale) {
  given <- !is.null(from_scale) || !is.null(to_scale)
  if (!is.null(attr(points, "link"))) {
    abort_if(given, "The points of kits_link_points() record their scales; ",
             "give `from_scale` and `to_scale` only with your own points (a data.frame).")
    return(points)
  }
  if (!given) {
    return(points)
  }
  from_scale <- if (is.null(from_scale)) NA else from_scale
  to_scale <- if (is.null(to_scale)) NA else to_scale
  attr(points, "link") <- link_ends(from_scale, to_scale)
  points
}

#' @rdname kits_link_steps
#' @export
kits_link_line <- function(screened) {
  abort_if(!inherits(screened, "kits_link_screen") || is.null(attr(screened, "link")),
           "`screened` must be a kits_link_screen() result, e.g. kits_link_line(kits_link_screen(points)).")
  line_step(screened, call_text(sys.call()))
}

#' @rdname kits_link_steps
#' @export
kits_link_error <- function(link, n_boot = 200, seed = 1) {
  abort_if(!inherits(link, "kits_link"), "`link` must be a link: a kits_link_line() or kits_linking() result.")
  need_n_boot(n_boot)
  need_seed(seed)
  error_step(link, n_boot, seed)
}
