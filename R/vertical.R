# 垂直量尺的統計（Kolen & Brennan, 2014, pp. 461-469）：每個世代 x 波次的平均、SD、中位數、真分數的 SD，
# 和前一個波次比的成長、效果量（Eq. 9.47）、分佈的分離（P(上一波 < 這一波)）；兩個世代都有的波次另外比較。

# 這一波比前一波高的機率（相同時算一半）：Mann-Whitney，用中位秩
p_superior <- function(upper, lower) {
  r <- rank(c(upper, lower))
  nu <- length(upper)
  (sum(r[seq_len(nu)]) - nu * (nu + 1) / 2) / (nu * length(lower))
}

# 一組的平均的連結誤差：每個孩子的連結誤差都來自同一串連結，所以平均分數的連結誤差是 se_link 在 theta 的平均上的值。
# * 分數記得自己的連結（kits_subj() 的分數用 kits_convert() 換過，attr(, "raw_theta") 記下每一步）時，直接算：
#   link_moments(links, x0)，x0 是平均換回起點量尺的值（和 kits_convert() 的 se_link 同一個計算）。
# * 否則從每個孩子的 se_link 反推（近似）：se_link^2 是 theta 的二次式，用最小平方配二次式，取在平均上的值
#   （不同的 theta 少於 3 個時用 se_link^2 的平均，比較保守）。
mean_link_se <- function(theta, se_link, links = NULL) {
  if (!is.null(links) && length(theta) > 0 && all(is.finite(theta)) && all(is.finite(se_link))) {
    x0 <- algo_line(links, mean(theta), invert = TRUE)
    return(sqrt(link_moments(links, x0)$var))
  }
  ok <- is.finite(theta) & is.finite(se_link)
  t <- theta[ok] - mean(theta[ok])
  v <- se_link[ok]^2
  if (length(v) == 0) {
    return(NA_real_)
  }
  if (length(unique(t)) < 3) {
    return(sqrt(mean(v)))
  }
  fit <- stats::lm.fit(cbind(1, t, t^2), v)
  sqrt(max(unname(fit$coefficients[1]), 0))
}


#' Statistics of a vertical scale: growth, variability and separation by wave
#'
#' For scores of every cohort on one scale (for example D, after [kits_convert()]), gives per
#' cohort x wave the location and spread (with `sd_true`, the SD without measurement error), the
#' growth from the previous wave, and its separation (`effect_size` of Yen, 1986; `p_sup`), the
#' properties of Kolen and Brennan (2014). `attr(, "cohorts")` compares cohorts at shared waves.
#' `plot()` draws the mean +/- SD (or the median and quartiles) at the median age of each wave.
#' Details: tutorial_algorithms/06_scale.html.
#'
#' @param s Scores on one scale: a [kits_subj()] result (converted) with `wave_m` and optionally `cohort`, `age_m`,
#'   `se`, `se_link` and `flag`; or a data.frame with `theta` and `wave_m` and any of those.
#' @param extremes Use the children with a zero or a perfect score (`flag`) in `sd_true`? `TRUE` (the default) uses
#'   every child; `FALSE` leaves them out of both the SD and the mean error variance.
#' @return A table (a data.frame) with one row per cohort x wave, in the order of cohort and wave:
#'   `cohort` (when given), `wave_m`, `n`, `age_med` (median age, when `age_m` is given), `mean`, `sd`, `median`,
#'   `P25`, `P75`,
#'   `sd_true` (NA without `se`, or with fewer than 2 children without a zero or perfect score),
#'   `growth`, `effect_size`, `p_sup` (NA for the first wave of a cohort),
#'   and `se_link` (the linking error of the mean, when the scores have `se_link`).
#'   `attr(, "cohorts")` compares cohorts at shared waves; diagnostics and [kits_info()] too.
#' @references
#' Kolen, M. J., & Brennan, R. L. (2014). *Test equating, scaling, and linking* (3rd ed.). Springer.
#' \doi{10.1007/978-1-4939-0317-7}
#'
#' Yen, W. M. (1986). The choice of scale for educational measurement: An IRT perspective.
#' *Journal of Educational Measurement, 23*(4), 299-325. \doi{10.1111/j.1745-3984.1986.tb00252.x}
#' @examples
#' library(dplyr)
#' set.seed(1)
#' s <- tibble(cohort = rep(c("M03", "M36"), c(600, 400)),
#'             wave_m = c(rep(c(12, 24, 36), each = 200), rep(c(36, 48), each = 200))) |>
#'   mutate(theta = 40 + 8 * log(wave_m) + rnorm(n(), 0, 4), se = 2)
#' v <- kits_vertical(s)
#' v
#' attr(v, "cohorts")
#' p <- plot(v, type = "median")
#' if (interactive()) p
#'
#' # Two domains stacked: one panel per domain
#' both <- list(lan = v, motor = kits_vertical(mutate(s, theta = theta + 5))) |>
#'   purrr::list_rbind(names_to = "domain")
#' p <- plot(both)
#' @export
kits_vertical <- function(s, extremes = TRUE) {
  call <- call_text(sys.call())
  abort_if(!is.data.frame(s), "`s` must be a kits_subj() result or a data.frame with theta and wave_m.")
  need_cols(s, c("theta", "wave_m"), "s")
  abort_if(!is.numeric(s$theta) || !is.numeric(s$wave_m), "`s$theta` and `s$wave_m` must be numbers.")
  abort_if(!isTRUE(extremes) && !isFALSE(extremes), "`extremes` must be TRUE or FALSE.")
  d <- as.data.frame(s) |>
    filter(is.finite(theta), !is.na(wave_m))
  has_cohort <- "cohort" %in% names(d)
  cohort <- if (has_cohort) as.character(d$cohort) else rep(NA_character_, nrow(d))
  d <- mutate(d, cohort = cohort, extreme = !extremes & flag_extreme(d))
  chains <- score_link_chains(s, d)
  d$.chain <- chains$row
  tab <- d |>
    summarise(wave_stats(pick(everything()), chains$links), thetas = list(theta),
              .by = c(cohort, wave_m)) |>
    arrange(sort_rank(cohort), wave_m) |>
    mutate(growth = mean - lag(mean),
           effect_size = growth / sqrt((sd^2 + lag(sd)^2) / 2),
           p_sup = map2_dbl(thetas, lag(thetas), \(th, prev) if (is.null(prev)) NA_real_ else p_superior(th, prev)),
           .by = cohort) |>
    select(-thetas) |>
    relocate(se_link, .after = p_sup)
  cohorts <- vertical_cohorts(tab)
  if (!has_cohort) tab$cohort <- NULL
  if (!"age_m" %in% names(d)) tab$age_med <- NULL
  if (!"se_link" %in% names(d)) tab$se_link <- NULL
  out <- kits_table(tab, "kits_vertical")
  attr(out, "cohorts") <- cohorts
  attr(out, "scale") <- scale_of(s)
  attr(out, "diagnostics") <- vertical_diagnostics(out, cohorts, "se" %in% names(d))
  attr(out, "info") <- make_info(call, list(scale = scale_of(s), extremes = extremes), values_fingerprint(d$theta, d$wave_m,
      cohort))
  out
}

# 0 分和滿分的人次（kits_subj 的 flag）：他們的 WLE se 很大，扣掉測量誤差的 SD 不用他們。沒有 flag 時都不算。
flag_extreme <- function(d) {
  if (!"flag" %in% names(d)) {
    return(rep(FALSE, nrow(d)))
  }
  d$flag %in% c("zero", "perfect")
}

# 每一列分數換算時經過的連結（kits_subj() 的分數記在 attr(, "raw_theta")：每個校準的人次和之後的換算）。
# 回傳 links（每個來源一串連結；有一條沒有重抽，或沒有換算時 NULL）和 row（每一列是哪個來源；不知道時 NA）
score_link_chains <- function(s, d) {
  sources <- attr(s, "raw_theta")
  none <- list(links = list(), row = rep(NA_integer_, nrow(d)))
  if (is.null(sources) || !"visit_id" %in% names(d) || !"se_link" %in% names(d)) {
    return(none)
  }
  links <- map(sources, \(src) {
    l <- unlist(map(src$steps, "links"), recursive = FALSE)
    if (length(l) == 0 || !all(map_lgl(l, link_has_draws))) NULL else l
  })
  owner <- sources |>
    imap(\(src, i) data.frame(visit_id = src$persons, source = i)) |>
    list_rbind()
  owner <- owner[!owner$visit_id %in% owner$visit_id[duplicated(owner$visit_id)], ]   # 兩個來源都有：不知道
  list(links = links, row = owner$source[match(as.character(d$visit_id), owner$visit_id)])
}

# 一個世代 x 波次的位置和分散：sd 用所有孩子；sd_true 只用不是 0 分或滿分的孩子（SD 和 mean(se^2) 都是）。
# links：每個來源的連結（score_link_chains()）；這組的列都來自同一個有連結的來源時，平均的連結誤差直接算
wave_stats <- function(g, links = list()) {
  src <- unique(g$.chain)
  chain <- if (length(src) == 1 && !is.na(src)) links[[src]]
  inner <- filter(g, !extreme)
  data.frame(n = nrow(g),
             age_med = if ("age_m" %in% names(g)) stats::median(g$age_m, na.rm = TRUE) else NA_real_,
             mean = mean(g$theta),
             sd = if (nrow(g) > 1) stats::sd(g$theta) else NA_real_,
             median = stats::median(g$theta),
             P25 = stats::quantile(g$theta, 0.25, names = FALSE),
             P75 = stats::quantile(g$theta, 0.75, names = FALSE),
             sd_true = if ("se" %in% names(g)) true_sd(inner$theta, inner$se) else NA_real_,
             se_link = if ("se_link" %in% names(g)) {
               mean_link_se(g$theta, g$se_link, chain)
             } else {
               NA_real_
             })
}

# 扣掉測量誤差的 SD：sqrt(max(var(theta) - mean(se^2), 0))；少於 2 人時 NA
true_sd <- function(theta, se) {
  if (length(theta) < 2) {
    return(NA_real_)
  }
  sqrt(max(stats::var(theta) - mean(se^2), 0))
}

# 兩個世代都有的波次：平均和 SD 的差（後一個世代減前一個）
vertical_cohorts <- function(tab) {
  if (all(is.na(tab$cohort))) {
    return(NULL)
  }
  pairs <- tab |>
    select(wave_m, cohort, n, mean, sd, se_link) |>
    mutate(cohort_1 = lag(cohort), n_1 = lag(n), mean_1 = lag(mean), sd_1 = lag(sd), se_link_1 = lag(se_link), .by = wave_m) |>
    filter(!is.na(n_1)) |>
    arrange(wave_m)
  if (nrow(pairs) == 0) {
    return(NULL)
  }
  out <- pairs |>
    transmute(wave_m, cohort_1, cohort_2 = cohort, n_1, n_2 = n,
              mean_diff = mean - mean_1, sd_diff = sd - sd_1, sd_ratio = sd / sd_1,
              se_diff = sqrt(sd_1^2 / n_1 + sd^2 / n),
              se_link = pair_link_se(se_link_1, se_link),
              z = mean_diff / sqrt(se_diff^2 + coalesce(se_link, 0)^2),
              z_link = mean_diff / se_link)
  if (all(is.na(out$se_link))) {
    out <- select(out, -se_link, -z_link)
  }
  kits_table(out)
}

# 兩個世代的平均之差的連結誤差。一個世代沒有 se_link（例如本來就在這條量尺上的參照世代）時當成 0；兩個都沒有時是 NA
pair_link_se <- function(se_1, se_2) {
  map2_dbl(se_1, se_2, \(a, b) if (is.na(a) && is.na(b)) NA_real_ else sqrt(sum(c(a, b)^2, na.rm = TRUE)))
}

vertical_diagnostics <- function(tab, cohorts, has_se) {
  cohort <- if (is.null(tab$cohort)) rep(NA_character_, nrow(tab)) else tab$cohort
  rows <- c(shrinkage_rows(tab, cohort, has_se), growth_rows(tab, cohort), zero_true_sd_rows(tab, cohort, has_se),
            cohort_rows(cohorts))
  rows <- compact(rows)
  if (length(rows) == 0) rows <- list(diag_row("waves", nrow(tab), "ok", "one wave per cohort: nothing to compare"))
  diag_table(rows)
}

# 量尺縮小（p. 468）：最大的波次的 sd_true（沒有 se 時 sd）和前面最大的比；每個世代一列，有好幾個世代時再加全部一起
shrinkage_rows <- function(tab, cohort, has_se) {
  spread <- if (has_se) tab$sd_true else tab$sd
  label <- if (has_se) "sd_true" else "sd"
  groups <- map(unique(cohort), \(co) which(if (is.na(co)) is.na(cohort) else cohort %in% co))
  names(groups) <- map_chr(unique(cohort), \(co) if (is.na(co)) "all waves" else co)
  if (length(groups) > 1) groups[["all cohorts"]] <- seq_len(nrow(tab))
  imap(groups, \(idx, what) shrinkage_row(idx[order(tab$wave_m[idx])], spread, label, what))
}

shrinkage_row <- function(idx, spread, label, what) {
  if (length(idx) < 2) {
    return(NULL)
  }
  last <- spread[idx[length(idx)]]
  earlier <- spread[idx[-length(idx)]]
  if (all(is.na(earlier))) {   # 例如每個波次只有一個孩子：沒有 SD
    return(NULL)
  }
  before <- max(earlier, na.rm = TRUE)
  if (!is.finite(last) || !is.finite(before) || before <= 0) {
    return(NULL)
  }
  ratio <- last / before
  shrinks <- ratio < 2 / 3
  status <- if (shrinks) "note" else "ok"
  message <- if (shrinks) {
    sprintf(paste("%s of the oldest wave is %.0f%% of the largest earlier one: the scale shrinks with age",
                  "(scale shrinkage, Kolen & Brennan p. 468); check the measurement error and the links"),
            label, 100 * ratio)
  } else {
    sprintf("%s of the oldest wave / the largest earlier %s", label, label)
  }
  diag_row(paste("scale shrinkage,", what), ratio, status, message)
}

# 兩個世代的平均差，以連結誤差為單位；超過 2 註記
link_z_row <- function(wave, zl) {
  status <- if (abs(zl) > 2) "note" else "ok"
  detail <- if (abs(zl) > 2) {
    paste(" (beyond 2 linking errors: the cohorts differ more than the link's error explains;",
          "a cohort effect or a link to check; Kolen & Brennan p. 333)")
  } else {
    ""
  }
  diag_row(sprintf("cohorts at wave %s vs linking error", wave), zl, status,
           sprintf("mean difference / linking error = %.2f%s", zl, detail))
}

# 平均比前一波低的波次
growth_rows <- function(tab, cohort) {
  down <- which(is.finite(tab$growth) & tab$growth < 0)
  map(down, \(i) diag_row(paste("mean below the previous wave,", vertical_label(cohort[i], tab$wave_m[i])), tab$growth[i],
                          "note", "the mean is lower than at the previous wave"))
}

# sd_true 是 0 的波次：測量誤差至少和觀察到的變異一樣大
zero_true_sd_rows <- function(tab, cohort, has_se) {
  zero <- which(has_se & is.finite(tab$sd_true) & tab$sd_true == 0)
  map(zero, \(i) diag_row(paste("sd_true is 0,", vertical_label(cohort[i], tab$wave_m[i])), 0, "note",
                          "the measurement error is at least the observed variance: the wave cannot separate children"))
}

# 兩個世代都有的波次：平均差的 z、和連結誤差比、SD 的比
cohort_rows <- function(cohorts) {
  if (is.null(cohorts)) {
    return(list())
  }
  rows <- map(seq_len(nrow(cohorts)), \(i) cohort_pair_rows(cohorts[i, ]))
  list_flatten(rows)
}

cohort_pair_rows <- function(pair) {
  wave <- format(pair$wave_m)
  z <- pair$z
  z_status <- if (abs(z) > 2) "note" else "ok"
  z_detail <- if (abs(z) > 2) " (beyond 2 SE: check the cohort link)" else ""
  z_message <- sprintf("%s minus %s: mean %.2f, SD %.2f; z = %.2f%s", pair$cohort_2, pair$cohort_1,
                       pair$mean_diff, pair$sd_diff, z, z_detail)
  rows <- list(diag_row(sprintf("cohorts at wave %s", wave), z, z_status, z_message))
  zl <- pair$z_link
  if (!is.null(zl) && is.finite(zl)) {
    rows <- c(rows, list(link_z_row(wave, zl)))
  }
  r <- pair$sd_ratio
  if (is.finite(r) && (r < 0.8 || r > 1.2)) {
    rows <- c(rows, list(diag_row(sprintf("cohorts at wave %s, SD ratio", wave), r, "note",
      sprintf(paste("SD of %s / SD of %s = %.2f, outside 0.8 to 1.2: the slope of the cohort link may be off",
                    "(Kolen & Brennan p. 301 use this range for group differences)"), pair$cohort_2, pair$cohort_1, r))))
  }
  rows
}

vertical_label <- function(cohort, wave) if (is.na(cohort)) sprintf("wave %s", format(wave)) else sprintf("%s wave %s", cohort,
    format(wave))

#' @rdname kits_vertical
#' @param x A `kits_vertical`.
#' @param digits Decimals to print.
#' @param ... Not used.
#' @export
print.kits_vertical <- function(x, digits = 2, ...) {
  need_digits(digits)
  cat(sprintf("<kits_vertical> %d cohort x wave group(s); scale: %s\n", nrow(x), scale_label(scale_of(x))))
  print(format_table(x, digits), row.names = FALSE, right = TRUE)
  co <- attr(x, "cohorts")
  if (!is.null(co)) {
    cat("Cohorts at the same wave (cohort_2 minus cohort_1):\n")
    print(format_table(co, digits), row.names = FALSE, right = TRUE)
  }
  checks_line(attr(x, "diagnostics"))
  invisible(x)
}


#' @rdname kits_vertical
#' @param type What `plot()` draws: `"mean"` (the mean +/- SD) or `"median"` (the median, bars from
#'   P25 to P75).
#' @param cohorts Which cohorts to draw, in this order (the first blue circles, the second orange
#'   squares);
#'   by default all of them.
#' @param age_ticks Age ticks.
#' @param y_range The range of the y axis, the same in every panel; by default from the data.
#' @export
plot.kits_vertical <- function(x, type = c("mean", "median"), cohorts = NULL, age_ticks = AGE_TICKS,
                               y_range = NULL, ...) {
  type <- match_option(type, c("mean", "median"), "type")
  d <- vertical_rows(x)
  if (is.null(cohorts)) cohorts <- unique(d$cohort)
  abort_if(!any(cohorts %in% d$cohort), "None of `cohorts` (", paste(cohorts, collapse = ", "), ") is in the table: ",
           paste(unique(d$cohort), collapse = ", "), ".")
  d <- filter(d, cohort %in% cohorts)
  use_age <- all(is.finite(d$age_med))
  d$x <- if (use_age) d$age_med else d$wave_m
  domains <- unique(d$domain)
  titled <- length(domains) > 1
  first <- seq_along(domains) == 1   # 圖例只在第一格
  panels <- map2(domains, first, \(sc, f) {
    vertical_panel(filter(d, domain %in% sc), sc, cohorts, type, titled, f)
  })
  y_axis <- list(title = vertical_y_title(type, scale_of(x)), range = y_range)
  title <- paste0(if (type == "mean") LAB$title_vertical else LAB$title_vertical_median,
                  age_note(d))
  x_axis <- list(title = if (use_age) LAB$age else LAB$wave_axis, tickvals = age_ticks)
  if (!titled) {
    p <- plotly::layout(panels[[1]], xaxis = x_axis, yaxis = y_axis)
    return(kits_style(p, title))
  }
  panels <- map(panels, \(p) plotly::layout(p, yaxis = y_axis))
  rows <- ceiling(length(panels) / 2)
  plotly::subplot(panels, nrows = rows, shareX = TRUE, titleY = TRUE, margin = 0.05) |>
    plotly::layout(xaxis = x_axis, xaxis2 = x_axis) |>
    kits_style(title)
}

# 實際月齡的中位數和波次差超過 1 個月的波次：標題下的註記（沒有時空字串）
age_note <- function(d) {
  late <- d |>
    filter(is.finite(age_med), abs(age_med - wave_m) > 1) |>
    summarise(age_med = stats::median(age_med), .by = wave_m) |>
    arrange(wave_m)
  if (nrow(late) == 0) {
    return("")
  }
  note <- paste0(late$wave_m, LAB$wave, " ", sprintf("%.1f", late$age_med), LAB$month,
                 collapse = " / ")
  paste0("<br>", LAB$median_age, ": ", note)
}

# 圖用的列：沒有 domain 欄（一個領域）時 NA，沒有 cohort 欄時 "all"，沒有 age_med 時 NA
vertical_rows <- function(x) {
  d <- as.data.frame(x)
  if (!"domain" %in% names(d)) d$domain <- NA_character_
  if (!"cohort" %in% names(d)) d$cohort <- "all"
  if (!"age_med" %in% names(d)) d$age_med <- NA_real_
  d
}

# 一個領域的一格：每個世代一條線；圖例只在第一格
vertical_panel <- function(d, sc, cohorts, type, titled, legend) {
  p <- plotly::plot_ly()
  for (i in seq_along(cohorts)) {
    g <- filter(d, cohort == cohorts[i])
    if (nrow(g) > 0) p <- add_vertical_trace(p, g, cohorts[i], i, type, legend)
  }
  if (titled) p <- plotly::layout(p, annotations = panel_title(lab_name("domain", sc)))
  p
}

# 一個世代的線：平均和對稱的 SD，或中位數和不對稱的 P25、P75
add_vertical_trace <- function(p, g, cohort, i, type, legend) {
  col <- COL_MANY[(i - 1) %% length(COL_MANY) + 1]
  age <- ifelse(is.finite(g$age_med), sprintf(" (%s %.1f)", LAB$median_age, g$age_med), "")
  wave <- paste0(cohort, ", ", format(g$wave_m), LAB$wave, age)
  if (type == "mean") {
    y <- g$mean
    bars <- list(type = "data", array = g$sd, color = col)
    hover <- sprintf("%s: mean %.2f, SD %.2f, n = %d", wave, g$mean, g$sd, as.integer(g$n))
  } else {
    y <- g$median
    bars <- list(type = "data", symmetric = FALSE, array = g$P75 - g$median,
                 arrayminus = g$median - g$P25, color = col)
    hover <- sprintf("%s: median %.2f (P25 %.2f, P75 %.2f), n = %d", wave, g$median, g$P25, g$P75,
                     as.integer(g$n))
  }
  marker <- list(color = col, size = MARKER_SIZE, symbol = SYM_TWO[(i - 1) %% 2 + 1])
  line <- list(color = col, width = LINE_W$main)
  plotly::add_trace(p, x = g$x, y = y, type = "scatter", mode = "lines+markers", name = cohort,
                    legendgroup = cohort, showlegend = legend, line = line, marker = marker,
                    error_y = bars, hovertext = hover, hoverinfo = "text")
}

# y 軸的標題：統計量和量尺
vertical_y_title <- function(type, scale) {
  if (type == "median") {
    return(if (is.na(scale)) LAB$median_iqr_plain else sprintf(LAB$median_iqr, scale))
  }
  if (is.na(scale)) LAB$mean_sd_plain else sprintf(LAB$mean_sd, scale)
}
