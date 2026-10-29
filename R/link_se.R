# 連結誤差在每個 theta 的大小，以及一串連結的連結誤差（Kolen & Brennan, 2014, pp. 248-258, 269-270）。
# 連結用 n_boot > 0 做時，留下每次重抽的截距和斜率（link$draws 的 intercept、slope 欄）。
# 一條線：SE_link(theta) = 重抽之間 intercept + slope theta 的 SD；
#   也就是 sqrt(Var(intercept) + theta^2 Var(slope) + 2 theta Cov(intercept, slope))。
# 一串連結：不同的連結用不同的題時，重抽彼此獨立（有共同的題時是近似）；把每條線的重抽和下一條線的每個重抽配起來
# （全部的配法），用動差遞迴算，不用亂數；中心放在估計的線上（m 是換算後的值，和 Braun-Holland 的 delta 法相同）：
#   Y_0 = x（m = x，變異 v = 0）
#   Y_k = intercept_k + slope_k Y_{k-1}：
#     v_k = Var(intercept) + 2 Cov(intercept, slope) m + Var(slope) (v + m^2) + slope_k^2 v，
#     m_k = intercept_k + slope_k m
# Var、Cov 用 n - 1 當分母（和 intercept_error、slope_error 相同），所以一條線時剛好是 sd(intercept + slope x)。
# 因為中心是換算後的值，分兩次換算（先 co，再 d_link）和一次換算的 se_link 完全相同。

# 這條連結有沒有可以用的重抽
link_has_draws <- function(link) {
  d <- link$draws
  is.matrix(d) && all(c("intercept", "slope") %in% colnames(d)) && nrow(d) >= 2
}

# 一串連結在 x 的平均和變異（x 是第一條線起點量尺上的值；v0 是 x 本身已經有的連結變異）
link_moments <- function(links, x, v0 = 0) {
  m <- x
  v <- v0 + 0 * x
  for (link in links) {
    intercept <- link$draws[, "intercept"]
    slope <- link$draws[, "slope"]
    v_int <- stats::var(intercept)
    v_slope <- stats::var(slope)
    cab <- stats::cov(intercept, slope)
    v <- v_int + 2 * cab * m + v_slope * (v + m^2) + link$slope^2 * v
    m <- line_apply(link, m)
  }
  list(mean = m, var = pmax(v, 0))
}

# kits_convert() 用：每條連結都有重抽時加上 se_link（和 se_total）；有一條沒有時不加（原本有的也拿掉）。
# earlier：前面的換算裡誤差沒有帶下來的連結（no_se_record()）；有的話也不加
add_link_se <- function(x, x0, v0, links, earlier = NULL) {
  if (!is.data.frame(x) || is.null(x0)) {
    return(x)
  }
  se_col <- if ("theta" %in% names(x)) "se" else "thr_se"
  earlier_missing <- !is.null(earlier) && nrow(earlier) > 0
  if (!all(map_lgl(links, link_has_draws)) || earlier_missing) {
    x$se_link <- NULL
    x$se_total <- NULL
    return(x)
  }
  x$se_link <- sqrt(link_moments(links, x0, v0)$var)
  if (se_col %in% names(x)) {
    x$se_total <- sqrt(x[[se_col]]^2 + x$se_link^2)
  } else {
    x$se_total <- NULL
  }
  x
}
