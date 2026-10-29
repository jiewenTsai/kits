# 直線代數（內部用）：一條線 l = list(intercept, slope)，新值 = intercept + slope x 舊值。
# kits 裡所有的連結、換算、篩點的殘差、連結的檢查都用這幾個函數，不再各寫一份。
#   line_apply(l, x)          套用：intercept + slope x
#   line_resid(l, x, y)       殘差：y - intercept - slope x（篩點、kits_drift() 的殘差）
#   line_compose(l1, l2)      先 l1 再 l2：slope = slope_2 slope_1，
#                             intercept = intercept_2 + slope_2 intercept_1
#   line_minus(l1, l2)        兩條線的差也是一條線：intercept_1 - intercept_2、slope_1 - slope_2
#   line_stretch(l, d)        差、SE、SD 只乘 slope（截距抵消）
#   mean_sigma_moments(...)   mean/sigma 的線：slope = s_to / s_from，intercept = m_to - slope m_from
# 連結物件（kits_link）也有 intercept、slope，可以直接給這些函數。反函數和一串線的合成在 algo_line()（algo.R）。

line_apply <- function(l, x) {
  l$intercept + l$slope * x
}

line_resid <- function(l, x, y) {
  y - l$intercept - l$slope * x
}

line_compose <- function(l1, l2) {
  list(intercept = line_apply(l2, l1$intercept), slope = l2$slope * l1$slope)
}

# 兩條線在同一個 x 的差（l1 - l2）
line_minus <- function(l1, l2) {
  list(intercept = l1$intercept - l2$intercept, slope = l1$slope - l2$slope)
}

# 換算後的差、SE、SD：只乘 slope
line_stretch <- function(l, d) {
  l$slope * d
}

# mean/sigma（Kolen & Brennan 的 B、A），從兩邊的平均和 SD
mean_sigma_moments <- function(m_from, s_from, m_to, s_to) {
  slope <- s_to / s_from
  list(intercept = m_to - slope * m_from, slope = slope)
}


# ---- 換算過的分數：一串連結合成的線（summary() 用）-----------------------------------------------------

# kits_subj() 的分數記下每次 kits_convert()（attr(, "raw_theta") 的 steps）。每個校準（來源）一份：
#   chain：每一步 step、from、to、intercept、slope（照套用的順序；from 是實際經過的量尺）
#   intercept、slope：合成的線（line_compose()），from_scale、to_scale
#   theta_mean：這些列的 theta 平均換回起點量尺的值；value_mean：它換算後的值（現在的量尺）
#   se_link：在 theta_mean 的連結誤差，link_moments()，和 kits_convert() 的 se_link 相同；
#            有一條連結沒有重抽時 NA
# 沒有換算的紀錄時 NULL。疊起來的表有好幾個換過尺的來源時，一個來源一份（list）。
conversion_summary <- function(x) {
  sources <- attr(x, "raw_theta")
  if (is.null(sources) || !is.data.frame(x) || !"theta" %in% names(x)) {
    return(NULL)
  }
  out <- compact(map(sources, \(src) source_conversion(x, src)))
  if (length(out) == 0) {
    return(NULL)
  }
  if (length(out) == 1) out[[1]] else out
}

source_conversion <- function(x, src) {
  links <- unlist(map(src$steps, "links"), recursive = FALSE)
  if (length(links) == 0) {
    return(NULL)
  }
  rows <- if ("visit_id" %in% names(x)) as.character(x$visit_id) %in% src$persons else TRUE
  theta <- x$theta[rows & is.finite(x$theta)]
  to <- map_chr(links, \(l) scale_label(l$to_scale))
  from <- c(scale_label(src$scale), utils::head(to, -1))   # 實際經過的量尺
  chain <- data.frame(step = seq_along(links), from = from,
                      to = to, intercept = map_dbl(links, "intercept"),
                      slope = map_dbl(links, "slope"))
  line <- algo_line(links)
  theta_mean <- if (length(theta) > 0) algo_line(line, mean(theta), invert = TRUE) else NA_real_
  draws <- all(map_lgl(links, link_has_draws))
  se <- if (draws && is.finite(theta_mean)) sqrt(link_moments(links, theta_mean)$var) else NA_real_
  list(chain = kits_table(chain), intercept = line$intercept, slope = line$slope,
       from_scale = src$scale, to_scale = to[length(to)], theta_mean = theta_mean,
       value_mean = line_chain_apply(links, theta_mean), se_link = se)
}

# 一步一步套用（和 kits_convert() 一樣）
line_chain_apply <- function(links, x) {
  for (l in links) x <- line_apply(l, x)
  x
}

# summary() 印的換算：每一步、合成的線、連結誤差
print_conversion <- function(conv, digits) {
  if (is.null(conv$chain)) {
    for (one in conv) print_conversion(one, digits)
    return(invisible())
  }
  cat(sprintf("Conversion %s -> %s (kits_convert(), the links in the order applied):\n",
              scale_label(conv$from_scale), scale_label(conv$to_scale)))
  print(format_table(conv$chain, digits), row.names = FALSE, right = TRUE)
  cat(sprintf("  composed: new value = %s + %s x old value\n", fmt(conv$intercept, digits),
              fmt(conv$slope, digits)))
  cat("  Chained links compose as lines: slope = slope_2 x slope_1,\n",
      "  intercept = intercept_2 + slope_2 x intercept_1.\n", sep = "")
  if (is.na(conv$se_link)) {
    cat("  linking error: not available (a link has no bootstrap draws, n_boot = 0)\n")
    return(invisible())
  }
  cat(sprintf(paste0("  linking error at the mean theta (%s on %s, %s on %s): se_link %s\n",
                     "  (it accumulates through the chain: link_moments(), as kits_convert())\n"),
              fmt(conv$theta_mean, digits), scale_label(conv$from_scale),
              fmt(conv$value_mean, digits), scale_label(conv$to_scale),
              fmt(conv$se_link, digits + 1)))
  invisible()
}
