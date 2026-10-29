# PCM 的單一核心（內部用）：一題的類別機率 P(X_i = k | theta) 和它的動差只在這裡算。
#
# eta_k = sum_{g <= k} (theta - beta_g)，k = 0..m（Masters, 1982 的寫法；eta_0 = 0）；
# 每列先減掉最大值再取 exp（數值穩定），S = 列和，P = exp(eta) / S，log S = max + log(列和)。
# eta 每列是 cumsum(theta - beta) 的值（0.25.2 起用 rowSums() 逐欄算，見 pcm_kernel()）：
# 和 0.23.0 以前的門檻（step_threshold()）、kits_sim() 的算法一個位元都不差，
# 所以校準、連結、常模的數字不動。
#
# 用到它的地方：
#   pcm_kernel()       ：p_at_least()（calibrate.R：門檻和它的 SE）、kits_sim() 的抽樣（simulate.R）
#   pcm_probs()        ：algo_lw()（algo.R）、類別機率圖
#   algo_pcm()         ：WLE 的加權概似（ability.R）、item_moments()（fit.R：個人適配、
#                        kits_itemfit(type = "wle")）、tcc()（stats.R：kits_information()、
#                        kits_expected_score()）、每題期望分數圖（plots_results.R）、
#                        kits_pcm_prob()、latent_normal()（link_check.R）
#
# theta 或 beta 是 RTMB 的 advector 時（WLE、門檻、同齡線的自動微分），走 pcm_kernel_ad()：同樣的 eta，
# 但不比大小（RTMB 不能對「比大小」自動微分）。數值的輸入走原來的路，數字一個位元都不變。

# 一題的核心：exp(eta - 列最大值)（ex，未標準化）、列和 S、log S。列是 theta，欄是分數 0..m
pcm_kernel <- function(theta, beta_i) {
  if (is_ad(theta) || is_ad(beta_i)) {
    return(pcm_kernel_ad(theta, beta_i))
  }
  n <- length(theta)
  m <- length(beta_i)
  d <- outer(theta, beta_i, "-")
  # eta_k = sum_{g <= k} (theta - beta_g)：第 k 欄是前 k 欄的 rowSums()。rowSums() 和 cumsum() 一樣
  # 從 0 起、照順序用 long double 累加，所以和每列 cumsum() 一個位元都不差，但不用逐列 apply()
  steps <- vapply(seq_len(m), \(k) rowSums(d[, seq_len(k), drop = FALSE]), numeric(n))
  eta <- cbind(0, matrix(steps, nrow = n, ncol = m))
  top <- do.call(pmax, map(seq_len(m + 1), \(k) eta[, k]))   # 每列的最大值（逐欄比）
  ex <- exp(eta - top)
  S <- rowSums(ex)
  list(ex = ex, S = S, log_S = top + log(S))
}

# RTMB 版：eta_k = eta_(k-1) + (theta - beta_k)（和 cumsum 同一個加法順序）；
# log S 用 logspace_add() 逐欄累加（不用列最大值），ex = exp(eta - log S)，所以 S = 1
pcm_kernel_ad <- function(theta, beta_i) {
  col <- 0 * theta + 0 * beta_i[1]          # eta_0 = 0
  eta <- cbind(col)
  log_S <- col
  for (k in seq_along(beta_i)) {
    col <- col + (theta - beta_i[k])
    eta <- cbind(eta, col)
    log_S <- RTMB::logspace_add(log_S, col)
  }
  list(ex = exp(eta - log_S), S = rep(1, length(theta)), log_S = log_S)
}

# x 是不是 RTMB 的 advector（自動微分中）
is_ad <- function(x) inherits(x, "advector")

# 列和：RTMB 的 advector 不能用 rowSums()，改用矩陣乘法；數值照舊用 rowSums()
row_sum <- function(M) {
  if (is_ad(M)) drop(M %*% rep(1, ncol(M))) else rowSums(M)
}

# 類別機率 P(X_i = k | theta)：列是 theta，欄是分數 0..m
pcm_probs <- function(theta, beta_i) {
  q <- pcm_kernel(theta, beta_i)
  q$ex / q$S
}

# 類別機率和中央動差（P、E、Var、log_S，third = TRUE 時再加 M3）：algo_pcm()（algo.R）

# 測驗訊息的斜率 I'(theta)：只寫 I(theta)（tcc() 的 slope），導數由 RTMB 給。
# 每個 I 只和自己的 theta 有關，所以 sum(I) 的梯度就是每個 theta 的 I'
information_slope <- function(theta, parms, items) {
  tape <- RTMB::MakeTape(\(t) sum(tcc(t, parms, items)$slope), theta)
  drop(tape$jacobian(theta))
}

# 每題的步驟參數，照步驟排；用題號取，例如 beta_of[["a1"]]（內部用）
betas_by_item <- function(parms) {
  by_item <- parms |>
    arrange(item_score) |>
    summarise(beta = list(beta), .by = item_id)
  stats::setNames(by_item$beta, by_item$item_id)
}

# 多項式相乘，一列一個多項式（第 j 欄是 z^(j - 1) 的係數）：g(z) sum_k coef[[k + 1]] z^k，留 width 欄；
# 乘上 z^k 就是整列往右移 k 欄。CML 的 gamma（times_item_poly()；g 可以是 RTMB 的 advector）和
# algo_lw() 共用。每一項乘完馬上加：tape 的順序和以前相同，CML 的估計值一個位元都不變
poly_times <- function(g, coef, width = ncol(g) + length(coef) - 1) {
  out <- shift_cols(g, 0, width) * coef[[1]]
  for (k in seq_along(coef)[-1] - 1) {
    out <- out + shift_cols(g, k, width) * coef[[k + 1]]
  }
  out
}

# g(z) z^k：整列往右移 k 欄（左邊補 0），留 width 欄（多的截掉，不夠的補 0）
shift_cols <- function(g, k, width) {
  keep <- min(ncol(g), width - k)
  cbind(matrix(0, nrow(g), k), g[, seq_len(keep), drop = FALSE], matrix(0, nrow(g), width - k - keep))
}


#' PCM category probabilities, expected scores and information given theta
#'
#' A teaching view of the partial credit model (Masters, 1982) at given abilities: for each item and
#' theta, the category probabilities P(X_i = x | theta), the expected score E(X_i | theta) and the
#' item information Var(X_i | theta), from the same core as the rest of kits; their sums over the
#' items equal [kits_expected_score()] and [kits_information()].
#'
#' Details: the models tutorial (tutorial_models/kits_models.html), sections 1 and 9.1.
#'
#' @param item The item parameters: a [kits_item()] result (or its `fit`, a [kits()], or a
#'   data.frame with `item_id`, `item_score`, `beta`).
#' @param theta The abilities (a numeric vector), on the calibration scale of `item` (logits).
#' @param items Optional item ids (default: all items of `item`). For `plot()`: the items to draw
#'   (default: all items of `x`); `type = "category"` draws the first of them.
#' @return A `kits_pcm_prob` object, a list with
#'   * `prob`: long, one row per item, theta and score: `item_id`, `theta`, `x`, `prob`;
#'   * `moments`: one row per item and theta: `item_id`, `theta`, `expected` (E(X_i | theta)),
#'     `variance` (Var(X_i | theta), the item information), `third` (the third central moment);
#'   * `test`: one row per theta: `theta`, `expected` (the expected raw score), `information`,
#'     `information_slope` (I'(theta)); also `as.data.frame(x)`;
#'   * `items`, `steps` (the betas used), `scale`, `info`.
#' @references
#' Masters, G. N. (1982). A Rasch model for partial credit scoring. *Psychometrika, 47*(2),
#' 149-174. \doi{10.1007/BF02296272}
#' @seealso [kits_esf()] (P(X_i = x | r)), [kits_information()], [kits_expected_score()],
#'   [kits_algorithm()].
#' @examples
#' # One item with two steps, by hand
#' parms <- data.frame(item_id = "a1", item_score = 1:2, beta = c(-1, 0.5))
#' pp <- kits_pcm_prob(parms, theta = c(-2, 0, 2))
#' pp
#' pp$prob
#'
#' # Two items of the language calibration of kitm03
#' item <- kitm03$responses |>
#'   dplyr::filter(domain == "lan") |>
#'   kits_resp() |>
#'   kits_item(scale = "lan M03")
#' pp <- kits_pcm_prob(item, seq(-6, 6, by = 0.1), items = c("lana01", "lana02"))
#' p <- plot(pp, items = "lana02")
#' if (interactive()) p
#' @export
kits_pcm_prob <- function(item, theta, items = NULL) {
  call <- call_text(sys.call())
  need_theta(theta)
  theta <- as.numeric(theta)
  sel <- function_steps(item, items, NULL)
  its <- sel$items
  per_item <- map(its, \(i) algo_pcm(theta, sel$beta_of[[i]], third = TRUE))
  test <- tcc(theta, sel$steps, its)   # 和 kits_expected_score()、kits_information() 同一個函數
  structure(list(prob = kits_table(prob_long(map(per_item, "P"), its, theta, "theta")),
                 moments = pcm_moment_table(per_item, its, theta),
                 test = kits_table(data.frame(theta = theta, expected = test$tau,
                                              information = test$slope,
                                              information_slope = information_slope(theta, sel$steps, its))),
                 items = its, steps = kits_table(sel$steps),
                 scale = beta_scale_of(item),
                 info = make_info(call, list(items = its, theta = range(theta)),
                                  item_data = source_of(item))),
            class = "kits_pcm_prob")
}

# 每題的動差（長表）
pcm_moment_table <- function(per_item, its, theta) {
  map2(per_item, its, \(m, i) {
    data.frame(item_id = i, theta = theta, expected = m$E, variance = m$Var, third = m$M3)
  }) |>
    list_rbind() |>
    kits_table()
}

# ---- kits_pcm_prob() 和 kits_esf() 共用 ----

# 這些題（items 或 pattern_id 的題，或全部）的步驟，照題、分數排；beta_of：一題一個 beta 向量，照 items 的順序
function_steps <- function(item, items, pattern_id) {
  parms <- param_frame(item)
  its <- function_items(item, parms, items, pattern_id)
  steps <- arrange(filter(parms, item_id %in% its), item_id, item_score)
  list(items = its, steps = steps, beta_of = betas_by_item(steps)[its])
}

# 機率的長表：P_of 一題一個矩陣（列是 at 的值，欄是分數 0..m），變成 item_id、at（欄名 at_name）、x、prob
prob_long <- function(P_of, its, at, at_name) {
  map2(P_of, its, \(P, i) {
    d <- data.frame(item_id = i, at = rep(at, each = ncol(P)),
                    x = rep(seq_len(ncol(P)) - 1L, length(at)), prob = as.vector(t(P)))
    names(d)[2] <- at_name
    d
  }) |>
    list_rbind()
}

# print 的前兩行：標題（後面接量尺）和題目（超過 8 題時印前 6 題、...、最後 1 題）
items_header <- function(title, items, scale) {
  shown <- if (length(items) <= 8) items else c(utils::head(items, 6), "...", utils::tail(items, 1))
  cat(title, "scale ", scale_label(scale), "\n", sep = "")
  cat("Items: ", paste(shown, collapse = ", "), "\n", sep = "")
}

# as.data.frame() 給結果裡的一張表：只拿掉 kits_table 的類別
drop_table_class <- function(d) {
  class(d) <- "data.frame"
  d
}

#' @rdname kits_pcm_prob
#' @param x A `kits_pcm_prob` object.
#' @param ... Not used.
#' @export
print.kits_pcm_prob <- function(x, ...) {
  items_header(sprintf("<kits_pcm_prob> PCM category probabilities: %d item(s), %d theta value(s); ",
                       length(x$items), nrow(x$test)), x$items, x$scale)
  d <- x$test
  rows <- unique(round(seq(1, nrow(d), length.out = min(nrow(d), 9))))
  tab <- data.frame(theta = fmt(d$theta[rows], 2), expected = fmt(d$expected[rows], 3),
                    information = fmt(d$information[rows], 3),
                    slope = fmt(d$information_slope[rows], 3))
  print(tab, row.names = FALSE, right = TRUE)
  if (length(rows) < nrow(d)) cat(sprintf("(%d of %d theta values; all in x$test)\n", length(rows), nrow(d)))
  cat("P(X = x | theta) of the PCM; information = Var(X | theta), slope = I'(theta).\n",
      "  expected and information are kits_expected_score(item)(theta) and\n",
      "  kits_information(item)(theta), from the same core.\n", sep = "")
  cat("kits_esf() is the other way round: P(X = x | r), conditional on the raw score r,\n",
      "  with theta gone; here the probabilities are conditional on theta, not on r.\n", sep = "")
  cat("Next: plot(x) (category curves), plot(x, type = \"expected\"), ",
      "plot(x, type = \"information\"),\n  x$prob, x$moments, kits_esf(item).\n", sep = "")
  invisible(x)
}

#' @rdname kits_pcm_prob
#' @export
as.data.frame.kits_pcm_prob <- function(x, ...) {
  drop_table_class(x$test)
}

#' @rdname kits_pcm_prob
#' @param type `"category"` (the default): the category curves P(X = x | theta) of one item;
#'   `"expected"`: the expected test score (blue) and item scores (gray);
#'   `"information"`: the test information (blue), its slope I'(theta) (orange, dashed) and the
#'   item information (gray).
#' @export
plot.kits_pcm_prob <- function(x, type = c("category", "expected", "information"), items = NULL, ...) {
  type <- match.arg(type)
  shown <- pick_items(x$items, items)
  if (type == "category") {
    return(pcm_category_plot(x, shown[1]))
  }
  info <- type == "information"
  test_y <- if (info) x$test$information else x$test$expected
  item_y <- if (info) "variance" else "expected"
  what <- if (info) "Information" else "Expected score"
  p <- plotly::plot_ly()
  for (i in shown) {
    m <- x$moments[x$moments$item_id == i, ]
    hover <- sprintf("%s<br>theta %.2f: %s %.3f", i, m$theta, tolower(what), m[[item_y]])
    p <- plotly::add_lines(p, x = m$theta, y = m[[item_y]], name = i, legendgroup = "items",
                           showlegend = i == shown[1], line = list(color = PAL$gray, width = LINE_W$thin),
                           hovertext = hover, hoverinfo = "text")
  }
  hover <- sprintf("Test<br>theta %.2f: %s %.3f", x$test$theta, tolower(what), test_y)
  p <- plotly::add_lines(p, x = x$test$theta, y = test_y, name = "Test", hovertext = hover,
                         hoverinfo = "text", line = list(color = PAL$blue, width = LINE_W$main))
  if (info) {
    hover <- sprintf("theta %.2f: I'(theta) %.3f", x$test$theta, x$test$information_slope)
    p <- plotly::add_lines(p, x = x$test$theta, y = x$test$information_slope, name = "I'(theta)",
                           hovertext = hover, hoverinfo = "text",
                           line = list(color = PAL$orange, width = LINE_W$main, dash = "dash"))
  }
  title <- if (info) "Information = Var(X | theta) under the PCM" else "Expected scores under the PCM"
  plotly::layout(p, xaxis = list(title = "theta"), yaxis = list(title = what)) |>
    kits_style(title)
}

# 一題的類別機率曲線：x$prob 的這一題變回矩陣（列是 theta，欄是分數），每個分數一條線（category_lines()）
pcm_category_plot <- function(x, item) {
  d <- x$prob[x$prob$item_id == item, ]
  P <- matrix(d$prob, nrow = nrow(x$test), byrow = TRUE)
  category_lines(item, x$test$theta, P) |>
    plotly::layout(xaxis = list(title = "theta"), yaxis = list(title = "P(X = x | theta)", range = c(0, 1))) |>
    kits_style(sprintf("Category probabilities of %s under the PCM", item))
}
