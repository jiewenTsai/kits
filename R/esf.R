# CML 的基本對稱函數（elementary symmetric functions，gamma_r）和它們給的條件分布。
# kits_esf() 是教學用的外觀；algo_esf()（algo.R）是核心，條件適配（fit.R）、題目–其餘分數
# （restscore.R）都用它。不是第二份實作：gamma 直接用 CML 的 cml_book_design() 和
# cml_gamma_matrix()（cml.R）算，同一個多項式卷積、同一個置中（減掉這個題目組合的步驟平均）。
#   gamma_r = [z^r] prod_i (1 + e_i1 z + e_i2 z^2 + ...)，e_ix = exp(-sum_{g <= x} beta_ig)
#   gamma^(-i)_r：拿掉第 i 題的乘積
#   P(X_i = x | r) = e_ix gamma^(-i)_{r - x} / gamma_r（theta 消掉了）
#   E(X_i | r) = sum_x x P，Var(X_i | r) = sum_x (x - E)^2 P

# 一個題目組合的 gamma、gamma^(-i)、條件分布和條件期望、變異數：algo_esf()（algo.R）

# 一題的 P(X_i = x | r) = e_ix gamma^(-i)_{r - x} / gamma_r，用置中的值算（c 在分子分母消掉）：
# e~_ix = exp(-sum_{g <= x} (beta_ig - c))，和 times_item_poly() 的係數相同
esf_item_prob <- function(beta_i, centre, g_rest, g_all, m_i, top) {
  e <- c(1, exp(-cumsum(beta_i - centre)))
  r <- 0:top
  P <- matrix(0, top + 1, m_i + 1)
  for (x in 0:m_i) {
    k <- r - x   # 其他題的原始分
    ok <- k >= 0 & k <= top - m_i
    P[ok, x + 1] <- e[x + 1] * g_rest[k[ok] + 1] / g_all[ok]
  }
  P
}


#' Elementary symmetric functions of CML and the conditional item-score distributions
#'
#' A teaching view of the core of conditional maximum likelihood (CML): for one set of items (one
#' item pattern), the elementary symmetric functions gamma_r and the distribution of each item
#' score given the raw score r that they imply, with its mean and variance. The numbers come from
#' the same internal code as [kits_item()], [kits_itemfit()] and [kits_item_restscore()].
#'
#' Details: the models tutorial (tutorial_models/kits_models.html), sections 2 and 6.
#'
#' @param item The item parameters: a [kits_item()] result (or its `fit`, a [kits()], or a
#'   data.frame with `item_id`, `item_score`, `beta`).
#' @param items Optional item ids: the items of the pattern (default: all items). For `plot()`:
#'   which items to draw (default: all items of `x`).
#' @param pattern_id Optional: one item pattern of `item$design` (a [kits_item()] or [kits_cml()]
#'   result); its items. Give `items` or `pattern_id`, not both.
#' @param ... Not used: an argument the function does not have stops with an error (class
#'   `kits_error_input`).
#' @return A `kits_esf` object, a list with
#'   * `gamma`: one row per raw score: `r`, `gamma`, `log_gamma` (also `as.data.frame(x)`);
#'   * `prob`: the conditional distributions, long: `item_id`, `r`, `x`, `prob`
#'     (P(X_i = x | r); 0 where x cannot occur given r);
#'   * `moments`: `item_id`, `r`, `expected` (E(X_i | r)), `variance` (Var(X_i | r));
#'   * `items`, `max_score` (the highest raw score), `steps` (the betas used), `scale`, `info`.
#' @seealso [kits_itemfit()] (conditional infit and outfit), [kits_item_restscore()],
#'   [kits_algorithm()] (the topic "cml").
#' @examples
#' # Two items by hand: a1 with 2 steps, a2 with 1 step
#' parms <- data.frame(item_id = c("a1", "a1", "a2"), item_score = c(1, 2, 1),
#'                     beta = c(-1, 0.5, 0.3))
#' esf <- kits_esf(parms)
#' esf
#' esf$prob
#'
#' # The language calibration of kitm03; bk03 is the form of 12 to 24 months
#' item <- kitm03$responses |>
#'   dplyr::filter(domain == "lan") |>
#'   kits_resp() |>
#'   kits_item(scale = "lan M03")
#' esf <- kits_esf(item, pattern_id = "bk03")
#' esf
#' p <- plot(esf)
#' if (interactive()) p
#' @export
kits_esf <- function(item, items = NULL, pattern_id = NULL, ...) {
  no_extra_args(..., fun = "kits_esf")
  call <- call_text(sys.call())
  sel <- function_steps(item, items, pattern_id)
  its <- sel$items
  esf <- algo_esf(sel$beta_of)
  structure(list(gamma = esf_gamma_table(esf), prob = esf_prob_table(esf),
                 moments = esf_moment_table(esf), items = its, max_score = esf$max_score,
                 steps = kits_table(sel$steps), scale = beta_scale_of(item),
                 info = make_info(call, list(items = its, pattern_id = pattern_id),
                                  item_data = source_of(item))),
            class = "kits_esf")
}

# 三張表：gamma；條件分布（長）；條件期望、變異數
esf_gamma_table <- function(esf) {
  kits_table(data.frame(r = 0:esf$max_score, gamma = exp(esf$log_gamma), log_gamma = esf$log_gamma))
}

esf_prob_table <- function(esf) {
  d <- prob_long(esf$prob[esf$items], esf$items, 0:esf$max_score, "r")
  d$x <- as.numeric(d$x)   # 這張表的 x 一直是 double（kits_pcm_prob() 的是 integer）
  kits_table(d)
}

esf_moment_table <- function(esf) {
  n_r <- esf$max_score + 1
  kits_table(data.frame(item_id = rep(esf$items, each = n_r),
                        r = rep(0:esf$max_score, length(esf$items)),
                        expected = as.vector(esf$expected), variance = as.vector(esf$variance)))
}

#' @rdname kits_esf
#' @param x A `kits_esf` object.
#' @export
print.kits_esf <- function(x, ...) {
  items_header(sprintf("<kits_esf> Elementary symmetric functions of CML: %d item(s), raw score 0..%d; ",
                       length(x$items), x$max_score), x$items, x$scale)
  cat("gamma_r = sum over the score vectors with raw score r of exp(-sum of the reached steps):\n")
  g <- x$gamma
  tab <- data.frame(r = g$r, gamma = formatC(g$gamma, digits = 4, format = "g"),
                    log_gamma = fmt(g$log_gamma, 3))
  if (nrow(tab) > 12) {
    gap <- data.frame(r = "...", gamma = "", log_gamma = "")
    tab <- rbind(utils::head(tab, 6), gap, utils::tail(tab, 3))
  }
  print(tab, row.names = FALSE, right = TRUE)
  cat("P(response pattern | raw score r) = exp(-sum of reached steps) / gamma_r:\n",
      "  theta cancels, so CML needs no ability distribution.\n", sep = "")
  cat("Conditional infit/outfit (kits_itemfit()) and kits_item_restscore() are built on\n",
      "  x$moments$expected = E(X_i | r) and x$moments$variance = Var(X_i | r),\n",
      "  from P(X_i = x | r) = e_ix gamma^(-i)_{r-x} / gamma_r (x$prob).\n", sep = "")
  cat("Next: plot(x) (conditional expected scores), x$prob, x$moments,\n",
      "  kits_itemfit(data, item), kits_item_restscore(data, item).\n", sep = "")
  invisible(x)
}

#' @rdname kits_esf
#' @export
as.data.frame.kits_esf <- function(x, ...) {
  drop_table_class(x$gamma)
}

#' @rdname kits_esf
#' @export
plot.kits_esf <- function(x, items = NULL, ...) {
  shown <- pick_items(x$items, items)
  d <- x$moments
  cols <- c(PAL$blue, PAL$orange, PAL$gray, "#009E73", "#CC79A7", "#E69F00", "#56B4E9", PAL$black)
  dashes <- c("solid", "dash", "dot")
  p <- plotly::plot_ly()
  for (j in seq_along(shown)) {
    e <- d[d$item_id == shown[j], ]
    hover <- sprintf("%s: r = %d, E(X | r) = %.3f, Var(X | r) = %.3f", shown[j], as.integer(e$r),
                     e$expected, e$variance)
    look <- list(color = cols[(j - 1) %% length(cols) + 1], width = LINE_W$main,
                 dash = dashes[((j - 1) %/% length(cols)) %% length(dashes) + 1])
    p <- plotly::add_trace(p, x = e$r, y = e$expected, type = "scatter", mode = "lines+markers",
                           name = shown[j], line = look,
                           marker = list(color = look$color, size = MARKER_SIZE - 2),
                           hovertext = hover, hoverinfo = "text")
  }
  p |>
    plotly::layout(xaxis = list(title = "Raw score r on these items"),
                   yaxis = list(title = "E(X_i | r) (conditional expected item score)")) |>
    kits_style("Conditional expected item scores given the raw score (CML)")
}
