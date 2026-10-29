# 題目–其餘分數的關聯（Kreiner, 2011）：觀察到的 Goodman–Kruskal gamma 對比 PCM 下的期望 gamma。
# 和 iarm::item_restscore() 的算法相同（完整資料時數字相同）：
#   觀察值：題目分數 x 其餘分數（總分減這題）的交叉表的 gamma，
#     SE 是它的漸近標準誤（vcdExtra::GKgamma 的 sigma）。
#   期望值：同一張表在模型下的期望次數的 gamma。總分 t 的人數 N_t 照觀察的，
#     給定總分時這題的分數照 CML 的
#     條件分布 P(X_i = x | T = t) = e_ix gamma^(-i)_{t-x} / gamma_t（不需要 theta）。
#     期望次數 E[r, x] = N_{r+x} P(X_i = x | T = r + x)。
#   z = (觀察 - 期望) / SE，雙尾常態 p 值，所有列一起做 BH 校正。
# 不完整設計：每個題目組合（作答型態）各算一次。其餘分數是同一個組合裡其他題的總分，
# 題目和滿分都不同，所以不同組合的其餘分數不能放在一起比；條件分布也只看這個組合的題。
# 人數少於 min_n 的組合跳過（gamma 的 SE 是漸近的）。

#' Item-restscore association: observed against expected gamma
#'
#' Kreiner's (2011) test of item fit: for each item, the Goodman-Kruskal gamma of item score x rest
#' score (the sum of the other items) is compared with the gamma the PCM implies given the item
#' parameters. Each row is one item in one item pattern with at least `min_n` persons.
#' `gamma_obs` below `gamma_exp` means the item separates persons less than the model implies,
#' above means more (for example local dependence); with large N, look at the size of the
#' difference.
#'
#' Details: the models tutorial (tutorial_models/kits_models.html), sections 6 and 9.2.
#'
#' @inheritParams kits_itemfit
#' @param min_n Smallest number of persons in an item pattern for it to be tested (skipped
#'   patterns get a note in the diagnostics).
#' @param ... Not used: an argument the function does not have stops with an error (class
#'   `kits_error_input`).
#' @return A data.frame (class `kits_item_restscore`) with one row per item and item pattern:
#'   * `item_id`, `pattern_id` (as in [kits_subj()]), `n` (persons);
#'   * `gamma_obs`, `gamma_exp`: the observed gamma and the gamma expected under the PCM;
#'   * `se` (asymptotic SE of `gamma_obs`; Goodman & Kruskal, 1972), `z`, `p` (two-sided, normal);
#'   * `p_adj`: Benjamini-Hochberg (1995) over all rows.
#'
#'   NA when the gamma cannot be computed (the item or the rest score is constant in the pattern).
#'   The skipped patterns are in `attr(x, "skipped")`.
#' @references
#' Kreiner, S. (2011). A note on item-restscore association in Rasch models.
#' *Applied Psychological Measurement, 35*, 557-561.
#'
#' Benjamini, Y., & Hochberg, Y. (1995). Controlling the false discovery rate: a practical and
#' powerful approach to multiple testing. *Journal of the Royal Statistical Society B, 57*,
#' 289-300.
#'
#' Goodman, L. A., & Kruskal, W. H. (1972). Measures of association for cross classifications,
#' IV: Simplification of asymptotic variances. *Journal of the American Statistical Association,
#' 67*, 415-421.
#' @seealso [kits_itemfit()] for infit, outfit and the item-rest correlation.
#' @examples
#' betas <- list(a1 = c(-1, 0), a2 = c(0, 1), a3 = 0.5, a4 = -0.5, a5 = 0)
#' X <- kits_sim(seq(-3, 3, length.out = 500), betas, seed = 1)
#' X[, "a5"] <- rep(0:1, 250)   # a5 is unrelated to the others
#' item <- kits_item(X)
#' rs <- kits_item_restscore(X, item)
#' rs
#' @export
kits_item_restscore <- function(data, item, min_n = 100, ...) {
  no_extra_args(..., fun = "kits_item_restscore")
  call <- call_text(sys.call())
  need_min_n(min_n)
  resp <- resp_long(data, gaps = FALSE)
  parms <- param_frame(item)
  need_items(resp, parms)
  persons <- person_sets(resp, parms) |>
    add_booklet_ids()
  patterns <- persons |>
    summarise(n = n(), .by = c(pattern_id, set)) |>
    mutate(n_items = lengths(map(set, set_items))) |>
    arrange(pattern_id)
  used <- filter(patterns, n >= min_n, n_items >= 2)
  abort_if(nrow(used) == 0, sprintf("No item pattern has at least %s persons and 2 items ", format(min_n)),
           sprintf("(the largest has %d persons). Lower `min_n`.", max(patterns$n)), class = "kits_error_data")
  beta_of <- betas_by_item(parms)
  rows <- map(seq_len(nrow(used)), \(b) {
    visits <- persons$visit_id[persons$pattern_id == used$pattern_id[b]]
    restscore_pattern(resp, visits, set_items(used$set[b]), beta_of) |>
      mutate(pattern_id = used$pattern_id[b], n = used$n[b], .after = item_id)
  })
  out <- list_rbind(rows) |>
    mutate(z = (gamma_obs - gamma_exp) / se,
           p = 2 * stats::pnorm(-abs(z)),
           p_adj = stats::p.adjust(p, method = "BH")) |>
    arrange(sort_rank(item_id), pattern_id) |>
    kits_table("kits_item_restscore")
  skipped <- filter(patterns, !pattern_id %in% used$pattern_id) |>
    select(pattern_id, n, n_items)
  attr(out, "skipped") <- skipped
  attr(out, "min_n") <- min_n
  attr(out, "diagnostics") <- restscore_diagnostics(out, used, skipped, min_n)
  attr(out, "info") <- make_info(call, list(min_n = min_n, scale = beta_scale_of(item)),
                                 data_fingerprint(resp), item_data = source_of(item))
  out
}

need_min_n <- function(min_n) {
  need_number(min_n, "min_n", "one number of persons (0 or more), e.g. 100", lower = 0)
}

# 一個題目組合：每一題的觀察 gamma、它的 SE、期望 gamma（一列一題）。
# 條件分布 P(X_i = x | T = t)：algo_esf()（algo.R；和 CML 的 cml_gamma_matrix() 同一段程式）
restscore_pattern <- function(resp, visits, items, beta_of) {
  X <- score_matrix(resp, visits, items)
  esf <- algo_esf(beta_of[items])
  total <- rowSums(X)
  n_total <- tabulate(total + 1, nbins = esf$max_score + 1)   # N_t，t = 0..滿分
  cols <- map(seq_along(items), \(j) {
    obs <- item_rest_table(X[, j], total - X[, j], esf$m[j] + 1, esf$max_score - esf$m[j] + 1)
    exp <- expected_rest_table(esf$prob[[j]], n_total)
    g_obs <- gk_gamma(obs)
    data.frame(item_id = items[j], gamma_obs = g_obs$gamma, gamma_exp = gk_gamma(exp)$gamma,
               se = g_obs$se)
  })
  list_rbind(cols)
}

# 這些人次（列）x 這些題（欄）的分數矩陣
score_matrix <- function(resp, visits, items) {
  r <- resp[resp$visit_id %in% visits, ]
  X <- matrix(0L, length(visits), length(items))
  X[cbind(match(r$visit_id, visits), match(r$item_id, items))] <- as.integer(r$item_score)
  X
}

# 觀察的表：列是其餘分數 0..R，欄是題目分數 0..m（n_x、n_r 是欄數、列數）
item_rest_table <- function(x, rest, n_x, n_r) {
  counts <- tabulate(rest + 1 + n_r * x, nbins = n_r * n_x)
  matrix(counts, n_r, n_x)
}

# 期望的表：列是其餘分數 k = 0..R - m_i，欄是題目分數 x；總分 t = k + x：
# E[k, x] = N_t P(X_i = x | T = t)
expected_rest_table <- function(P, n_total) {
  rest <- 0:(nrow(P) - ncol(P))
  outer(rest, seq_len(ncol(P)) - 1, \(k, x) n_total[k + x + 1] * P[cbind(k + x + 1, x + 1)])
}

# Goodman–Kruskal gamma 和它的漸近 SE（和 vcdExtra::GKgamma() 相同）。只寫 gamma，導數由 RTMB 給：
#   n：格子的次數（一格一個）；C、D：一致、不一致的配對數，C = n' K_C n / 2、D = n' K_D n / 2
#   （K_C[a, b] = 1：格子 a、b 的列和欄同方向不同；K_D：反方向）；gamma = (C - D) / (C + D)
#   delta method（多項分布）：SE^2 = sum n J^2 - (sum n J)^2，J = d gamma / d n
gk_gamma <- function(tab) {
  n <- as.vector(tab)
  K <- concordance_signs(nrow(tab), ncol(tab))
  gamma_of <- function(n) {
    C <- sum(n * drop((K > 0) %*% n)) / 2
    D <- sum(n * drop((K < 0) %*% n)) / 2
    (C - D) / (C + D)
  }
  C_plus_D <- sum(n * drop((K != 0) %*% n)) / 2
  if (!is.finite(C_plus_D) || C_plus_D <= 0) {
    return(list(gamma = NA_real_, se = NA_real_))
  }
  J <- drop(RTMB::MakeTape(gamma_of, n)$jacobian(n))
  list(gamma = gamma_of(n), se = sqrt(max(sum(n * J^2) - sum(n * J)^2, 0)))
}

# 兩兩格子的方向（照 as.vector() 的格子順序）：+1 一致（列和欄同方向不同），-1 不一致，0 同列或同欄
concordance_signs <- function(n_row, n_col) {
  i <- rep(seq_len(n_row), n_col)
  j <- rep(seq_len(n_col), each = n_row)
  sign(outer(i, i, "-") * outer(j, j, "-"))
}

restscore_diagnostics <- function(out, used, skipped, min_n) {
  n_na <- sum(is.na(out$gamma_obs) | is.na(out$gamma_exp))
  skip_msg <- if (nrow(skipped) == 0) {
    "none"
  } else {
    sprintf("%s with fewer than %s persons or only 1 item (%d person-waves); not tested",
            plural(nrow(skipped), "item pattern"), format(min_n), sum(skipped$n))
  }
  diag_table(list(
    diag_row("item patterns tested", nrow(used), "ok",
             sprintf("%s with at least %s persons (%d person-waves)",
                     plural(nrow(used), "item pattern"), format(min_n), sum(used$n))),
    diag_row("item patterns skipped", nrow(skipped), if (nrow(skipped) > 0) "note" else "ok",
             skip_msg),
    diag_row("rows without gamma", n_na, if (n_na > 0) "note" else "ok",
             if (n_na > 0) "the item or the rest score is constant in the pattern" else "none")))
}

#' @rdname kits_item_restscore
#' @param x A `kits_item_restscore`.
#' @param digits Decimals to print.
#' @param n Rows to print in each of the two lists (significant rows; largest differences).
#' @export
print.kits_item_restscore <- function(x, digits = 3, n = 10, ...) {
  need_digits(digits)
  need_rows(n)
  d <- as.data.frame(x)
  n_pat <- length(unique(d$pattern_id))
  size <- sprintf("%s (%s x %s)", plural(nrow(d), "row"), plural(length(unique(d$item_id)), "item"),
                  plural(n_pat, "item pattern"))
  title <- "<kits_item_restscore> item-restscore gamma (Kreiner, 2011), observed vs expected"
  cat(sprintf("%s: %s\n", title, size))
  skipped <- attr(x, "skipped")
  if (!is.null(skipped) && nrow(skipped) > 0) {
    cat(sprintf("Skipped: %s with fewer than %s persons or a single item (%d person-waves).\n",
                plural(nrow(skipped), "item pattern"), format(attr(x, "min_n")), sum(skipped$n)))
  }
  sig <- which(d$p_adj < 0.05)
  cat(sprintf("p_adj < 0.05 (Benjamini-Hochberg): %d of %d rows\n", length(sig), nrow(d)))
  if (length(sig) > 0) {
    print_rows(x[sig[order(d$p_adj[sig])], ], digits, min(n, length(sig)), row.names = FALSE)
  }
  diff <- abs(d$gamma_obs - d$gamma_exp)
  top <- utils::head(order(-diff, na.last = NA), min(n, 5))
  if (length(top) > 0) {
    cat("Largest |gamma_obs - gamma_exp|:\n")
    print_rows(x[top, ], digits, length(top), row.names = FALSE)
  }
  cat("With large N, small differences become significant:",
      "look at the size of gamma_obs - gamma_exp.\n")
  cat("Next: plot(x), kits_itemfit(data, item), as.data.frame(x) for all rows.\n")
  checks_line(attr(x, "diagnostics"))
  invisible(x)
}

#' @rdname kits_item_restscore
#' @param items Which items to draw (item ids); by default all of them.
#' @examples
#'
#' # Observed against expected gamma; off the line = misfit
#' p <- plot(rs)
#' if (interactive()) p
#' @export
plot.kits_item_restscore <- function(x, items = NULL, ...) {
  d <- as.data.frame(x)
  shown <- pick_items(unique(d$item_id), items)
  d <- filter(d, item_id %in% shown, is.finite(gamma_obs), is.finite(gamma_exp))
  one_pattern <- length(unique(d$pattern_id)) <= 1
  d$label <- if (one_pattern) d$item_id else paste0(d$item_id, " (", d$pattern_id, ")")
  sig <- d$p_adj < 0.05
  d$group <- case_when(sig & d$gamma_obs < d$gamma_exp ~ "lower than expected (p_adj < 0.05)",
                       sig ~ "higher than expected (p_adj < 0.05)",
                       .default = "not significant")
  looks <- list(`lower than expected (p_adj < 0.05)` = PAL$orange,
                `higher than expected (p_adj < 0.05)` = PAL$blue, `not significant` = PAL$gray)
  ends <- range(c(d$gamma_obs - 1.96 * d$se, d$gamma_obs + 1.96 * d$se, d$gamma_exp), na.rm = TRUE)
  ends <- ends + c(-0.02, 0.02)
  p <- plotly::plot_ly() |>
    plotly::add_lines(x = ends, y = ends, name = "observed = expected", hoverinfo = "skip",
                      line = list(color = PAL$gray, width = LINE_W$thin, dash = "dash"))
  for (g in intersect(names(looks), d$group)) {
    e <- d[d$group == g, ]
    hover <- sprintf("%s: observed %.3f, expected %.3f (n = %d, p_adj = %s)", e$label, e$gamma_obs,
                     e$gamma_exp, as.integer(e$n), format.pval(e$p_adj, digits = 2))
    p <- plotly::add_markers(p, x = e$gamma_exp, y = e$gamma_obs, name = g,
                             error_y = list(type = "data", array = 1.96 * e$se, color = looks[[g]],
                                            thickness = LINE_W$thin, width = 0),
                             marker = list(color = looks[[g]], size = MARKER_SIZE + 1),
                             hovertext = hover, hoverinfo = "text")
  }
  p |>
    plotly::layout(xaxis = list(title = "Expected gamma (PCM)", range = ends),
                   yaxis = list(title = "Observed gamma (+/- 1.96 SE)", range = ends)) |>
    kits_style("Item-restscore association (Kreiner, 2011)")
}
