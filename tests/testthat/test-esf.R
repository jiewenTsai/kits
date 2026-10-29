# kits_esf()：CML 的基本對稱函數 gamma_r 和條件分布（教學用）；algo_esf() 是條件適配、
# 題目–其餘分數共用的核心。檢查：和直接列舉相同、和 CML 的概似用的 gamma 相同、機率和期望的恆等式。

# 三題的小例子：步驟數 2、1、3
tiny <- data.frame(item_id = c("a1", "a1", "a2", "a3", "a3", "a3"),
                   item_score = c(1, 2, 1, 1, 2, 3),
                   beta = c(-1.2, 0.4, 0.3, -0.5, 0.8, 1.6))

# 直接列舉：所有分數組合，權重 exp(-sum 達到的步驟)，按原始分加起來（和方程式教學的檢查相同）
enumerate_tiny <- function(parms) {
  beta_of <- split(parms$beta, parms$item_id)
  cum_of <- lapply(beta_of, \(b) c(0, cumsum(b)))
  grid <- expand.grid(lapply(beta_of, \(b) 0:length(b)))
  cum <- Reduce(`+`, Map(\(k, cb) cb[k + 1], grid, cum_of))
  list(grid = grid, r = rowSums(grid), w = exp(-cum))
}

test_that("gamma_r 和直接列舉所有分數組合相同", {
  esf <- kits_esf(tiny)
  en <- enumerate_tiny(tiny)
  gamma_enum <- as.vector(tapply(en$w, en$r, sum))
  expect_s3_class(esf, "kits_esf")
  expect_identical(esf$items, c("a1", "a2", "a3"))
  expect_identical(esf$max_score, 6L)
  expect_equal(esf$gamma$gamma, gamma_enum, tolerance = 1e-12)
  expect_equal(esf$gamma$log_gamma, log(gamma_enum), tolerance = 1e-12)
  expect_identical(esf$gamma$gamma[1], 1)   # gamma_0 = 1：只有全部 0 分
  # 兩題的子集合
  two <- kits_esf(tiny, items = c("a1", "a2"))
  en2 <- enumerate_tiny(tiny[tiny$item_id != "a3", ])
  expect_equal(two$gamma$gamma, as.vector(tapply(en2$w, en2$r, sum)), tolerance = 1e-12)
})

test_that("P(X_i = x | r) 和列舉相同，對 x 加起來是 1；E(X_i | r) 加起來是 r", {
  esf <- kits_esf(tiny)
  en <- enumerate_tiny(tiny)
  for (i in esf$items) {
    by_rx <- tapply(en$w, list(en$r, en$grid[[i]]), sum)
    by_rx[is.na(by_rx)] <- 0
    p_enum <- by_rx / rowSums(by_rx)   # 一列一個 r，一欄一個 x
    d <- esf$prob[esf$prob$item_id == i, ]
    expect_equal(d$prob, as.vector(t(p_enum)), tolerance = 1e-12)
  }
  sums <- tapply(esf$prob$prob, list(esf$prob$item_id, esf$prob$r), sum)
  expect_equal(as.vector(sums), rep(1, length(sums)), tolerance = 1e-12)
  e_sum <- tapply(esf$moments$expected, esf$moments$r, sum)
  expect_equal(as.vector(e_sum), 0:esf$max_score, tolerance = 1e-12)
  expect_true(all(esf$moments$variance >= -1e-15))
  # 0 分和滿分：每一題只有一個可能的分數，變異數是 0
  ends <- esf$moments[esf$moments$r %in% c(0, esf$max_score), ]
  expect_equal(ends$variance, rep(0, nrow(ends)), tolerance = 1e-12)
})

test_that("algo_esf()：gamma^(-i) 和拿掉這題後的 gamma 相同；free 是 x 有兩個以上的可能", {
  beta_of <- split(tiny$beta, tiny$item_id)
  m <- algo_esf(beta_of)
  for (i in names(beta_of)) {
    rest <- kits_esf(tiny, items = setdiff(names(beta_of), i))$gamma$log_gamma
    expect_equal(m$log_gamma_rest[i, seq_along(rest)], rest, tolerance = 1e-12)
    expect_true(all(m$log_gamma_rest[i, -seq_along(rest)] == -Inf))
  }
  n_x <- vapply(m$prob, \(P) rowSums(P > 1e-12), numeric(m$max_score + 1))
  expect_identical(m$free, n_x >= 2)
})

test_that("和 CML 的概似用的 gamma 完全相同（kitm03 12 個月的語言題本）", {
  item <- kitm03$responses |>
    dplyr::filter(domain == "lan", wave_m == 12) |>
    kits_resp() |>
    kits_item(scale = "lan M12") |>
    suppressMessages()
  pats <- item$stats$patterns
  form <- pats$pattern_id[which.max(pats$n)]   # 12 個月的標準題本：最多人的作答組合
  esf <- kits_esf(item, pattern_id = form)
  expect_identical(esf$items, sort(pats$items[[which.max(pats$n)]]))
  # CML 的概似的 log gamma_br（cml_nll() 的那幾行），在估計值
  inp <- cml_inputs(item$stats)
  design <- cml_book_design(inp)
  beta <- item$fit$est$beta
  centre <- drop(design$has_step %*% beta) / design$n_steps
  G <- cml_gamma_matrix(beta, centre, inp, design)
  log_gamma <- log(G[design$cells]) - design$raw * centre[design$book]
  b <- match(form, pats$pattern_id)
  expect_equal(esf$gamma$log_gamma, log_gamma[design$book == b], tolerance = 1e-12)
  # 用 kits_esf() 的 gamma 重新組出 -log L：和校準的最小值相同
  nll <- sum(beta * inp$n_ge) + sum(vapply(seq_along(inp$book_items), \(k) {
    sum(inp$n_raw[[k]] * kits_esf(item, pattern_id = pats$pattern_id[k])$gamma$log_gamma)
  }, numeric(1)))
  expect_equal(nll, item$fit$est$nll, tolerance = 1e-12)
  # 恆等式
  expect_equal(as.vector(tapply(esf$moments$expected, esf$moments$r, sum)), 0:esf$max_score,
               tolerance = 1e-10)
  # print、plot、as.data.frame
  out <- capture.output(print(esf))
  expect_true(any(grepl("theta cancels, so CML needs no ability distribution", out)))
  expect_true(any(grepl("kits_item_restscore()", out, fixed = TRUE)))
  expect_true(any(grepl("^Next:", out)))
  expect_true(any(grepl("...", out, fixed = TRUE)))   # 長的表只印頭尾
  expect_identical(as.data.frame(esf), structure(esf$gamma, class = "data.frame"))
  p <- plot(esf, items = esf$items[1:3])
  expect_s3_class(p, "plotly")
  expect_length(p$x$data, 3)
  expect_identical(p$x$data[[1]]$line$color, "#0072B2")
  expect_identical(p$x$data[[2]]$line$color, "#D55E00")
  expect_identical(p$x$data[[3]]$line$color, "#7F7F7F")
  expect_s3_class(kits_info(esf), "kits_info")
})

test_that("長的題本（80 題）不會下溢：機率和期望的恆等式還是成立", {
  set.seed(1)
  n <- 80
  m <- rep(1:3, length.out = n)
  ids <- sprintf("i%02d", seq_len(n))
  long <- data.frame(item_id = rep(ids, m), item_score = sequence(m),
                     beta = stats::rnorm(sum(m), sd = 2) + rep(seq(-4, 4, length.out = n), m))
  esf <- kits_esf(long)
  expect_true(all(is.finite(esf$gamma$log_gamma)))
  expect_true(all(is.finite(esf$prob$prob)))
  sums <- tapply(esf$prob$prob, list(esf$prob$item_id, esf$prob$r), sum)
  expect_equal(as.vector(sums), rep(1, length(sums)), tolerance = 1e-10)
  expect_equal(as.vector(tapply(esf$moments$expected, esf$moments$r, sum)), 0:esf$max_score,
               tolerance = 1e-8)
})

test_that("輸入錯誤時停下來", {
  expect_error(kits_esf(tiny, items = "a1", pattern_id = "bk1"), class = "kits_error_input")
  expect_error(kits_esf(tiny, items = "zz"), class = "kits_error_items")
  expect_error(kits_esf(tiny, pattern_id = "bk1"), class = "kits_error_input")
  expect_error(kits_esf(tiny, foo = 1), class = "kits_error_input")
  expect_error(plot(kits_esf(tiny), items = "zz"), class = "kits_error_input")
})
