# 0.21.0：WLE 改成 base R 的牛頓法（封閉式的導數），拿掉 RTMB；數字和 0.20.0 相同。
# CML 留在 RTMB（R 的解析梯度比較慢，見 NEWS）；這裡用列舉驗證它的梯度公式。

wle020 <- readRDS(test_path("fixtures", "wle020.rds"))   # kits 0.20.0（RTMB）的 WLE 表
# 0.23.0 的欄名：booklet_id -> pattern_id、booklet_score -> raw、person_id -> visit_id
new_names <- function(d) {
  if (!is.data.frame(d)) return(d)
  old <- c(booklet_id = "pattern_id", booklet_score = "raw", booklet_max = "raw_max",
           person_id = "visit_id")
  hit <- names(d) %in% names(old)
  names(d)[hit] <- old[names(d)[hit]]
  d
}
wle020 <- lapply(wle020, \(case) lapply(case, new_names))

max_diff <- function(a, b) max(abs(a - b))

# f'(theta) = r - E(theta) + I'(theta) / (2 I(theta))，和 SE = 1 / sqrt(I(theta))，照式子一題一題算
wle_score <- function(theta, raw, parms) {
  moments <- lapply(split(parms$beta, parms$item_id), \(b) {
    p <- kits:::pcm_probs(theta, b)
    k <- seq_along(p) - 1
    mu <- sum(k * p)
    c(mu = mu, var = sum((k - mu)^2 * p), m3 = sum((k - mu)^3 * p))
  })
  m <- Reduce(`+`, moments)
  c(score = raw - m[["mu"]] + m[["m3"]] / (2 * m[["var"]]),
    se = 1 / sqrt(m[["var"]]))
}

# f(theta) = r theta - sum_i log S_i(theta) + 1/2 log I(theta)
wle_target <- function(theta, raw, parms) {
  parts <- lapply(split(parms$beta, parms$item_id), \(b) {
    k <- seq_len(length(b) + 1) - 1
    p <- kits:::pcm_probs(theta, b)
    c(log_S = log(sum(exp(k * theta - c(0, cumsum(b))))), var = sum(k^2 * p) - sum(k * p)^2)
  })
  s <- Reduce(`+`, parts)
  raw * theta - s[["log_S"]] + 0.5 * log(s[["var"]])
}

test_that("WLE：theta 和 0.20.0（RTMB）相同，到 1e-10：多題本、固定參數、beta 整體偏移 60、兩個高點（far）", {
  for (case in names(wle020)) {   # far：0.20.0 和 0.25.0 都取比較高的高點
    ref <- wle020[[case]]
    tab <- kits:::wle_table(ref$parms, ref$design)
    expect_equal(tab$pattern_id, ref$table$pattern_id, info = case)
    expect_equal(tab$raw, ref$table$raw, info = case)
    expect_lt(max_diff(tab$theta, ref$table$theta), 1e-10)
    expect_true(all(tab$se <= ref$table$se + 1e-12))   # 1 / sqrt(I) <= Warm 的 sqrt(I + (I'/2I)^2) / I
  }
  s <- kits_subj(wle020$poly$X, wle020$poly$parms)
  expect_equal(s$visit_id, wle020$poly$subj$visit_id)
  expect_lt(max_diff(s$theta, wle020$poly$subj$theta), 1e-10)
})

test_that("WLE 的封閉式：r - E + I'/2I 在估計值是 0，SE 是 1 / sqrt(I)；0 分和滿分也是", {
  parms <- wle020$poly$parms
  tab <- kits_subj_table(parms)   # 所有題一個題本，0 到滿分
  top <- max(tab$raw)
  for (r in tab$raw) {
    th <- tab$theta[tab$raw == r]
    at <- wle_score(th, r, parms)
    expect_lt(abs(at[["score"]]), 1e-9)
    expect_equal(tab$se[tab$raw == r], at[["se"]], tolerance = 1e-10)
  }
  expect_true(all(is.finite(tab$theta)) && all(is.finite(tab$se)))
  expect_true(all(diff(tab$theta) > 0))   # 原始分越高，theta 越高，包括 0 分和滿分
  expect_lt(tab$theta[1], tab$theta[2])
  expect_gt(tab$theta[top + 1], tab$theta[top])
})

test_that("WLE：兩個高點時取比較高的（f 的值決定），一定是高點、不是中間的低點；一樣高時取 theta 小的", {
  ref <- wle020$far
  tab <- kits:::wle_table(ref$parms, ref$design)
  for (b in unique(ref$design$pattern_id)) {
    its <- ref$design$item_id[ref$design$pattern_id == b]
    p <- ref$parms[ref$parms$item_id %in% its, ]
    rows <- tab$pattern_id == b
    grid <- seq(min(p$beta) - 3, max(p$beta) + 3, by = 0.01)
    for (r in tab$raw[rows]) {
      th <- tab$theta[rows & tab$raw == r]
      expect_lt(abs(wle_score(th, r, p)[["score"]]), 1e-8)                                  # f'(theta) = 0
      expect_gte(wle_target(th, r, p), max(map_dbl(grid, wle_target, raw = r, parms = p)))   # 最高的高點
    }
  }
  # two items with steps far apart: f has two maxima of unequal height; 0.24.1 stopped at the lower one
  # (the one Newton reaches from TAM's start value), 0.25.0 takes the higher one
  cases <- list(
    list(parms = data.frame(item_id = rep(c("p1", "p2"), c(3, 2)), item_score = c(1:3, 1:2),
                            beta = c(-4, -3, -2, 3, 2)), raw = 3, low = -1.083667, high = 1.663137),
    list(parms = data.frame(item_id = c("d1", "p2", "p2", "p2"), item_score = c(1, 1:3),
                            beta = c(0, -1, 5, 5)), raw = 2, low = 1.433598, high = 3.932702))
  for (case in cases) {
    tab <- kits_subj_table(case$parms)
    th <- tab$theta[tab$raw == case$raw]
    expect_equal(th, case$high, tolerance = 1e-6)
    expect_lt(abs(wle_score(th, case$raw, case$parms)[["score"]]), 1e-9)
    low <- stats::optimize(wle_target, case$low + c(-0.5, 0.5), raw = case$raw, parms = case$parms,
                           maximum = TRUE, tol = 1e-10)
    expect_equal(low$maximum, case$low, tolerance = 1e-5)   # the other, lower maximum
    expect_gt(wle_target(th, case$raw, case$parms) - low$objective, 0.15)
  }
  # two dichotomous items 6 logits apart, raw score 1: two maxima equally high, the lower theta is kept
  sym <- data.frame(item_id = c("d1", "d2"), item_score = 1, beta = c(-3, 3))
  th <- kits_subj_table(sym)$theta[2]
  expect_lt(th, 0)
  expect_lt(abs(wle_score(th, 1, sym)[["score"]]), 1e-9)
  expect_equal(wle_target(th, 1, sym), wle_target(-th, 1, sym), tolerance = 1e-12)
})

test_that("WLE：固定參數（含沒放進校準的固定題）的換算表和 kits_subj() 一致", {
  ref <- wle020$fixed
  tab <- kits:::wle_table(ref$parms, ref$design)
  expect_true("a4" %in% ref$design$item_id)
  scores <- vapply(seq_len(nrow(tab)), \(j) {
    its <- ref$design$item_id[ref$design$pattern_id == tab$pattern_id[j]]
    p <- ref$parms[ref$parms$item_id %in% its, ]
    wle_score(tab$theta[j], tab$raw[j], p)[["score"]]
  }, numeric(1))
  expect_lt(max(abs(scores)), 1e-9)
})

test_that("WLE：theta 和 dexter 的 ability_tables(method = \"WLE\") 相同，30 個隨機的題目組合", {
  skip_if_not_installed("dexter")
  # dexter sets a whole booklet to NA when one of its WLEs does not converge; those designs are not compared
  diffs <- vapply(1:30, \(seed) {
    set.seed(seed)
    n_items <- sample(2:7, 1)
    spread <- if (seed %% 3 == 0) 4 else 1.5
    steps <- sample(1:3, n_items, replace = TRUE)
    parms <- data.frame(item_id = rep(sprintf("i%d", seq_len(n_items)), steps), item_score = sequence(steps),
                        beta = rep(rnorm(n_items, 0, spread), steps) + rnorm(sum(steps), 0, 0.7))
    k <- kits_subj_table(parms)
    design <- data.frame(pattern_id = "all", item_id = unique(parms$item_id))
    d <- suppressWarnings(dexter::ability_tables(parms, design = design, method = "WLE"))
    if (!all(is.finite(d$theta))) {
      return(NA_real_)
    }
    expect_true(all(k$se <= d$se + 1e-12))   # dexter 的 SE 是 Warm 的 sqrt(I + (I'/2I)^2) / I
    max(abs(k$theta - d$theta))
  }, numeric(1))
  expect_gte(sum(!is.na(diffs)), 15)
  expect_lt(max(diffs, na.rm = TRUE), 5e-8)
})

# ---- CML：梯度的公式 -------------------------------------------------------------------------------------------

test_that("CML 的梯度：RTMB 的梯度、數值梯度、和公式 n_ig - sum_br n_br P(X_i >= g | b, r)（列舉）相同", {
  betas <- list(a1 = c(-1, 0.5), a2 = 0.3, a3 = c(-0.5, 0.5, 1.5))
  X <- kits_sim(rnorm(300), betas, seed = 4)
  X[1:80, "a2"] <- NA   # two booklets
  stats <- kits_stats(X)
  inp <- kits:::cml_inputs(stats)
  con <- kits:::cml_constraint(inp, NULL)
  nll <- kits:::cml_nll(inp, con)
  obj <- RTMB::MakeADFun(nll, list(b = numeric(ncol(con$A))), silent = TRUE)
  set.seed(5)
  b <- rnorm(ncol(con$A), 0, 0.5)
  beta <- con$beta0 + drop(con$A %*% b)
  # the formula, with P(X_i >= g | b, r) from the enumeration of all response patterns of each booklet
  formula_grad <- inp$n_ge
  for (bk in seq_along(inp$book_items)) {
    its <- inp$book_items[[bk]]
    patterns <- as.matrix(expand.grid(lapply(inp$m[its], \(m) 0:m)))
    steps <- lapply(its, \(i) which(inp$step_item == inp$items[i]))
    weight <- exp(-rowSums(vapply(seq_along(its), \(j) {
      vapply(patterns[, j], \(x) sum(beta[steps[[j]]][seq_len(x)]), numeric(1))
    }, numeric(nrow(patterns)))))
    raw <- rowSums(patterns)
    gamma <- tapply(weight, raw, sum)
    n_r <- inp$n_raw[[bk]]
    for (j in seq_along(its)) {
      for (g in seq_along(steps[[j]])) {
        reach <- tapply(weight * (patterns[, j] >= g), raw, sum) / gamma   # P(X_i >= g | b, r)
        idx <- steps[[j]][g]
        formula_grad[idx] <- formula_grad[idx] - sum(n_r * reach)
      }
    }
  }
  numeric_grad <- vapply(seq_along(b), \(p) {
    h <- 1e-5
    up <- b
    down <- b
    up[p] <- b[p] + h
    down[p] <- b[p] - h
    (obj$fn(up) - obj$fn(down)) / (2 * h)
  }, numeric(1))
  rtmb_grad <- drop(obj$gr(b))
  expect_lt(max_diff(rtmb_grad, drop(crossprod(con$A, formula_grad))), 1e-9)
  expect_lt(max_diff(rtmb_grad, numeric_grad), 1e-5)
  # at the estimate the score equations hold: n_ig = sum_br n_br P(X_i >= g | b, r) for the free steps
  fit <- kits_cml(stats)
  expect_lt(max(abs(obj$gr(drop(solve(crossprod(con$A), crossprod(con$A, fit$est$beta)))))), 1e-6)
})
