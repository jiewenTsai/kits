# 演算法層（algo_）：kits_ 的函數呼叫它們，所以數字一個位元都不差（identical()）
betas <- list(a1 = c(-1, 0), a2 = c(0, 1), a3 = 0.5, a4 = -0.5)
X <- kits_sim(seq(-3, 3, length.out = 400), betas, seed = 1)
X[1:100, "a4"] <- NA   # 兩個題本
item <- kits_item(X)
cml <- kits_stats(X) |>
  algo_cml()

test_that("kits_stats() |> algo_cml() |> algo_thresholds() 就是 kits_item() 的參數和門檻", {
  expect_identical(cml$coef$beta, item$fit$coef$beta)
  expect_identical(cml$coef$se, item$fit$coef$SE_beta)
  expect_identical(cml$vcov, item$fit$est$acov.beta)
  expect_identical(cml$nll, item$fit$est$nll)
  thr <- kits_stats(X) |>
    algo_cml() |>
    algo_thresholds()
  expect_identical(thr, as.data.frame(coef(item))[names(thr)])
  expect_identical(names(cml$betas), c("a1", "a2", "a3", "a4"))
  expect_identical(algo_thresholds(cml, rp = 0.8)$thr, coef(item, rp = 0.8)$thr)
})

test_that("algo_jme() 就是 kits_subj() 和 kits_subj_table() 的 WLE", {
  s <- kits_subj(X, item)
  items <- map(seq_len(nrow(X)), \(p) names(X[p, ])[!is.na(X[p, ])])
  wle <- algo_jme(cml$betas, raw = s$raw, items = items)
  expect_identical(wle$theta, s$theta)
  expect_identical(wle$se, s$se)
  tab <- kits_subj_table(item)
  full <- tab[tab$pattern_id == tab$pattern_id[nrow(tab)], ]   # 四題的題本
  expect_identical(algo_jme(cml$betas)$theta, full$theta)
})

test_that("algo_pcm() 和 algo_esf() 就是 kits_pcm_prob() 和 kits_esf() 的數字", {
  th <- c(-2, 0, 1.5)
  pp <- kits_pcm_prob(item, th, items = "a1")
  p <- algo_pcm(th, cml$betas$a1, third = TRUE)
  expect_identical(pp$moments$expected, p$E)
  expect_identical(pp$moments$variance, p$Var)
  expect_identical(pp$prob$prob, as.vector(t(p$P)))
  esf <- algo_esf(cml$betas)
  expect_identical(kits_esf(item)$gamma$log_gamma, esf$log_gamma)
  expect_identical(kits_esf(item)$moments$expected, as.vector(esf$expected))
})

test_that("algo_screen() |> algo_mean_sigma()、algo_link_error() 就是 kits_linking() 的數字", {
  set.seed(11)
  b2 <- list(a1 = -2, a2 = -1, a3 = c(-0.5, 0.5), a4 = 0, a5 = c(0.5, 1.5), a6 = 2, a7 = 3)
  young <- kits_item(kits_sim(rnorm(500, 0), b2[1:5], seed = 1), scale = "young")
  old <- kits_item(kits_sim(rnorm(500, 2), b2[3:7], prefix = "q", seed = 2), scale = "old")
  for (spread in c("sd", "se")) {
    link <- kits_linking(old, young, n_boot = 50, seed = 3, spread = spread, cut = 1.8)
    pts <- link$points
    s <- algo_screen(pts$from, pts$to, pts$item, cut = 1.8, spread = spread,
                     se_from = pts$se_from, se_to = pts$se_to)
    expect_identical(c(s$intercept0, s$slope0), c(link$intercept0, link$slope0))
    expect_identical(s$resid, pts$resid)
    expect_identical(s$dropped, pts$dropped)
    expect_identical(algo_mean_sigma(s), list(intercept = link$intercept, slope = link$slope))
    err <- algo_link_error(s, n_boot = 50, seed = 3)
    expect_identical(err$draws, link$draws)
    expect_identical(err$intercept_error, link$intercept_error)
    expect_identical(err$slope_error, link$slope_error)
  }
  own <- data.frame(id = paste0("p", 1:8), item = rep(paste0("i", 1:4), 2),
                    from = c(-2, -1, 0, 1, 2, 3, 0.5, 1.5), to = c(-1, 0, 1, 2, 3, 4, 1.5, 5.5))
  sc <- kits_link_screen(own, spread = "mad", unit = "item")
  s <- algo_screen(own$from, own$to, own$item, spread = "mad", unit = "item")
  expect_identical(s$dropped, sc$dropped)
  expect_identical(algo_mean_sigma(s)$slope, kits_link_line(sc)$slope)
})

test_that("algo_line()：套用、反函數、合成", {
  l1 <- list(intercept = 1, slope = 2)
  l2 <- list(intercept = -1, slope = 0.5)
  expect_identical(algo_line(l1, x = c(0, 1)), c(1, 3))
  expect_identical(algo_line(l1, invert = TRUE), list(intercept = -0.5, slope = 0.5))
  expect_identical(algo_line(list(l1, l2)), list(intercept = -0.5, slope = 1))
  expect_identical(algo_line(list(l1, l2), x = 3), algo_line(l2, x = algo_line(l1, x = 3)))
  expect_identical(algo_line(list()), list(intercept = 0, slope = 1))
})

test_that("algo_shash() 就是 kits_norm_fit() 的模型；algo_shash_z() 是常模診斷的同齡 z", {
  set.seed(1)
  age <- runif(600, 3, 36)
  D <- 20 + 15 * log(age) + rnorm(600, 0, 3)
  norm_data <- kits_norm_data(D, age, ages = seq(3, 36, by = 3))
  fit <- kits_norm_fit(norm_data, degrees = c(3, 2, 1, 1), delta_starts = 0.5)   # delta0 = 1 停在平坦處
  model <- algo_shash(sqrt(norm_data$age_m), norm_data$d, degrees = c(3, 2, 1, 1), delta0 = 0.5)
  expect_identical(model, fit$model)
  z <- model |>
    algo_shash_z(sqrt(norm_data$age_m), norm_data$d)
  diag <- attr(kits_norm_table(fit), "diagnostics")
  worst <- max(abs(z_by_age_decile(z, norm_data$age_m)$mean))
  expect_identical(diag$value[diag$check == "same-age z: max |mean| by age decile"], worst)
})

test_that("algo_lw() 就是 kits_precision() 的原始分分布", {
  th <- seq(-3, 3, by = 0.5)
  pr <- kits_precision(item, c("a1", "a2", "a3"), theta = th)
  P <- algo_lw(th, cml$betas[c("a1", "a2", "a3")])
  wle <- algo_jme(cml$betas[c("a1", "a2", "a3")])
  expect_identical(pr$expected, drop(P %*% wle$theta))
  expect_equal(rowSums(P), rep(1, length(th)), tolerance = 1e-12)
})

test_that("Warm 的原路（分數方程式 + J / 2I = 0，uniroot）的根就是 algo_jme() 和 kits_subj() 的 WLE；SE = 1/sqrt(I)", {
  # 只用來比較：P'、P''、d log S / d theta 由 RTMB 的 tape（jacfun()）給，不手推
  tapes <- map(cml$betas, \(b) {
    k <- seq(0, length(b))
    eta <- function(t) t * k - c(0, cumsum(b))
    P <- RTMB::MakeTape(\(t) exp(eta(t)) / sum(exp(eta(t))), 0)
    P1 <- P$jacfun()
    log_S <- RTMB::MakeTape(\(t) log(sum(exp(eta(t)))), 0)
    list(P = P, P1 = P1, P2 = P1$jacfun(), dlog_S = log_S$jacfun())
  })
  warm_score <- function(theta, items, r) {
    parts <- map(tapes[items], \(tp) {
      p <- drop(tp$P(theta))
      p1 <- drop(tp$P1(theta))
      c(dlog_S = drop(tp$dlog_S(theta)), I = sum(p1^2 / p), J = sum(p1 * drop(tp$P2(theta)) / p))
    })
    s <- Reduce(`+`, parts)
    r - s[["dlog_S"]] + s[["J"]] / (2 * s[["I"]])
  }
  items <- list(c("a1", "a2", "a3"), c("a1", "a2", "a3", "a4"))
  combos <- data.frame(set = rep(1:2, c(6, 7)), raw = c(0:5, 0:6))
  warm <- map2_dbl(items[combos$set], combos$raw, \(its, r) {
    stats::uniroot(warm_score, c(-15, 15), items = its, r = r, tol = 1e-12)$root
  })
  jme <- algo_jme(cml$betas, raw = combos$raw, items = items[combos$set])
  expect_lt(max(abs(warm - jme$theta)), 1e-8)
  info <- map2_dbl(items[combos$set], jme$theta, \(its, t) {   # I = sum_ik P'_ik^2 / P_ik
    sum(map_dbl(tapes[its], \(tp) sum(drop(tp$P1(t))^2 / drop(tp$P(t)))))
  })
  expect_equal(jme$se, 1 / sqrt(info), tolerance = 1e-10)   # SE = 1 / sqrt(I(theta))
  s <- kits_subj(X, item)
  four <- !is.na(X[, "a4"])
  at <- match(s$raw, combos$raw[combos$set == 2]) + 6   # 四題的題本
  expect_lt(max(abs(warm[at[four]] - s$theta[four])), 1e-8)
})
