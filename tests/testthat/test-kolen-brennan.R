# Kolen & Brennan (2014) 的第一批功能：連結誤差（在每個 theta、一串連結）、斜率和 1 比、
# 垂直量尺的統計、量尺上的精準度、原始分對照（kits_subj_table() 的 WLE 和 kits_expected_score()）。
# 每一項都和手算或暴力法比。

set.seed(1)   # 測試檔可以分開跑（平行）：亂數從固定的起點開始
betas <- list(a1 = -2, a2 = -1, a3 = c(-0.5, 0.5), a4 = 0, a5 = c(0.5, 1.5), a6 = 2, a7 = 3)
X03 <- kits_sim(rnorm(500, 0), betas[1:5], seed = 1)
X36 <- kits_sim(rnorm(500, 2), betas[3:7], prefix = "q", seed = 2)
item03 <- kits_item(X03, scale = "lan M03")
item36 <- kits_item(X36, scale = "lan M36")
anchors <- data.frame(item_id = c("a1", "a2", "a3", "a4", "a5"), tau = c(40, 47, 51, 55, 60))
co <- kits_linking(item36, item03, n_boot = 200)
d_link <- kits_linking(item03, anchors, n_boot = 200)

# PCM 的類別機率（手算）
pcm_p <- function(theta, b) {
  eta <- c(0, cumsum(theta - b))
  exp(eta) / sum(exp(eta))
}

test_that("n_boot > 0：連結留下每次重抽的截距和斜率和它們的共變異數；n_boot = 0 時什麼都不變", {
  expect_true(is.matrix(co$draws))
  expect_equal(colnames(co$draws), c("intercept", "slope"))
  expect_equal(nrow(co$draws), 200 - co$n_failed)
  expect_equal(sd(co$draws[, "intercept"]), co$intercept_error)
  expect_equal(sd(co$draws[, "slope"]), co$slope_error)
  expect_equal(cov(co$draws[, "intercept"], co$draws[, "slope"]), co$vcov["intercept", "slope"])
  none <- kits_linking(item36, item03)
  expect_null(none$draws)
  expect_null(none$vcov)
  expect_false("slope vs 1" %in% none$diagnostics$check)
  # 每次重抽的線就是對重抽的點重做一次 mean/sigma（含篩點）：用同樣的亂數重做第一次（第一次沒有失敗，所以是 draws 的第 1 列）
  pts <- co$points
  items <- unique(pts$item)
  rows <- split(seq_along(pts$item), factor(pts$item, levels = items))
  set.seed(1)
  i <- unlist(rows[sample(items, replace = TRUE)], use.names = FALSE)
  first <- algo_screen(pts$from[i], pts$to[i], pts$item[i], cut = 2) |>
    algo_mean_sigma()
  expect_identical(unname(co$draws[1, ]), c(first$intercept, first$slope))
})

test_that("選定的 theta 的 se_link：一條線時是 a theta + b 在重抽之間的 SD（delta 法）；0 時是 intercept_error", {
  th <- c(-2, -0.5, 0, 1, 3)
  conv <- kits_convert(data.frame(theta = th), co, from = "lan M36")
  se <- link_se_of(co, theta = th, from = "lan M36")
  brute <- vapply(th, function(t) sd(co$draws[, "intercept"] + co$draws[, "slope"] * t), numeric(1))
  expect_equal(se$se_link, brute, tolerance = 1e-12)
  v_int <- var(co$draws[, "intercept"])
  v_slope <- var(co$draws[, "slope"])
  delta <- sqrt(v_int + th^2 * v_slope + 2 * th * co$vcov[1, 2])
  expect_equal(se$se_link, delta, tolerance = 1e-12)
  expect_equal(se$se_link[th == 0], co$intercept_error)
  expect_equal(se$value, co$intercept + co$slope * th)
  expect_equal(attr(conv, "scale"), "lan M03")
  no_draws <- kits_linking(item36, item03)
  expect_no_message(plain <- kits_convert(data.frame(theta = 0), no_draws, from = "lan M36"))
  expect_false("se_link" %in% names(plain))   # 沒有重抽：沒有 se_link
  expect_message(kits_convert(data.frame(theta = 0), no_draws, d_link, from = "lan M36"),
                 "no bootstrap draws",
                 class = "kits_message_scale")
  expect_error(link_se_of(co, d_link, theta = 0, from = "lan M03"), class = "kits_error_scale")
})

test_that("一串連結：每個重抽和下一條線的每個重抽配起來（全部的配法），和暴力法一致；分兩次換算結果相同", {
  for (x in c(-1, 0.7, 2.5)) {
    y1 <- co$draws[, "intercept"] + co$draws[, "slope"] * x                              # 每個 co 的重抽
    # 配上每個 d_link 的重抽
    y <- outer(y1, d_link$draws[, "slope"]) + rep(d_link$draws[, "intercept"], each = length(y1))
    expect_equal(link_se_of(co, d_link, theta = x, from = "lan M36")$se_link, sd(as.vector(y)),
                 tolerance = 0.03)
    # 手算的遞迴：中心在估計的線上
    v1 <- var(y1)
    m1 <- co$intercept + co$slope * x
    v_int <- var(d_link$draws[, "intercept"])
    v_slope <- var(d_link$draws[, "slope"])
    v2 <- v_int + 2 * d_link$vcov[1, 2] * m1 + v_slope * (v1 + m1^2) + d_link$slope^2 * v1
    expect_equal(link_se_of(co, d_link, theta = x, from = "lan M36")$se_link, sqrt(v2),
                 tolerance = 1e-12)
  }
  s <- kits_subj(X36, item36)
  one <- kits_convert(s, co, d_link)
  two <- kits_convert(kits_convert(s, co), d_link)
  expect_equal(two$se_link, one$se_link, tolerance = 1e-12)
  expect_equal(one$se_link, link_se_of(co, d_link, theta = s$theta, from = "lan M36")$se_link,
               tolerance = 1e-12)
  expect_equal(one$se_total, sqrt(one$se^2 + one$se_link^2))
})

test_that("kits_convert()：每條連結都有重抽時才加 se_link、se_total；門檻表用 thr_se；沒有時不加、不警告", {
  s <- kits_subj(X36, item36)
  expect_false("se_link" %in% names(kits_convert(s, kits_linking(item36, item03), d_link)))
  expect_warning(kits_convert(s, kits_linking(item36, item03), d_link), NA)
  expect_false("se_link" %in% names(kits_convert(kits_convert(s, co), kits_linking(item03, anchors))))   # 後面的連結沒有重抽：拿掉
  thr <- kits_convert(coef(item36), co)
  expect_equal(thr$se_link, link_se_of(co, theta = coef(item36)$thr, from = "lan M36")$se_link)
  expect_equal(thr$se_total, sqrt(thr$thr_se^2 + thr$se_link^2))
  tab <- kits_convert(subj_table_of(item36, X36), co, d_link)
  expect_true(all(c("se_link", "se_total") %in% names(tab)))
  expect_null(attr(kits_convert(c(-1, 0, 1), co, from = "lan M36"), "se_link"))
  # rbind：只有一部分的表有 se_link 時，其他的補 NA
  s03 <- kits_convert(kits_subj(X03, item03), kits_linking(item03, anchors))
  both <- rbind(s03, kits_convert(s, co, d_link))
  expect_equal(nrow(both), nrow(s03) + nrow(s))
  expect_true(all(is.na(both$se_link[seq_len(nrow(s03))])))
})

test_that("斜率和 1 比：(b - 1) / slope_error；|z| > 2 時註記；錨題（D）不比", {
  z <- (co$slope - 1) / co$slope_error
  row <- co$diagnostics[co$diagnostics$check == "slope vs 1", ]
  expect_equal(row$value, z)
  expect_equal(row$status, if (abs(z) > 2) "note" else "ok")
  expect_false("slope vs 1" %in% d_link$diagnostics$check)
  # 一個明顯伸縮的連結：to = 1.5 from
  set.seed(3)
  from <- rnorm(20)
  stretch <- link_of_points(from, 0.3 + 1.5 * from + rnorm(20, 0, 0.05), id = paste0("p", 1:20),
                            n_boot = 100)
  row <- stretch$diagnostics[stretch$diagnostics$check == "slope vs 1", ]
  expect_equal(row$status, "note")
  expect_match(row$message, "stretches or shrinks")
  same <- link_of_points(from, 0.3 + from + rnorm(20, 0, 0.05), id = paste0("p", 1:20),
                         n_boot = 100)
  expect_equal(same$diagnostics$status[same$diagnostics$check == "slope vs 1"], "ok")
})

test_that("Lord-Wingersky 遞迴：原始分的分佈和逐一列舉作答型態相同；E[theta_hat | theta] 也相同", {
  b <- list(x1 = c(-1, 0.5), x2 = 0.3, x3 = c(-0.5, 0, 1))
  parms <- data.frame(item_id = rep(names(b), lengths(b)), item_score = sequence(lengths(b)), beta = unlist(b))
  th <- c(-2, 0, 1.3)
  rec <- algo_lw(th, b)
  pats <- expand.grid(x1 = 0:2, x2 = 0:1, x3 = 0:3)
  wle <- kits_subj_table(parms)$theta
  for (t in seq_along(th)) {
    p <- apply(pats, 1, function(x) prod(pcm_p(th[t], b$x1)[x[1] + 1], pcm_p(th[t], b$x2)[x[2] + 1], pcm_p(th[t], b$x3)[x[3] + 1]))
    raw <- rowSums(pats)
    expect_equal(rec[t, ], as.vector(tapply(p, factor(raw, levels = 0:6), sum)), tolerance = 1e-12)
    pr <- kits_precision(parms, list(f = names(b)), theta = th)
    expect_equal(pr$expected[t], sum(p * wle[raw + 1]), tolerance = 1e-12)
    expect_equal(pr$sd[t], sqrt(sum(p * (wle[raw + 1] - sum(p * wle[raw + 1]))^2)), tolerance = 1e-10)
  }
  expect_equal(pr$rmse, sqrt(pr$bias^2 + pr$sd^2))
  expect_equal(pr$bias, pr$expected - pr$theta)
})

test_that("kits_precision()：WLE 和 kits_subj_table() 相同；連結換算；可用範圍是 RMSE < criterion 的連續區間", {
  b <- list(c1 = c(-3, -2), c2 = c(-2, -1), c3 = c(-1, 0), c4 = c(0, 1), c5 = c(1, 2), c6 = c(2, 3), c7 = -1, c8 = 1)
  Xb <- kits_sim(rnorm(800, 0, 1.5), b, seed = 5)
  it <- kits_item(Xb, scale = "lan M03")
  tab <- subj_table_of(it, Xb, wave_m = rep(c(12, 24), each = 400))
  pr <- kits_precision(it, tab, theta = seq(-6, 6, by = 0.05), criterion = 0.8)
  expect_equal(unique(pr$form), c("wave 12", "wave 24"))
  # E 只用到 subj_table 的 theta：theta 很大時 E 接近滿分的 WLE
  far <- kits_precision(it, tab, theta = c(0, 12))
  expect_equal(far$expected[2], max(tab$theta[tab$wave_m == 12]), tolerance = 1e-3)
  r <- attr(pr, "ranges")
  d <- pr[pr$form == "wave 12", ]
  inside <- d$theta >= r$from[1] & d$theta <= r$to[1]
  expect_true(all(d$rmse[inside] < 0.8))
  first <- which(d$theta == r$from[1])
  expect_gt(first, 1)
  expect_true(d$rmse[first - 1] >= 0.8)
  # 換到 D：theta 和 E 照線換，bias、SD、RMSE 乘上斜率；criterion 預設 1 logit（0.15.0 起；之前是 0.5）
  dl <- kits_linking(it, data.frame(item_id = paste0("c", 1:6), tau = 50 + 7 * coef(it)$thr[coef(it)$item_score == 1][1:6]))
  pd <- kits_precision(it, tab, theta = seq(-6, 6, by = 0.05), link = dl,
                       criterion = 0.8 * dl$slope)
  expect_equal(pd$theta, dl$intercept + dl$slope * pr$theta)
  expect_equal(pd$rmse, dl$slope * pr$rmse)
  expect_equal(attr(kits_precision(it, tab, link = dl), "criterion"), 1 * dl$slope)   # 預設 1 logit
  expect_equal(attr(pd, "ranges")$from, dl$intercept + dl$slope * r$from)
  expect_equal(attr(pd, "scale"), "D")
  expect_error(kits_precision(it, list(f = c("c1", "zz"))), class = "kits_error_items")
  expect_error(kits_precision(it, 3), class = "kits_error_input")
  expect_output(print(pr), "Usable range")
  expect_s3_class(plot(pr), "plotly")
})

test_that("原始分對照：kits_subj_table() 的 WLE 代進另一組題的 kits_expected_score()，和手算相同", {
  b <- list(a1 = c(-2, -1), a2 = c(-1, 0), a3 = c(0, 1), a4 = c(1, 2), a5 = 2.5)
  X <- kits_sim(rnorm(500, 0, 1.5), b, seed = 1)
  X[1:40, c("a4", "a5")] <- NA                     # 第二個作答組合：只答 a1-a3
  item <- kits_item(X)
  # 兩行：A = 只答 a1-a3 的作答組合，B = 全部 5 題
  tab <- kits_subj_table(item)
  short <- tab[tab$pattern_id == "bk1", ]
  eq <- kits_expected_score(item, items = names(b))(short$theta)
  # 手算：估計的 beta、PCM 的動差、WLE 的估計式 r - E + T3 / (2 I) = 0（Warm, 1989），再算 B 的期望分數
  pr <- coef(item)
  bh <- split(pr$beta, pr$item_id)
  mom <- function(t, its) {
    m <- vapply(bh[its], function(beta) {
      p <- pcm_p(t, beta)
      k <- seq_along(p) - 1
      e <- sum(k * p)
      c(E = e, I = sum((k - e)^2 * p), T3 = sum((k - e)^3 * p))
    }, numeric(3))
    rowSums(m)
  }
  wle <- vapply(0:6, function(r) {
    f <- function(t) {
      m <- mom(t, c("a1", "a2", "a3"))
      r - m[["E"]] + m[["T3"]] / (2 * m[["I"]])
    }
    uniroot(f, c(-15, 15), tol = 1e-12)$root
  }, numeric(1))
  expect_equal(short$raw, 0:6)
  expect_equal(short$theta, wle, tolerance = 1e-6)
  expect_equal(eq, vapply(wle, function(t) mom(t, names(b))[["E"]], numeric(1)), tolerance = 1e-6)
  expect_true(all(diff(eq) > 0) && eq[1] > 0 && eq[7] < 9)   # 0 分和滿分也有有限的 WLE
})

test_that("kits_vertical()：平均、SD、sd_true、成長、效果量（Eq. 9.47）、p_sup 都和手算相同；世代比較；量尺縮小的註記", {
  set.seed(2)
  s <- data.frame(cohort = rep(c("M03", "M36"), c(450, 300)), wave_m = c(rep(c(12, 24, 36), each = 150), rep(c(36, 48), each = 150)))
  s$theta <- 40 + 8 * log(s$wave_m) + rnorm(750, 0, 4)
  s$se <- runif(750, 1, 2)
  s$age_m <- s$wave_m + runif(750, -1, 1)
  v <- kits_vertical(s)
  expect_equal(v$cohort, c("M03", "M03", "M03", "M36", "M36"))
  g <- s[s$cohort == "M03" & s$wave_m == 24, ]
  p <- s[s$cohort == "M03" & s$wave_m == 12, ]
  r <- v[2, ]
  expect_equal(r$mean, mean(g$theta))
  expect_equal(r$sd, sd(g$theta))
  expect_equal(r$median, median(g$theta))
  expect_equal(r$sd_true, sqrt(var(g$theta) - mean(g$se^2)))
  expect_equal(r$growth, mean(g$theta) - mean(p$theta))
  expect_equal(r$effect_size, (mean(g$theta) - mean(p$theta)) / sqrt((var(g$theta) + var(p$theta)) / 2))
  expect_equal(r$p_sup, mean(outer(g$theta, p$theta, ">")) + 0.5 * mean(outer(g$theta, p$theta, "==")))
  expect_true(is.na(v$growth[4]))                         # 每個世代的第一波沒有前一波
  co36 <- attr(v, "cohorts")
  a <- s[s$cohort == "M03" & s$wave_m == 36, "theta"]
  b <- s[s$cohort == "M36" & s$wave_m == 36, "theta"]
  expect_equal(co36$mean_diff, mean(b) - mean(a))
  expect_equal(co36$z, (mean(b) - mean(a)) / sqrt(var(a) / 150 + var(b) / 150))
  # p_sup 只看順序：單調轉換不變，效果量會變
  s2 <- s
  s2$theta <- exp(s$theta / 20)
  expect_equal(kits_vertical(s2)$p_sup, v$p_sup)
  # 量尺縮小：最後一波的 sd_true 只有前面的一半
  s3 <- s
  last <- s3$cohort == "M36" & s3$wave_m == 48
  s3$theta[last] <- mean(s3$theta[last]) + (s3$theta[last] - mean(s3$theta[last])) / 2.5
  d3 <- attr(kits_vertical(s3), "diagnostics")
  expect_equal(d3$status[d3$check == "scale shrinkage, M36"], "note")
  expect_match(d3$message[d3$check == "scale shrinkage, M36"], "p. 468")
  expect_output(print(v), "Cohorts at the same wave")
  expect_s3_class(plot(v), "plotly")
  expect_error(kits_vertical(s[c("theta", "cohort")]), class = "kits_error_input")
})

test_that("kits_vertical() 的 se_link：平均的連結誤差就是在平均上的 se_link（se_link^2 是 theta 的二次式）", {
  persons <- data.frame(visit_id = rownames(X36), wave_m = rep(c(36, 48), each = 250))
  s <- kits_convert(kits_subj(X36, item36, persons = persons), co, d_link)
  v <- kits_vertical(s)
  for (w in c(36, 48)) {
    raw <- kits_subj(X36, item36)$theta[persons$wave_m == w]
    expect_equal(v$se_link[v$wave_m == w],
                 link_se_of(co, d_link, theta = mean(raw), from = "lan M36")$se_link,
                 tolerance = 1e-8)
  }
})
