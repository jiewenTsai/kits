# 物件、印出、輸入檢查（不需要 kitsPreview）
set.seed(1)

betas <- list(a1 = -2, a2 = -1, a3 = c(-0.5, 0.5), a4 = 0, a5 = c(0.5, 1.5), a6 = 2, a7 = 3)
X03 <- kits_sim(rnorm(500, 0), betas[1:5], seed = 1)
X36 <- kits_sim(rnorm(500, 2), betas[3:7], prefix = "q", seed = 2)
item03 <- kits_item(X03)
item36 <- kits_item(X36)

test_that("kits_item：門檻表、fit、人數；印出來簡短；coef() 是門檻表；原點是 beta 的和 0", {
  expect_s3_class(item03, "kits_item")
  expect_s3_class(item03$fit, "kits_fit")
  expect_equal(names(coef(item03)), c("item_id", "item_score", "beta", "beta_se", "thr", "thr_se"))
  expect_equal(coef(item03), item03$thresholds)
  expect_equal(item03$n_persons, 500)
  expect_equal(sum(coef(item03)$beta), 0, tolerance = 1e-10)
  expect_output(print(item03), "500 person-waves, 5 items, 7 steps")
  expect_output(print(item03$fit), "kits_fit")
  expect_error(kits_item(X03, method = "Bayes"), "CML")
  expect_error(kits_item(X03, item_id != "a1"), "predicate")
  expect_error(kits_item(X03, nDraws = 10), "unknown argument")
  expect_true(all(tapply(coef(item03)$thr, coef(item03)$item_id, function(t) all(diff(t) > 0))))
  two_point <- coef(item03)$item_id %in% c("a1", "a2", "a4")
  expect_equal(coef(item03)$thr[two_point], coef(item03)$beta[two_point], tolerance = 1e-8)
})

test_that("kits_item：分數不連續時停下來；沒有 0 分的題停下來；固定的題沒有 0 分時不放進校準", {
  gap <- X03
  gap[gap[, "a3"] == 1, "a3"] <- 2L
  expect_error(kits_item(gap), "consecutive")
  Y <- X03
  Y[, "a1"] <- pmax(Y[, "a1"], 1L)
  expect_error(kits_item(Y), "no 0 score")
  expect_message(upd <- kits_item(Y, fixed_params = item03), "a1")
  expect_equal(upd$left_out, "a1")
  expect_equal(coef(upd)[coef(upd)$item_id == "a1", "thr"], coef(item03)[coef(item03)$item_id == "a1", "thr"], tolerance = 1e-6)
  expect_true(all(is.finite(kits_subj(Y, upd)$theta)))
  z <- X03
  z[] <- 0L
  expect_error(suppressMessages(kits_item(z, fixed_params = item03)), "nothing to calibrate")
  fixed <- subset(coef(item03), item_id != "a3")
  upd2 <- kits_item(X03, fixed_params = fixed)
  expect_equal(upd2$thresholds$beta[upd2$thresholds$item_id == "a1"], item03$thresholds$beta[item03$thresholds$item_id == "a1"])
})

test_that("kits_subj：欄位、順序、raw_max、flag；收 kits_item、fit、門檻表；長格式也收", {
  s <- kits_subj(X03, item03)
  expect_equal(names(s),
               c("pattern_id", "visit_id", "raw", "raw_max", "theta", "se", "infit", "outfit",
                 "flag"))
  expect_equal(s$visit_id, rownames(X03))
  expect_equal(s$pattern_id, rep("bk1", 500))
  expect_equal(unique(s$raw_max), 7)
  expect_equal(kits_subj(X03, item03$fit)$theta, s$theta)
  expect_equal(kits_subj(X03, coef(item03))$theta, s$theta)
  expect_equal(s$flag[s$raw == 0][1], "zero")
  expect_true(all(s$flag[s$raw == s$raw_max] == "perfect"))
  expect_true(all(is.na(s$infit[s$flag %in% c("zero", "perfect")])))
  long <- data.frame(visit_id = rep(rownames(X03), ncol(X03)),
                     item_id = rep(colnames(X03), each = nrow(X03)), item_score = as.vector(X03))
  expect_equal(kits_subj(long, item03)$theta, s$theta)
  expect_error(kits_subj(data.frame(visit_id = "p1", item_id = "a1"), item03), "item_score",
               class = "kits_error_input")
  expect_error(kits_subj(unname(X03), item03), "row names")
  expect_error(kits_subj(cbind(X03, zz = 0L), item03), "zz")
  expect_error(kits_subj(X03[1, ], item03), "drop = FALSE")
  expect_equal(nrow(kits_subj(X03[1, , drop = FALSE], item03)), 1)
  p <- coef(item03)
  p$beta[1] <- NA
  expect_error(kits_subj(X03, p), "finite")
  m <- X03[1:6, ]
  m[1, "a3"] <- 5L
  expect_error(kits_subj(m, item03), "a3:5")
})

test_that("WLE 是 Warm 的加權概似最大值：梯度在估計值是 0", {
  s <- kits_subj(X03, item03)
  parms <- coef(item03)
  target <- function(theta, raw) {   # log L + 1/2 log I
    out <- raw * theta
    info <- 0
    for (i in unique(parms$item_id)) {
      b <- parms$beta[parms$item_id == i]
      p <- kits:::pcm_probs(theta, b)
      k <- seq_along(p) - 1
      out <- out - log(sum(exp(k * theta - c(0, cumsum(b)))))
      info <- info + sum(k^2 * p) - sum(k * p)^2
    }
    out + 0.5 * log(info)
  }
  for (raw in c(0, 3, 7)) {
    th <- s$theta[s$raw == raw][1]
    best <- stats::optimize(function(t) target(t, raw), c(th - 2, th + 2), maximum = TRUE, tol = 1e-10)$maximum
    expect_equal(th, best, tolerance = 1e-6)
  }
})

test_that("參數整體偏移很大時仍找到正確的 WLE；每個組合只算一次（很快）", {
  p <- coef(item03)
  p$beta <- p$beta + 60
  shifted <- kits_subj(X03[1:10, ], p)
  expect_equal(shifted$theta, kits_subj(X03[1:10, ], item03)$theta + 60, tolerance = 1e-6)
  big <- kits_sim(rnorm(3000), betas[1:5], seed = 9)
  expect_lt(system.time(kits_subj(big, p))[["elapsed"]], 20)
})

test_that("kits_subj_table(item)：每個題本 0 到滿分各一列，包括沒放進校準的固定題", {
  tab <- kits_subj_table(item03)
  expect_equal(names(tab), c("pattern_id", "raw", "theta", "se"))
  expect_equal(tab$raw, 0:7)
  expect_true(all(diff(tab$theta) > 0))
  parms <- data.frame(item_id = c("a1", "a1", "a2", "a3", "a3"), item_score = c(1, 2, 1, 1, 2), beta = c(-1, 0, 0.5, 0, 1))
  expect_equal(kits_subj_table(parms)$raw, 0:5)
  expect_error(kits:::wle_table(parms, design = data.frame(pattern_id = "b", item_id = "zz")), "zz")
  Y <- X03
  Y[, "a1"] <- 1L
  lo <- suppressMessages(kits_item(Y, fixed_params = item03))
  tab <- kits_subj_table(lo)
  s <- kits_subj(Y, lo)
  expect_equal(max(tab$raw), max(s$raw_max))
  expect_equal(tab$theta[match(s$raw, tab$raw)], s$theta, tolerance = 1e-8)
})

test_that("kits_subj_table：人數加總；wave_m 照列或用名字；只有標準題本", {
  tab <- subj_table_of(item03, X03, wave_m = rep(12, 500))
  expect_equal(sum(tab$n), 500)
  expect_equal(tail(tab$cum_pct, 1), 100)
  expect_false("form" %in% names(tab))
  expect_true(all(is.na(subj_table_of(item03, X03)$wave_m)))   # 沒給 wave_m：所有列是一組
  Y <- X03
  Y[1:100, "a4"] <- NA
  std <- subj_table_of(item03, Y, wave_m = rep(12, 500))
  expect_equal(sum(std$n), 400)
  expect_equal(unique(std$form_pct), 80)
  expect_equal(attr(std, "forms")$n, c(400, 100))
  expect_equal(attr(std, "forms")$type, c("standard", "missing responses"))
})

test_that("kits_itemfit：雜訊題的 infit、outfit 比較大，rir 比較小；多題本也可以；一題時停下來", {
  b <- list(a1 = c(-1, 0), a2 = c(0, 1), a3 = 0.5, a4 = -0.5, a5 = 0)
  Y <- kits_sim(seq(-3, 3, length.out = 400), b, seed = 4)
  Y[, "a5"] <- sample(0:1, 400, replace = TRUE)
  Y[1:100, "a4"] <- NA
  f <- kits_itemfit(Y, kits_item(Y))
  expect_equal(names(f), c("item_id", "n", "pvalue", "rir", "infit", "outfit", "flag"))
  expect_equal(f$item_id, names(b))
  expect_equal(f$item_id[which.max(f$infit)], "a5")
  expect_equal(f$item_id[which.min(f$rir)], "a5")
  m <- X03
  m[500, c("a1", "a2")] <- NA   # 一個孩子自己一個題本
  expect_true(all(is.finite(kits_itemfit(m, item03)$rir)))
  expect_error(kits_itemfit(X03[, "a1", drop = FALSE], item03), "at least 2")
})

test_that("kits：校準加計分；印出來簡短；kits 物件可以當 item", {
  k <- kits(X03)
  expect_s3_class(k, "kits")
  expect_equal(k$subj, kits_subj(X03, k$item), ignore_attr = "info")   # info$call 是 kits(X03)
  expect_output(print(k), "Abilities")
  expect_equal(kits_subj(X03, k)$theta, k$subj$theta)
  expect_equal(kits_linking(item36, k, n_boot = 0)$intercept,
               kits_linking(item36, item03, n_boot = 0)$intercept)
})

test_that("報表印小數點下 2 位，digits 可以改；物件裡的數值不變", {
  thr <- coef(item03)
  expect_s3_class(thr, "kits_table")
  printed <- capture.output(print(thr))
  expect_true(any(grepl(sprintf("%.2f", thr$thr[1]), printed, fixed = TRUE)))
  expect_false(any(grepl(sprintf("%.4f", thr$thr[1]), printed, fixed = TRUE)))
  expect_true(any(grepl(sprintf("%.4f", thr$thr[1]), capture.output(print(thr, digits = 4)), fixed = TRUE)))
  expect_error(print(thr, digits = -1), "digits")
  link <- kits_linking(item36, item03, n_boot = 0)
  expect_output(print(link), "new value = intercept + slope x old value", fixed = TRUE)
  expect_output(print(link), sprintf("intercept = %.2f\n", link$intercept), fixed = TRUE)
  expect_output(print(link), sprintf("slope     = %.2f\n", link$slope), fixed = TRUE)
  expect_output(print(link, digits = 4), sprintf("slope     = %.4f", link$slope), fixed = TRUE)
  expect_output(print(item03, digits = 3), sprintf("%.3f", thr$thr[1]), fixed = TRUE)
  expect_error(capture.output(print(item03, digits = NA)), "digits")
})


test_that("kits_subj_table：其他題目組合分成少答（NA）和其他題；其他題是診斷裡的警告；印出來有註記", {
  set.seed(4)
  X <- kits_sim(rnorm(300), list(a1 = c(-1, 0), a2 = 0.5, a3 = -0.5, a4 = c(0, 1)), seed = 4)
  item <- kits_item(X)
  Y <- X
  Y[1:20, "a4"] <- NA                       # 少答一題
  Y <- cbind(Y, a5 = NA)
  Y[21:30, "a5"] <- 1                        # 多答一題（標準題本以外）
  Y[21:30, "a5"][1:3] <- 0
  item5 <- kits_item(Y)
  tab <- subj_table_of(item5, Y, wave_m = rep(12, 300))
  sets <- attr(tab, "forms")
  expect_equal(sets$type, c("standard", "missing responses", "other items"))
  expect_equal(sets$n, c(270, 20, 10))
  expect_equal(sum(tab$n), 270)
  out <- capture.output(print(tab))
  expect_true(any(grepl("complete responses to the standard form only", out)))
  expect_true(any(grepl("Wave 12: 20 of 300 person-waves with missing responses are not in the table", out)))
  expect_true(any(grepl("Wave 12: 10 of 300 person-waves answered items outside the standard form \\(a5\\)", out)))
  expect_false(any(grepl("forms = ", out)))
  expect_true(any(grepl("Checks: 1 warning", out)))
  dg <- attr(tab, "diagnostics")
  expect_equal(dg$status[dg$check == "other items, wave 12"], "warn")
  expect_match(dg$message[dg$check == "other items, wave 12"], "check their wave and item pattern")
  expect_equal(kits_check(tab)$status, "warn")
})

test_that("長的表只印前 10 列，最後一行說還有幾列；短的表全印；n 可以改；head() 和 [ 照舊", {
  s <- kits_subj(X03, item03)   # 500 列
  out <- capture.output(print(s))
  expect_true(any(grepl("^10 ", out)))
  expect_false(any(grepl("^11 ", out)))
  expect_true(any(grepl("... 490 more rows (more: print(x, n = 50); all rows: as.data.frame(x))", out, fixed = TRUE)))
  expect_true(any(grepl("Checks:", out)))
  out3 <- capture.output(print(s, n = 3))
  expect_true(any(grepl("... 497 more rows", out3, fixed = TRUE)))
  expect_false(any(grepl("more rows", capture.output(print(s, n = Inf)))))
  expect_false(any(grepl("more rows", capture.output(print(s[1:30, ])))))   # 30 列以內全印
  expect_true(any(grepl("... 21 more rows", capture.output(print(s[1:31, ])), fixed = TRUE)))   # 31 列：只印前 10 列
  expect_equal(nrow(head(s, 40)), 40)
  expect_s3_class(head(s, 40), "kits_subj")
  expect_equal(attr(s[1:5, ], "scale"), attr(s, "scale"))
  expect_error(print(s, n = -1), class = "kits_error_input")
  # 其他的長表：kits_table、kits_norms、kits_subj_table 用同一個規則（直接叫 kits 的方法：kitsPreview 也有 print.kits_table）
  pts <- kits_linking(item36, item03, n_boot = 0)$points
  long <- rbind(pts, pts, pts, pts, pts, pts, pts)
  expect_true(any(grepl("more rows", capture.output(kits:::print.kits_table(long)), fixed = TRUE)))
  expect_false(any(grepl("more rows", capture.output(kits:::print.kits_table(coef(item03))))))
  norm <- structure(data.frame(age_m = 1:40, mu = 1, sigma = 1, epsilon = 0, delta = 1), class = c("kits_norms", "kits_table", "data.frame"))
  expect_true(any(grepl("... 30 more rows", capture.output(print(norm)), fixed = TRUE)))
})

test_that("fixed_params：固定題這次有更高的分數時，舊的步驟固定、新的步驟估計（一波接一波）", {
  betas <- list(a1 = c(-1, 0), a2 = c(0, 1), a3 = 0.5, a4 = -0.5)
  X <- kits_sim(seq(-3, 3, length.out = 400), betas, seed = 1)
  item <- kits_item(X, scale = "lan M03")
  betas_2 <- c(betas[c("a1", "a2", "a4")], list(a3 = c(0.5, 2)))
  X_2 <- kits_sim(seq(-1, 5, length.out = 600), betas_2, prefix = "m", seed = 3)
  item_2 <- kits_item(X_2, fixed_params = coef(item))
  thr <- coef(item_2)
  old <- merge(coef(item), thr, by = c("item_id", "item_score"))
  expect_equal(old$beta.x, old$beta.y)   # 舊的 beta 一點都不變
  a3 <- thr[thr$item_id == "a3", ]
  expect_identical(a3$item_score, 1:2)
  expect_gt(abs(a3$thr[1] - old$thr.x[old$item_id == "a3"]), 0.01)   # 門檻跟著所有的步驟變
  expect_equal(a3$beta_se[1], 0)
  expect_gt(a3$beta_se[2], 0)
  expect_lt(abs(a3$beta[2] - 2), 0.5)   # 新的步驟估得到真值附近
  expect_identical(item_2$scale, "lan M03")
  # 和舊的波次一起計分：舊的人次照樣算得出來（a3 只有 0、1 分）
  expect_no_error(kits_subj(X, item_2))
})

test_that("fixed_params：新的最高分的固定題沒有 0 分也校準；新的步驟分數有缺時停下來", {
  betas <- list(a1 = c(-1, 0), a2 = c(0, 1), a3 = 0.5, a4 = -0.5)
  item <- kits_item(kits_sim(seq(-3, 3, length.out = 400), betas, seed = 1), scale = "lan M03")
  X_2 <- kits_sim(seq(-1, 5, length.out = 600), c(betas[-3], list(a3 = c(0.5, 2))), seed = 3)
  Y <- X_2
  Y[, "a3"] <- pmax(Y[, "a3"], 1L)   # a3 沒有 0 分，有新的 2 分
  item_y <- kits_item(Y, fixed_params = coef(item))
  expect_length(item_y$left_out, 0)
  expect_gt(coef(item_y)$beta_se[coef(item_y)$item_id == "a3"][2], 0)
  expect_true(all(is.finite(kits_subj(Y, item_y)$theta)))
  only_2 <- Y
  only_2[, "a3"] <- 2L   # 沒有人在新的步驟以下：估不出來
  expect_error(kits_item(only_2, fixed_params = coef(item)), "a3:2", class = "kits_error_items")
  gap <- X_2
  gap[gap[, "a3"] == 2, "a3"] <- 3L   # 固定 1，資料 0、1、3：缺 2
  expect_error(kits_item(gap, fixed_params = coef(item)), "a3 is fixed up to step 1.*missing 2",
               class = "kits_error_data")
  covered <- X_2
  covered[covered[, "a3"] == 1, "a3"] <- 0L   # 資料 0、2：缺的 1 是固定的步驟
  item_c <- kits_item(covered, fixed_params = coef(item))
  expect_identical(coef(item_c)$item_score[coef(item_c)$item_id == "a3"], 1:2)
})
