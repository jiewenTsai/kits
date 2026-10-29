# 連結、偏移、常模（不需要 kitsPreview）
set.seed(1)

test_that("kits_mean_sigma 剔除偏移的點，並找回原本的線；輸入有問題時說清楚", {
  from <- seq(-2, 3, length.out = 12)
  to <- 1 + 2 * from
  to[5] <- to[5] + 4
  link <- link_of_points(from, to, id = paste0("p", 1:12))
  expect_s3_class(link, "kits_link")
  expect_equal(link$dropped, "p5")
  expect_equal(c(link$intercept, link$slope), c(1, 2))
  expect_true(is.na(link$intercept_error))
  expect_output(print(link), "dropped: 1 \\(p5\\)")
  expect_length(link_of_points(from, 1 + from, id = paste0("p", 1:12))$dropped, 0)
  expect_error(link_of_points(c(1, NA, 3), 1:3, letters[1:3]), "NA")
  expect_error(link_of_points(1:2, 1:2, letters[1:2]), "3 points")
  expect_error(link_of_points(rep(1, 4), 1:4, letters[1:4]), "same")
  expect_error(link_of_points(from, to, paste0("p", 1:12), cut = 0), "positive")
  expect_error(link_of_points(from, to, paste0("p", 1:12), n_boot = 1), "at least 2")
  expect_warning(link_of_points(from[1:5], -from[1:5], letters[1:5]), "negatively")
})

test_that("連結誤差：算不出來的重抽另外計數；不改變使用者的亂數狀態", {
  from <- c(-2, -1, 0.5, 2)
  to <- 1 + from + c(0.1, -0.2, 0.1, 0)
  link <- link_of_points(from, to, paste0("p", 1:4), n_boot = 300)
  expect_true(is.finite(link$intercept_error))
  expect_gt(link$n_failed, 0)
  set.seed(42)
  before <- runif(1)
  set.seed(42)
  link_of_points(from, to, paste0("p", 1:4), n_boot = 10)
  expect_equal(runif(1), before)
})

betas <- list(a1 = -2, a2 = -1, a3 = c(-0.5, 0.5), a4 = 0, a5 = c(0.5, 1.5), a6 = 2, a7 = 3)
young <- kits_item(kits_sim(rnorm(500, 0), betas[1:5], seed = 1), scale = "young")
old <- kits_item(kits_sim(rnorm(500, 2), betas[3:7], prefix = "q", seed = 2), scale = "old")

test_that("kits_linking 用共同題的所有步驟；共同題太少、表重複時停下來", {
  link <- kits_linking(old, young, n_boot = 100)
  expect_equal(link$n, 5)
  expect_equal(link$n_items, 3)
  expect_true(abs(link$slope - 1) < 0.3)
  expect_equal(kits_linking(coef(old), coef(young), n_boot = 0)$intercept,
               kits_linking(old, young, n_boot = 0)$intercept)
  expect_error(kits_linking(coef(old)[coef(old)$item_id == "a7", ], young), "common")
  expect_error(kits_linking(rbind(coef(old), coef(old)), young), "repeated")
})

test_that("kits_linking(steps, anchors)：錨題用 milestone step；step 參數；有 domain 欄時給各領域截距", {
  steps <- data.frame(item_id = paste0("a", 1:8), item_score = 1, thr = c(-1, 0, 1, 2, -1.5, 0.5, 1, 2.5),
                      domain = rep(c("lan", "cog"), each = 4))
  steps <- rbind(steps, data.frame(item_id = "a1", item_score = 2, thr = 9, domain = "lan"))
  anchors <- data.frame(item_id = paste0("a", 1:8), tau = 64 + 7 * steps$thr[1:8] + ifelse(1:8 > 4, 1, 0))
  link <- kits_linking(steps, anchors, n_boot = 0)
  expect_equal(link$n_items, 8)
  expect_equal(link$domains$slope, c(7, 7), tolerance = 1e-8)
  intercepts <- stats::setNames(link$domains$intercept, link$domains$domain)
  expect_equal(intercepts[["cog"]] - intercepts[["lan"]], 1, tolerance = 1e-8)
  expect_output(print(link), "domains")
  expect_error(kits_linking(steps, anchors[1:2, ]), "at least 3")
  expect_error(kits_linking(steps, rbind(anchors, anchors[1, ])), "repeated")
  expect_error(kits_linking(steps, anchors, step = 0), "step")
  X_cog <- kits_sim(rnorm(300), stats::setNames(betas[1:4], paste0("ca", 1:4)), seed = 3)
  two <- kits_linking(list(lan = young, cog = kits_item(X_cog)),
                        data.frame(item_id = c("a1", "a2", "a3", "ca1", "ca2"), tau = c(50, 57, 60, 50, 57)), n_boot = 0)
  expect_equal(two$domains$domain, c("cog", "lan"))
  expect_error(kits_linking(list(young, old), anchors), "name each element")
})

test_that("kits_linking 兩個校準之間：step = 1 只用 milestone step；記下 kind 和 step；預設的 to_scale", {
  all_steps <- kits_linking(old, young, n_boot = 0)
  first <- kits_linking(old, young, step = 1, n_boot = 0)
  expect_equal(first$n, 3)
  expect_equal(first$points$id, c("a3:1", "a4:1", "a5:1"))
  expect_equal(c(first$intercept, first$slope),
               unname(unlist(kits:::algo_mean_sigma(first$points$from, first$points$to))))
  expect_equal(all_steps$settings$kind, "calibration")
  expect_equal(all_steps$settings$step, "all")
  expect_equal(first$settings$step, 1)
  expect_equal(kits_info(first)$settings$kind, "calibration")
  expect_error(kits_linking(old, young, step = 2, n_boot = 0), "at step 2")
  expect_error(kits_linking(old, young, step = 1.5), class = "kits_error_input")
  # to_scale：兩個校準之間是 to 的量尺；錨題是 "D"
  o <- kits_item(kits_sim(rnorm(500, 2), betas[3:7], prefix = "q", seed = 2), scale = "lan M36")
  y <- kits_item(kits_sim(rnorm(500, 0), betas[1:5], seed = 1), scale = "lan M03")
  co <- kits_linking(o, y, n_boot = 0)
  expect_equal(c(co$from_scale, co$to_scale), c("lan M36", "lan M03"))
  expect_equal(kits_linking(o, y, n_boot = 0, to_scale = "lan M03")$to_scale, "lan M03")
  expect_error(kits_linking(o, y, n_boot = 0, to_scale = "other"), class = "kits_error_scale")
  anchors <- data.frame(item_id = c("a1", "a2", "a3", "a4", "a5"), tau = c(40, 47, 51, 55, 60))
  dl <- kits_linking(y, anchors, n_boot = 0)
  expect_equal(dl$to_scale, "D")
  expect_equal(dl$settings$kind, "anchors")
  expect_equal(dl$settings$step, 1)
  expect_equal(dl$points$step, rep(1L, 5))
  expect_equal(kits_linking(y, anchors, n_boot = 0, to_scale = "GSED D")$to_scale, "GSED D")
  expect_output(kits:::print.kits_link(first), "points from: step 1 of each common item")
  expect_output(kits:::print.kits_link(dl), "points from: the pass step of each anchor item")
})

test_that("kits_linking 兩個校準之間也能分領域：共用斜率、各領域截距；一個領域的點都被篩掉時警告", {
  b2 <- list(c1 = -1.5, c2 = -0.5, c3 = c(0, 1), c4 = 0.5, c5 = 1.5)
  cog_y <- kits_item(kits_sim(rnorm(500, 0), b2[1:4], seed = 5), scale = "cog young")
  cog_o <- kits_item(kits_sim(rnorm(500, 1.5), b2[2:5], prefix = "r", seed = 6), scale = "cog old")
  link <- kits_linking(list(lan = old, cog = cog_o), rbind(as.data.frame(coef(young)), as.data.frame(coef(cog_y))), n_boot = 0)
  expect_equal(sort(unique(link$points$domain)), c("cog", "lan"))
  expect_equal(link$points$domain, ifelse(substr(link$points$item, 1, 1) == "a", "lan", "cog"))
  expect_equal(link$domains$domain, c("cog", "lan"))
  keep <- !link$points$dropped
  shared <- kits:::mean_sigma_shared_slope(link$points$from[keep], link$points$to[keep], link$points$domain[keep])
  expect_equal(link$domains$slope, rep(shared$slope, 2))
  expect_equal(link$domains$intercept, unname(shared$intercept[c("cog", "lan")]))
  expect_true("max |domain intercept - intercept|" %in% link$diagnostics$check)
  # 印出來時說清楚：領域的表是檢查用，共用的斜率不是上面那條線的 b
  out <- capture.output(kits:::print.kits_link(link))
  expect_true(any(grepl("domain check: one slope shared by the domains, fitted with domain intercepts (kept points);", out, fixed = TRUE)))
  expect_true(any(grepl("it is not the slope of the line above", out, fixed = TRUE)))
  expect_false(any(grepl("(shared slope, same kept points)", out, fixed = TRUE)))
  # 一個領域的點全部被篩掉
  steps <- data.frame(item_id = paste0("a", 1:9), item_score = 1, thr = c(-1, 0, 1, 2, -1.5, 0.5, 1, 2.5, 0),
                      domain = c(rep("lan", 8), "cog"))
  anchors <- data.frame(item_id = steps$item_id, tau = 64 + 7 * steps$thr + c(rep(0, 8), 30))
  expect_warning(kits_linking(steps, anchors, n_boot = 0), class = "kits_warning_scale")
})

test_that("kits_drift 找到偏移的題，資訊不足的題不判斷；drift$fix 可以直接固定", {
  betas <- list(a1 = c(-1, 0), a2 = c(0, 1), a3 = 0.5, a4 = -0.5, a5 = c(-0.5, 0.5), a6 = 1.5, a7 = -3.5)
  ref <- kits_item(kits_sim(seq(-3, 3, length.out = 400), betas, seed = 1))
  new_betas <- c(modifyList(betas, list(a3 = 2)), list(a8 = c(1, 2)))
  new <- kits_sim(seq(-2, 4, length.out = 400), new_betas, prefix = "n", seed = 2)
  chk <- kits_drift(new, ref)
  expect_equal(names(chk$items), c("item_id", "max_abs_resid", "drift", "suggest_fix", "reason"))
  expect_true(chk$items$drift[chk$items$item_id == "a3"])
  expect_true(is.na(chk$items$drift[chk$items$item_id == "a7"]))
  expect_equal(chk$items$reason[chk$items$item_id == "a7"], "too little information")
  expect_equal(sum(chk$items$drift, na.rm = TRUE), 1)
  expect_false("a8" %in% chk$items$item_id)
  upd <- kits_item(new, fixed_params = chk$fix)
  expect_equal(coef(upd)$beta[1:2], coef(ref)$beta[1:2])
  expect_true(all(c("a3", "a8") %in% coef(upd)$item_id))
  wider <- new
  wider[300:400, "a3"] <- 2L
  a3 <- kits_drift(wider, ref)$items
  expect_equal(a3$reason[a3$item_id == "a3"], "score range differs")
  expect_true(a3$suggest_fix[a3$item_id == "a3"])   # 0.25.0：舊的步驟可以固定，新的步驟估計
  up3 <- kits_item(wider, fixed_params = kits_drift(wider, ref)$fix)
  expect_equal(coef(up3)$beta[coef(up3)$item_id == "a3"][1], coef(ref)$beta[coef(ref)$item_id == "a3"])
  bank <- coef(ref)[c("item_id", "item_score", "beta")]   # 只有 beta 的表也可以當 ref
  expect_equal(kits_drift(new, bank)$items, chk$items)
  rp80 <- coef(ref, rp = 0.8)[c("item_id", "item_score", "beta")]   # beta 和 RP 無關
  expect_equal(kits_drift(new, rp80)$items, chk$items)
  odd <- bank
  odd$item_score[odd$item_id == "a3"] <- 2L   # 步驟不是 1..m：門檻算不出來
  expect_error(kits_drift(new, odd), "steps 1, 2", class = "kits_error_input")
  no0 <- wider
  no0[, "a3"] <- pmax(no0[, "a3"], 1L)   # 沒有 0 分、有新的最高分：建議固定，也校準得出來
  fix0 <- kits_drift(no0, ref)$fix
  expect_true("a3" %in% fix0$item_id)
  expect_true(all(is.finite(kits_subj(no0, kits_item(no0, fixed_params = fix0))$theta)))
  expect_error(kits_drift(new, ref, min_n = 1000), "min_n")
  expect_error(kits_drift(new, coef(ref)[coef(ref)$item_id == "zz", ]), "No item")
})

test_that("shash_z 和 shash_quantile 互為反函數", {
  z <- c(-2, -0.5, 0, 1.5)
  X <- kits:::shash_quantile(z, mu = 50, sigma = 8, epsilon = 0.3, delta = 1.4)
  expect_equal(kits:::shash_z(X, mu = 50, sigma = 8, epsilon = 0.3, delta = 1.4), z)
})

test_that("kits_norms 找回模擬的中位數；幾個起始值裡留概似最大的解；kits_norm_z 認得常模", {
  set.seed(1)
  age <- runif(2000, 3, 36)
  D <- 20 + 15 * log(age) + rnorm(2000, 0, 3)
  norm <- kits_norms(D, age, ages = 6:24, degrees = c(3, 2, 1, 1))
  expect_s3_class(norm, "kits_norms")
  expect_equal(names(norm)[1:7], c("age_m", "mu", "sigma", "epsilon", "delta", "n", "neg_loglik"))
  expect_true(all(paste0("P", NORM_PCT) %in% names(norm)))
  expect_equal(norm$P50[norm$age_m %in% c(6, 12, 24)], 20 + 15 * log(c(6, 12, 24)), tolerance = 0.02)
  expect_equal(attr(norm, "kits_args")$age_x, "sqrt")
  z <- kits_norm_z(D, age, transform(norm, segment = "all"))
  expect_equal(mean(z[age >= 6 & age <= 24]), 0, tolerance = 0.1)
  one <- suppressWarnings(kits_norms(D[1:800], age[1:800], ages = 6:24, degrees = c(3, 2, 1, 2), delta_starts = 1))
  best <- suppressWarnings(kits_norms(D[1:800], age[1:800], ages = 6:24, degrees = c(3, 2, 1, 2)))
  expect_lte(best$neg_loglik[1], one$neg_loglik[1] + 1e-6)
  expect_warning(kits_norms(D, age, ages = 6:24, degrees = c(3, 2, 1, 2), delta_max = 0), "delta reaches")
  expect_error(kits_norms(D, age, ages = 3:36, degrees = c(1, 1, 0, 0)), "at least 1")
  expect_error(kits_norms(D, rep(12, 2000), ages = 3:36, degrees = c(2, 1, 1, 1)), "fewer than 2")
})

test_that("kits_norm_z 在常模範圍外回傳 NA；缺欄位時說清楚；一列的常模停下來", {
  norm <- data.frame(age_m = c(6, 12), mu = c(40, 60), sigma = 5, epsilon = 0, delta = 1)
  z <- kits_norm_z(c(40, 60, 50), age_m = c(6, 12, 30), norm = norm)
  expect_equal(z[1:2], c(0, 0))
  expect_true(is.na(z[3]))
  expect_error(kits_norm_z(1, 6, norm[, 1:3]), "epsilon")
  expect_error(kits_norm_z(40, 6, norm[1, ]), "at least 2 rows")
})

test_that("kits_norm_screen 找出堆在最低分的波次；normed 每個人次一個旗標", {
  D <- c(10, 10, 10, 12, 15, 20, 22, 25, 30, 31, NA)
  wave_m <- c(rep(c(6, 12), each = 5), 12)
  out <- kits_norm_screen(D[1:10], wave_m[1:10])$waves
  expect_equal(names(out), c("wave_m", "n", "floor_pct", "normed"))
  expect_equal(out$floor_pct, c(60, 20))
  expect_equal(out$normed, c(FALSE, FALSE))
  screen <- kits_norm_screen(D, wave_m, cut = 50)
  expect_equal(screen$normed, c(rep(FALSE, 5), rep(TRUE, 5), FALSE))
  by_cohort <- kits_norm_screen(D[1:10], wave_m[1:10], cohort = rep(c("M03", "M36"), times = 5))$waves
  expect_equal(nrow(by_cohort), 4)
  expect_output(print(kits_norm_screen(D[1:10], wave_m[1:10])), "0 of 10")
  expect_error(kits_norm_screen(D, wave_m[-1]), "same length")
})

test_that("kits_norm_ages：少數特別小的孩子不會把範圍往前拉", {
  age <- c(seq(11, 13, length.out = 500), 7, 8, seq(23, 25, length.out = 500))
  wave <- c(rep(12, 502), rep(24, 500))
  ages <- kits_norm_ages(age, wave_m = wave)
  expect_equal(range(ages), c(11, 25))
  expect_equal(diff(ages[1:2]), 0.5)
  expect_error(kits_norm_ages(c(5, 6, 11, 12), c(6, 6, 12, 12), by = 0), "positive")
  expect_error(kits_norm_ages(c(30, 31, 5, 6), c(6, 6, 12, 12)), "not above")
})

test_that("kits_norm_ages 收 kits_subj：用它的 age_m、wave_m 欄，和向量的寫法相同", {
  set.seed(2)
  X <- kits_sim(rnorm(400), list(a1 = c(-1, 0), a2 = 0.5, a3 = -0.5), seed = 2)
  persons <- data.frame(visit_id = rownames(X), wave_m = rep(c(12, 24), each = 200))
  persons$age_m <- persons$wave_m + rnorm(400, 0, 0.6)
  s <- kits_subj(X, kits_item(X), persons = persons)
  expect_equal(kits_norm_ages(s), kits_norm_ages(persons$age_m, persons$wave_m))
  expect_equal(kits_norm_ages(s, by = 1), kits_norm_ages(persons$age_m, persons$wave_m, by = 1))
  keep <- s$wave_m == 24 | s$age_m > 12
  expect_equal(kits_norm_ages(s[keep, ]), kits_norm_ages(s$age_m[keep], s$wave_m[keep]))
  expect_error(kits_norm_ages(s, persons$wave_m), "wave_m comes from its column", class = "kits_error_input")
  no_age <- kits_subj(X, kits_item(X), persons = persons[c("visit_id", "wave_m")])
  expect_error(kits_norm_ages(no_age), "age_m", class = "kits_error_input")
})

test_that("kits_norms：nlminb 收斂時不警告；收斂的標準和 kits_item() 相同；真的沒收斂時照樣警告", {
  # 0.11.0 在這份資料（和 c(4, 2, 1, 1)、c(3, 2, 2, 1) 等）誤報沒有收斂（optim code 52）；
  # kitsPreview（cNORM）不警告，P50 相同到 1e-6
  set.seed(1)
  age <- runif(2000, 3, 36)
  D <- 20 + 15 * log(age) + rnorm(2000, 0, 3)
  for (deg in list(c(3, 2, 1, 1), c(3, 2, 2, 1))) {
    expect_no_warning(norm <- kits_norms(D, age, ages = seq(4, 35, by = 0.5), degrees = deg, pad = 0))
    dg <- attr(norm, "diagnostics")
    expect_equal(dg$status[dg$check == "optimizer convergence"], "ok")
    expect_match(dg$message[dg$check == "optimizer convergence"], "^converged \\(nlminb code")
    expect_lte(dg$value[dg$check == "max |gradient|"], 1e-4)
    expect_equal(dg$status[dg$check == "max |gradient|"], "ok")
  }
  # 收斂規則本身：不在界上的係數，最大 |梯度| <= 1e-4，它們的 Hessian 正定
  lower <- c(-Inf, -10)
  upper <- c(Inf, 10)
  small <- function(p) c(1e-5, -3)                                # 第 2 個係數在上界：不算它的梯度
  he <- function(p) diag(c(1, -1))                                # 在界上的那一格不算
  conv <- kits:::shash_convergence(c(1, 10), small, he, lower, upper)
  expect_true(conv$converged)
  expect_equal(conv$grad_max, 1e-5)
  expect_identical(conv$at_bound, 2L)
  expect_false(kits:::shash_convergence(c(1, 5), small, he, lower, upper)$converged)     # 都自由：梯度 3
  expect_false(kits:::shash_convergence(c(1, 10), \(p) c(1e-3, 0), he, lower, upper)$converged)
  expect_false(kits:::shash_convergence(c(1, 10), small, \(p) -diag(2), lower, upper)$hessian_pd)

  # 真的沒收斂（梯度大）：照樣警告，診斷也是 warn
  set.seed(2)
  age2 <- runif(500, 3, 36)
  D2 <- 20 + 15 * log(age2) + rnorm(500, 0, 3)
  warns <- list()
  bad <- withCallingHandlers(kits_norms(D2, age2, ages = seq(4, 35, by = 0.5), degrees = c(3, 2, 2, 1), pad = 0),
                             warning = function(w) {
                               warns[[length(warns) + 1]] <<- w
                               invokeRestart("muffleWarning")
                             })
  dg <- attr(bad, "diagnostics")
  if (dg$value[dg$check == "max |gradient|"] > 1e-4) {
    expect_true(any(vapply(warns, inherits, logical(1), what = "kits_warning_convergence")) || dg$value[1] == 0)
    expect_equal(dg$status[dg$check == "max |gradient|"], "warn")
  }
})

test_that("age_pass() 的 p50、kits_item_positions、kits_ladder_marks", {
  age <- rep(c(6, 12), each = 4)
  P <- cbind(a1 = c(0, 0, 0, 1, 1, 1, 1, 0), a2 = rep(1, 8), a3 = NA)
  out <- age50_of(P, age, age)
  expect_equal(names(out), c("a1", "a2", "a3"))
  expect_equal(unname(is.na(out)), c(FALSE, TRUE, TRUE))
  expect_equal(age50_of(P * 2, age, age), out)              # scores: at least step 1 is a pass
  expect_equal(age50_of(P * 2, age, age, step = 2), out)
  expect_equal(age50_of(P == 1, age, age, step = 2), out)    # logical: used as it is
  norm <- data.frame(score = "lan", age_m = c(6, 12, 18, 24), P50 = c(40, 52, 60, 66))
  items <- data.frame(item_id = c("a1", "a2", "a3"), score = "lan", D = c(45, 61, 30), age50 = c(8.5, NA, NA))
  pos <- old_item_positions(items, norm)
  expect_equal(pos$x, c(8.5, 24, 6))
  expect_equal(pos$below_range, c(FALSE, FALSE, TRUE))
  norm2 <- data.frame(score = rep(c("lan", "cog"), each = 2), age_m = c(12, 24, 12, 24), P50 = c(50, 60, 54, 66))
  expect_equal(old_ladder_marks(norm2, scores = c("lan", "cog"), ages = c(12, 24))$D, c(52, 63))
})
