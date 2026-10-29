# 0.18.0：連結和常模的分步函數、as.data.frame()、專家審查 Part B（B1 重複的 id、B2 CML 的速度、B4 完全分離、
# B5 參數改名和 S3、B6 風格）。

betas7 <- list(a1 = -2, a2 = -1, a3 = c(-0.5, 0.5), a4 = 0, a5 = c(0.5, 1.5), a6 = 2, a7 = 3)
young <- kits_sim(rnorm(500, 0), betas7[1:5], prefix = "y", seed = 1)
old <- kits_sim(rnorm(500, 2), betas7[3:7], prefix = "o", seed = 2)
item_y <- kits_item(young, scale = "young")
item_o <- kits_item(old, scale = "old")
anchors <- data.frame(item_id = paste0("a", 1:5), tau = c(40, 47, 51, 55, 60))
set.seed(5)
age <- runif(1500, 3, 36)
D <- 20 + 15 * log(age) + rnorm(1500, 0, 3)
persons <- data.frame(visit_id = rownames(old), age_m = runif(500, 30, 40), wave_m = 36)
s_old <- kits_subj(old, item_o, persons = persons)

# ---- 連結的四步 -------------------------------------------------------------------------------------------

test_that("連結的四步（points -> screen -> line -> error）和 kits_linking() 完全相同", {
  for (args in list(list(), list(step = 1), list(spread = "mad", unit = "item"), list(spread = "se", cut = 10))) {
    link <- do.call(kits_linking, c(list(item_o, item_y, n_boot = 50), args))
    points <- kits_link_points(item_o, item_y, step = args$step)
    screened <- kits_link_screen(points, cut = link$settings$cut, spread = link$settings$spread, unit = link$settings$unit)
    steps <- kits_link_error(kits_link_line(screened), n_boot = 50, seed = 1)
    for (n in setdiff(names(link), "info")) expect_identical(steps[[n]], link[[n]], info = n)
    expect_identical(names(steps), names(link))
  }
  # anchors with domains and the dscore key
  steps_tab <- data.frame(item_id = paste0("a", 1:12), item_score = 1, thr = seq(-2, 2, length.out = 12),
                          domain = rep(c("lan", "cog"), each = 6))
  to <- data.frame(item_id = steps_tab$item_id, tau = 60 + 7 * steps_tab$thr + c(0.3, -0.2), key = "gsed2406")
  link <- kits_linking(steps_tab, to, n_boot = 40)
  steps <- kits_link_points(steps_tab, to) |>
    kits_link_screen() |>
    kits_link_line() |>
    kits_link_error(n_boot = 40)
  for (n in setdiff(names(link), "info")) expect_identical(steps[[n]], link[[n]], info = n)
  expect_identical(steps$settings$key, "gsed2406")
  expect_s3_class(kits_link_points(steps_tab, to), "kits_link_points")
})

test_that("連結的第 2-4 步：點也可以是自己的 data.frame；n_boot = 0 拿掉連結誤差", {
  from <- c(-2, -1, 0, 1, 2, 3, 0.5, 1.5)
  to <- 1 + from + c(0, 0, 0, 0, 0, 0, 0, 3)
  points <- data.frame(id = paste0("p", 1:8), item = paste0("p", 1:8), from = from, to = to)
  ms <- kits_link_error(kits_link_line(kits_link_screen(points)), n_boot = 30)
  none <- kits_link_error(ms, n_boot = 0)
  expect_true(is.na(none$intercept_error))
  expect_null(none$draws)
  expect_identical(none$intercept, ms$intercept)
  expect_error(kits_link_line(points), class = "kits_error_input")
  expect_error(kits_link_error(points), class = "kits_error_input")
  expect_error(kits_link_screen(points[c("id", "from")]), class = "kits_error_input")
  expect_error(kits_link_error(ms, n_boot = 1), class = "kits_error_input")
  expect_output(print(kits_link_screen(points)), "resid")
})

# ---- 常模的三步、S3 -----------------------------------------------------------------------------------------

test_that("常模的三步（data -> fit -> table）和 kits_norms() 相同；kits_norms(s, ages) 用 S3 派送", {
  norm <- suppressWarnings(kits_norms(D, age, ages = 6:30, degrees = c(3, 2, 1, 1), pad = 0))
  norm_data <- kits_norm_data(D, age, ages = 6:30, pad = 0)
  fit <- suppressWarnings(kits_norm_fit(norm_data, degrees = c(3, 2, 1, 1)))
  steps <- suppressWarnings(kits_norm_table(fit))
  no_info <- function(x) {
    attr(x, "info") <- NULL
    as.data.frame(x)
  }
  expect_identical(no_info(steps), no_info(norm))
  expect_identical(attr(steps, "kits_args"), attr(norm, "kits_args"))
  expect_identical(attr(steps, "diagnostics"), attr(norm, "diagnostics"))
  expect_identical(kits_info(steps)$data, kits_info(norm)$data)
  expect_equal(kits_info(norm)$call,
               "kits_norms(D, age, ages = 6:30, degrees = c(3, 2, 1, 1), pad = 0)")
  expect_s3_class(norm_data, "kits_norm_data")
  expect_equal(names(norm_data), c("d", "age_m"))
  expect_output(print(norm_data), "kits_norm_data")
  expect_output(print(fit), "kits_norm_fit")
  expect_error(kits_norm_fit(data.frame(d = 1, age_m = 2)), class = "kits_error_input")
  expect_error(kits_norm_table(norm_data), class = "kits_error_input")
  expect_error(kits_norm_data(s_old, ages = 30:40), "vectors", class = "kits_error_input")
})

test_that("B5 kits_norms(s, ages) 和 kits_norm_z(s, norm)：S3 派送；s 給了 age_m 時清楚地停下來（不再移到 ages）", {
  n1 <- suppressMessages(suppressWarnings(kits_norms(s_old, 30:40, degrees = c(1, 1, 1, 1))))
  n2 <- suppressMessages(suppressWarnings(kits_norms(s_old$theta, s_old$age_m, 30:40, degrees = c(1, 1, 1, 1))))
  expect_identical(n1$P50, n2$P50)
  expect_equal(attr(n1, "scale"), "old")
  expect_error(kits_norms(s_old, age_m = 30:40), "age_m", class = "kits_error_input")
  expect_error(kits_norms(s_old, ages = 30:40, age_m = s_old$age_m), "age_m", class = "kits_error_input")
  # ages：kits_norm_ages()
  auto <- suppressMessages(suppressWarnings(kits_norms(s_old, degrees = c(1, 1, 1, 1))))
  expect_identical(kits_info(auto)$settings$ages_from, "kits_norm_ages()")
  z1 <- kits_norm_z(s_old, n1)
  expect_identical(z1$z, kits_norm_z(s_old$theta, s_old$age_m, n1))
  expect_identical(z1$pct, 100 * pnorm(z1$z))
  expect_error(kits_norm_z(s_old, age_m = s_old$age_m, norm = n1), "age_m", class = "kits_error_input")
  expect_error(kits_norm_z(s_old), "norm", class = "kits_error_input")
  expect_error(kits_norms(D, age, ages = 6:30, degree = 3), "degree", class = "kits_error_input")   # a misspelled name
})

# ---- B5 改名：新的名字、舊的名字（警告、結果相同） ----------------------------------------------------------------

test_that("B5 dataSrc -> data（0.23.0 拿掉舊的名字）：舊的名字停下來，說它不是這個函數的參數", {
  expect_identical(coef(kits_item(data = young, scale = "young")), coef(item_y))
  expect_error(kits_item(dataSrc = young, scale = "young"), "data")
  expect_error(kits_subj(young, item_y, NULL, "extra"), "does not have", class = "kits_error_input")
  expect_error(kits_dif(young, rep(c("boy", "girl"), 250), wave = 1), "wave", class = "kits_error_input")
})

test_that("B5 D -> d（0.23.0 拿掉舊的名字）：收 D 時停下來", {
  norm <- suppressWarnings(kits_norms(d = D, age_m = age, ages = 6:30, degrees = c(2, 1, 1, 1)))
  wave <- round(age / 6) * 6
  expect_error(kits_norm_z(D = D[1:5], age_m = age[1:5], norm = norm), "D",
               class = "kits_error_input")
  expect_error(kits_norm_screen(d = D, wave_m = wave, D = D), "D", class = "kits_error_input")
  expect_error(kits_norm_data(d = D, age_m = age, ages = 6:30, D = D), "D",
               class = "kits_error_input")
  expect_error(plot(norm, D = D, age_m = age), "does not have", class = "kits_error_input")
  expect_error(kits_norms(D, age, ages = 6:30, de = 2), "multiple")   # R's own error, as in 0.17.0
  expect_no_warning(suppressMessages(kits_norm_screen(D, wave)))
})

# ---- as.data.frame() -----------------------------------------------------------------------------------------

test_that("as.data.frame()：連結一列（要求的欄都在），幾個連結 bind_rows() 疊起來", {
  link <- kits_linking(item_o, item_y, n_boot = 50)
  d_link <- kits_linking(item_y, dplyr::mutate(anchors, key = "gsed2406"), n_boot = 20)
  ms <- link_of_points(c(-2, -1, 0, 1, 2), c(-1, 0, 1, 2, 3.2), id = paste0("p", 1:5))
  row <- as.data.frame(link)
  need <- c("from", "to", "intercept", "slope", "intercept_error", "slope_error",
            "cov_intercept_slope", "n_points", "n_dropped", "cut", "spread", "unit",
            "step", "n_boot", "seed", "key")
  expect_true(all(need %in% names(row)))
  expect_equal(nrow(row), 1)
  expect_equal(row$intercept, link$intercept)
  expect_equal(row$cov_intercept_slope, link$vcov[1, 2])
  expect_equal(row$from, "old")
  expect_equal(row$to, "young")
  expect_equal(row$step, "all")
  expect_equal(row$n_dropped, length(link$dropped))
  all_links <- dplyr::bind_rows(row, as.data.frame(d_link), as.data.frame(ms))
  expect_equal(nrow(all_links), 3)
  expect_equal(all_links$key, c(NA, "gsed2406", NA))
  expect_equal(all_links$kind, c("calibration", "anchors", NA))
  expect_true(is.na(all_links$cov_intercept_slope[3]))
})

test_that("as.data.frame()：drift、link_check、info、item、kits、norm_screen、norm_fit、stats", {
  drift <- kits_drift(old, item_y)
  dd <- as.data.frame(drift)
  expect_equal(nrow(dd), 1)
  expect_equal(dd$intercept, drift$intercept)
  expect_equal(dd$n_fix, sum(drift$items$suggest_fix))
  info <- as.data.frame(kits_info(item_y))
  expect_equal(info$call, "kits_item(young, scale = \"young\")")
  expect_equal(info$data_md5, kits_info(item_y)$data$scores_md5)
  expect_equal(nrow(dplyr::bind_rows(info, as.data.frame(kits_info(item_o)))), 2)
  di <- as.data.frame(item_y)
  expect_equal(di$n_persons, item_y$n_persons)
  expect_equal(di$neg_loglik, item_y$fit$est$nll)
  expect_equal(names(as.data.frame(kits(young, scale = "y"))),
               c(names(di), "theta_mean", "theta_sd", "n_zero", "n_perfect", "n_misfit"))
  screen <- kits_norm_screen(D, round(age / 6) * 6)
  expect_identical(as.data.frame(screen)$floor_pct, screen$waves$floor_pct)
  fit <- kits_norm_fit(kits_norm_data(D, age, ages = 6:30), degrees = c(2, 1, 1, 1))
  expect_equal(as.data.frame(fit)$neg_loglik, fit$model$value)
  bk <- as.data.frame(item_y$stats)
  expect_equal(sum(bk$n), item_y$n_persons)
})

# ---- B1 kits_link_check() 的 person id 重複 -------------------------------------------------------------------

test_that("B1 kits_link_check()：兩個世代的 person id 重複時，用 persons_from、persons_to；結果和不重複時相同", {
  betas <- lapply(seq(-3, 3, length.out = 10), function(b) b + c(-0.5, 0.5))
  names(betas) <- paste0("a", 1:10)
  set.seed(1)
  wave <- rep(c(24, 36, 36, 48), each = 200)
  theta <- 2 * log(wave / 24) + rnorm(800)
  yy <- kits_sim(theta[1:400], betas[1:7], prefix = "y", seed = 1)
  oo <- kits_sim(theta[401:800], betas[4:10], prefix = "o", seed = 2)
  item_yy <- kits_item(yy, scale = "young")
  item_oo <- kits_item(oo, scale = "old")
  link <- kits_linking(item_oo, item_yy, n_boot = 30)
  p_all <- data.frame(visit_id = c(rownames(yy), rownames(oo)), wave_m = wave)
  chk <- kits_link_check(link, item_oo, item_yy, oo, yy, persons = p_all)
  # the same children, the same ids in both cohorts
  oo2 <- oo
  rownames(oo2) <- rownames(yy)
  p_old <- data.frame(visit_id = rownames(oo2), wave_m = wave[401:800])
  p_young <- data.frame(visit_id = rownames(yy), wave_m = wave[1:400])
  chk2 <- kits_link_check(link, item_oo, item_yy, oo2, yy, persons_from = p_old, persons_to = p_young)
  expect_equal(chk2$methods$intercept, chk$methods$intercept, tolerance = 1e-10)
  expect_equal(chk2$waves$mean_from, chk$waves$mean_from, tolerance = 1e-10)
  expect_equal(chk2$waves$wave_m, 36)
  overlapping <- dplyr::bind_rows(p_old, p_young)
  expect_error(kits_link_check(link, item_oo, item_yy, oo2, yy, persons = overlapping), "persons_from",
               class = "kits_error_input")
  expect_error(kits_link_check(link, item_oo, item_yy, oo2, yy, persons_from = p_old), "both", class = "kits_error_input")
  expect_error(kits_link_check(link, item_oo, item_yy, oo2, yy, persons = p_all, persons_to = p_young), "not both",
               class = "kits_error_input")
})

# ---- B2 CML：所有題本一起算的概似和一個題本一個迴圈的相同 ------------------------------------------------------

test_that("B2 CML：所有題本一起算的 -log L 和梯度，和一個題本一個迴圈（教學程式 cml.R 的寫法）相同", {
  set.seed(2)
  X <- kits_sim(rnorm(800), betas7, seed = 3)
  X[sample(length(X), 1200)] <- NA
  X <- X[rowSums(!is.na(X)) > 1, ]
  stats <- kits_stats(X)
  inp <- kits:::cml_inputs(stats)
  con <- kits:::cml_constraint(inp, NULL)
  poly_conv <- function(g, e) {   # 兩個多項式的乘積（係數的卷積）；e 可以是 RTMB 的 advector
    "[<-" <- RTMB::ADoverload("[<-")
    out <- numeric(length(g) + length(e) - 1)
    for (k in seq_along(e)) {
      idx <- seq_along(g) + k - 1
      out[idx] <- out[idx] + g * e[k]
    }
    out
  }
  loop_nll <- function(par) {   # 0.17.0 的寫法：一個題本一個迴圈
    beta <- con$beta0 + drop(con$A %*% par$b)
    steps_of <- split(seq_along(inp$step_item), factor(inp$step_item, levels = inp$items))
    loglik <- -sum(beta * inp$n_ge)
    for (b in seq_along(inp$book_items)) {
      its <- inp$book_items[[b]]
      centre <- mean(beta[unlist(steps_of[its])])
      g <- 1
      for (i in its) {
        e <- exp(-cumsum(beta[steps_of[[i]]] - centre))
        g <- poly_conv(g, c(1 + 0 * e[1], e))   # (1, exp(-beta_1), exp(-beta_1 - beta_2), ...)
      }
      loglik <- loglik - sum(inp$n_raw[[b]] * (log(g) - (seq_along(g) - 1) * centre))
    }
    -loglik
  }
  new <- RTMB::MakeADFun(kits:::cml_nll(inp, con), list(b = numeric(ncol(con$A))), silent = TRUE)
  ref <- RTMB::MakeADFun(loop_nll, list(b = numeric(ncol(con$A))), silent = TRUE)
  par <- rnorm(ncol(con$A), 0, 0.5)
  expect_lt(abs(new$fn(par) - ref$fn(par)) / abs(ref$fn(par)), 1e-13)
  expect_lt(max(abs(new$gr(par) - ref$gr(par))), 1e-9)
  # the estimates: the same maximum as the per-booklet likelihood with nlminb and its Hessian
  fit <- kits_cml(stats)
  opt <- nlminb(ref$par, ref$fn, ref$gr, ref$he, control = list(rel.tol = 1e-14, eval.max = 1000, iter.max = 1000))
  expect_lt(max(abs(fit$est$beta - (con$beta0 + drop(con$A %*% opt$par)))), 1e-7)
  expect_lt(fit$optim$max_gradient, 1e-6)
})

# ---- B4 完全分離 -----------------------------------------------------------------------------------------------

test_that("B4 age_pass()：完全分離的題有標記和月齡的範圍；月齡不變", {
  age <- rep(c(6, 12, 18), each = 20)
  set.seed(1)
  P <- cbind(sep = rep(c(0, 1, 1), each = 20),                              # every 6-month child fails, every older passes
             ok = rbinom(60, 1, plogis((age - 12) / 3)))
  ap <- kits:::age_pass(P, age, age)
  expect_equal(ap$separated, c(TRUE, FALSE))
  expect_equal(ap$age_last_fail, c(6, max(age[P[, "ok"] == 0])))
  expect_equal(ap$age_first_pass, c(12, min(age[P[, "ok"] == 1])))
  fit <- suppressWarnings(glm(P[, "sep"] ~ age, family = binomial))
  expect_equal(ap$p50[1], unname(-coef(fit)[1] / coef(fit)[2]))     # the same number as before, flagged
  expect_equal(kits_check(ap, all = TRUE)$value, 1)
  expect_equal(kits_check(ap)$status, "note")
  # 題目表 data = ：同一個診斷
  norm <- data.frame(domain = NA_character_, age_m = 6:18, P10 = -4, P50 = 0, P90 = 4)
  thr <- data.frame(item_id = c("sep", "ok"), item_score = 1, thr = 0)
  tab <- kits_item_table(thr, norm, data = data.frame(visit_id = 1:60, age_m = age, wave_m = age, P))
  expect_equal(kits_check(tab)$status, "note")
  expect_identical(tab$age50, ap$p50)
  expect_no_warning(kits:::age_pass(P[, "ok", drop = FALSE], age, age))
})
