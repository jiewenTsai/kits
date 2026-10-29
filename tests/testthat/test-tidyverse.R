# 0.15.0：tibble 的輸入、dplyr 的動詞、和這一版的修正（sd_true、precision 的 criterion、spread = "se" 的警告）

betas_t <- list(a1 = c(-1, 0), a2 = c(0, 1), a3 = 0.5, a4 = -0.5, a5 = c(-0.5, 0.5))
X_t <- kits_sim(seq(-3, 3, length.out = 400), betas_t, seed = 1)
X_t[sample(length(X_t), 100)] <- NA
X_t <- X_t[rowSums(!is.na(X_t)) > 0, ]
persons_t <- data.frame(visit_id = rownames(X_t), wave_m = rep(c(12, 24), length.out = nrow(X_t)))
persons_t$age_m <- persons_t$wave_m + seq(-0.5, 0.5, length.out = nrow(X_t))
item_t <- kits_item(X_t, scale = "lan M03")

test_that("data 和 persons 可以是 tibble：寬的（visit_id 欄）、長的；結果和矩陣相同", {
  wide <- dplyr::as_tibble(X_t, rownames = "visit_id")
  long <- tidyr::pivot_longer(wide, -visit_id, names_to = "item_id", values_to = "item_score",
                              values_drop_na = TRUE)
  for (d in list(wide, long)) {
    it <- kits_item(d, scale = "lan M03")
    expect_equal(coef(it), coef(item_t), ignore_attr = TRUE)
    expect_equal(kits_subj(d, it, persons = dplyr::as_tibble(persons_t)), kits_subj(X_t, item_t, persons = persons_t),
                 ignore_attr = TRUE)
    expect_equal(kits_itemfit(d, it), kits_itemfit(X_t, item_t), ignore_attr = TRUE)
    wave <- stats::setNames(persons_t$wave_m, persons_t$visit_id)
    expect_equal(kits_subj_table(kits_subj(d, it, persons = persons_t)),
                 subj_table_of(item_t, X_t, wave_m = wave), ignore_attr = TRUE)
  }
  s <- kits_subj(wide, item_t, persons = dplyr::as_tibble(persons_t))
  expect_s3_class(s, "kits_subj")
  expect_false(inherits(s, "tbl_df"))
  sex <- rep(c("boy", "girl"), length.out = nrow(X_t))
  expect_equal(suppressWarnings(kits_dif(wide, sex))$dif, suppressWarnings(kits_dif(X_t, sex))$dif)
  expect_equal(kits_drift(wide, item_t)$items, kits_drift(X_t, item_t)$items, ignore_attr = TRUE)
})

test_that("表格的輸入可以是 tibble：固定參數、門檻表、錨題、常模、題目表、垂直量尺、參數表", {
  tb <- dplyr::as_tibble
  fixed <- dplyr::filter(tb(as.data.frame(coef(item_t))), item_id %in% c("a1", "a3"))
  expect_equal(coef(kits_item(X_t, fixed_params = fixed)), coef(kits_item(X_t, fixed_params = as.data.frame(fixed))), ignore_attr = TRUE)
  anchors <- dplyr::tibble(item_id = c("a1", "a2", "a3", "a4", "a5"), tau = c(40, 47, 51, 55, 60))
  steps <- tb(as.data.frame(coef(item_t)))
  expect_equal(kits_linking(steps, anchors)$intercept,
               kits_linking(item_t, as.data.frame(anchors))$intercept)
  norm <- dplyr::tibble(age_m = c(6, 12, 18), mu = c(40, 55, 62), sigma = 5, epsilon = 0, delta = 1, P50 = c(40, 55, 62),
                         score = "lan")
  expect_equal(kits_norm_z(c(40, 50), c(6, 9), norm), kits_norm_z(c(40, 50), c(6, 9), as.data.frame(norm)))
  items <- dplyr::tibble(item_id = c("a1", "a2"), score = "lan", D = c(45, 61), age50 = c(8, NA))
  expect_equal(as.data.frame(old_item_positions(items, norm)),
               as.data.frame(old_item_positions(as.data.frame(items), norm)))
  expect_s3_class(old_ladder_marks(norm, "lan", ages = c(6, 12)), "kits_table")
  d1 <- dplyr::tibble(item_id = c("a1", "a2"), dif = c(0.5, 0.1))
  expect_equal(kits_dif_combine(d1, d1)$DIF, c("favors girls", "none"))
  v <- dplyr::tibble(wave_m = rep(c(12, 24), each = 50), theta = rep(c(0, 1), each = 50) + rep(seq(-1, 1, length.out = 50), 2), se = 0.3)
  expect_equal(kits_vertical(v)$sd_true, kits_vertical(as.data.frame(v))$sd_true)
  expect_s3_class(kits_subj_table(tb(as.data.frame(coef(item_t)))), "kits_table")
})

test_that("dplyr 的動詞：kits 的表保留類別和量尺，拿掉描述整張表的診斷和來源", {
  s <- kits_subj(X_t, item_t, persons = persons_t)
  ok <- dplyr::filter(s, flag == "ok")
  expect_s3_class(ok, "kits_subj")
  expect_equal(attr(ok, "scale"), "lan M03")
  expect_null(attr(ok, "diagnostics"))
  expect_null(attr(ok, "info"))
  expect_equal(nrow(ok), sum(s$flag == "ok"))
  m <- dplyr::mutate(s, z = theta * 2)
  expect_equal(attr(m, "scale"), "lan M03")
  expect_s3_class(dplyr::arrange(coef(item_t), thr), "kits_table")
  fix <- dplyr::filter(coef(item_t), item_id %in% c("a1", "a2"))
  expect_equal(attr(fix, "scale"), "lan M03")
  expect_equal(kits_item(X_t, fixed_params = fix)$scale, "lan M03")   # 固定參數的量尺跟著走
  expect_output(print(ok), "scale: lan M03")
})

test_that("kits_vertical()：sd_true 預設用所有人次；extremes = FALSE 不用 0 分和滿分的人次", {
  s <- kits_subj(X_t, item_t, persons = persons_t)
  v <- kits_vertical(s, extremes = FALSE)
  g <- s[s$wave_m == 12, ]
  inner <- g[!g$flag %in% c("zero", "perfect"), ]
  expect_equal(v$sd[1], sd(g$theta))
  expect_equal(v$sd_true[1], sqrt(max(var(inner$theta) - mean(inner$se^2), 0)))
  # 沒有 flag 欄：所有人次
  plain <- as.data.frame(s)[c("theta", "se", "wave_m")]
  expect_equal(kits_vertical(plain)$sd_true[1], sqrt(max(var(g$theta) - mean(g$se^2), 0)))
  # 預設 extremes = TRUE：所有人次（0.14.0 的算法）
  expect_equal(kits_vertical(s)$sd_true[1], sqrt(max(var(g$theta) - mean(g$se^2), 0)))
  expect_equal(kits_vertical(s, extremes = TRUE)$sd, v$sd)
  expect_error(kits_vertical(s, extremes = NA), class = "kits_error_input")
})

test_that("kits_precision()：預設的 criterion 是 1 logit；短的題本也有可用的範圍", {
  pr <- kits_precision(item_t, list(f = c("a1", "a2", "a3", "a4", "a5")))
  expect_equal(attr(pr, "criterion"), 1)
  r <- attr(pr, "ranges")
  expect_true(r$min_rmse > 0.5)       # 5 題：0.5 logit 達不到
  expect_false(is.na(r$from))
  expect_equal(kits_precision(item_t, list(f = c("a1", "a2")), criterion = 0.5) |> attr("criterion"), 0.5)
})

test_that("篩點留下不到一半的點時，診斷警告；spread = \"se\" 時說明它只適合小樣本", {
  set.seed(2)
  from <- seq(-2, 2, length.out = 12)
  to <- from + rnorm(12, 0, 0.3)
  l <- link_of_points(from, to, id = paste0("p", 1:12), spread = "se", se_from = rep(0.03, 12),
                      se_to = rep(0.03, 12))
  row <- l$diagnostics[l$diagnostics$check == "screening keeps fewer than half", ]
  expect_equal(row$status, "warn")
  expect_match(row$message, "spread = \"sd\"")
  plain <- link_of_points(from, to, id = paste0("p", 1:12))
  expect_false("screening keeps fewer than half" %in% plain$diagnostics$check)
})
