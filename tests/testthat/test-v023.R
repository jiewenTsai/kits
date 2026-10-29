# kits 0.23.0：連結的 a、b 改名 intercept、slope（數字完全相同；link$a 停下來說新的名字）；編號欄改名（visit_id、
# subj_id、pattern_id、raw、raw_max），舊的欄名這一版照常認得。

set.seed(1)
betas <- list(a1 = -2, a2 = -1, a3 = c(-0.5, 0.5), a4 = 0, a5 = c(0.5, 1.5), a6 = 2, a7 = 3)
X03 <- kits_sim(rnorm(400, 0), betas[1:5], seed = 1)
X36 <- kits_sim(rnorm(400, 2), betas[3:7], prefix = "q", seed = 2)
item03 <- kits_item(X03, scale = "lan M03")
item36 <- kits_item(X36, scale = "lan M36")
anchors <- data.frame(item_id = c("a1", "a2", "a3", "a4", "a5"), tau = c(40, 47, 51, 55, 60))
co <- kits_linking(item36, item03, n_boot = 50)
d_link <- kits_linking(item03, anchors, n_boot = 50)
s36 <- kits_convert(kits_subj(X36, item36), co, d_link)
old <- readRDS(test_path("fixtures", "ab022.rds"))   # kits 0.22.0 算的，同樣的資料和亂數

# ---- a、b 改名 intercept、slope ----------------------------------------------------------------------

test_that("intercept、slope：0.22.0 的 a（截距）、b（斜率），數字完全相同；誤差、共變異數、重抽跟著改名", {
  for (n in c("co", "d_link")) {
    new <- get(n)
    was <- old[[n]]
    expect_identical(new$intercept, was$a, info = n)
    expect_identical(new$slope, was$b, info = n)
    expect_identical(new$intercept_error, was$a_error, info = n)
    expect_identical(new$slope_error, was$b_error, info = n)
    expect_identical(new$vcov[["intercept", "slope"]], was$ab_cov, info = n)
    expect_identical(new$vcov[["slope", "intercept"]], was$ab_cov, info = n)
    expect_identical(new$intercept0, was$a0, info = n)
    expect_identical(new$slope0, was$b0, info = n)
    expect_identical(unname(new$draws[, "intercept"]), unname(was$draws[, "a"]), info = n)
    expect_identical(unname(new$draws[, "slope"]), unname(was$draws[, "b"]), info = n)
  }
  expect_identical(colnames(co$draws), c("intercept", "slope"))
  expect_identical(dimnames(co$vcov), list(c("intercept", "slope"), c("intercept", "slope")))
  expect_identical(coef(co), c(intercept = co$intercept, slope = co$slope))
  expect_identical(confint(co)$term, c("intercept", "slope"))
})

test_that("沒有元素以 a 或 b 開頭，所以 link$a、link$b 不會部分比對到別的元素（0.24.0 拿掉舊名字的檢查）", {
  expect_false(any(substr(names(unclass(co)), 1, 1) %in% c("a", "b")))
  expect_null(co$a)
  expect_null(co$b)
  expect_identical(co[["slope"]], co$slope)
  expect_identical(co$points, unclass(co)$points)
})

test_that("換算後的分數和 0.22.0 相同：new = intercept + slope x old，SE 乘 slope", {
  got <- as.data.frame(s36)[c("theta", "se", "se_link", "se_total")]
  # 0.23.1 起 WLE 由 RTMB 的牛頓法解（0.22.0 是 base R 的封閉式），theta 差約 3e-12；D 量尺上再乘斜率（約 7）。
  # 連結和換算的算式不變
  # 0.25.0 起人的 SE 是 1 / sqrt(I)（0.22.0 是 Warm 的），所以 se 和 se_total 不比；se_total 的算式不變
  for (col in c("theta", "se_link")) expect_lt(max(abs(got[[col]] - old$s36[[col]])), 1e-10)
  expect_equal(got$se_total, sqrt(got$se^2 + got$se_link^2))
  raw <- kits_subj(X36, item36)
  expect_identical(kits_convert(raw, co)$theta, co$intercept + co$slope * raw$theta)
  expect_identical(kits_convert(raw, co)$se, co$slope * raw$se)
  expect_output(print(co), "new value = intercept + slope x old value", fixed = TRUE)
  expect_output(print(co), "intercept = .*linking error")
  expect_output(print(co), "slope     = ")
  expect_identical(as.data.frame(co)$intercept, co$intercept)
})

# ---- 編號欄改名 ---------------------------------------------------------------------------------------

test_that("分數表的欄：visit_id、pattern_id、raw、raw_max；persons 的 subj_id 留下", {
  persons <- data.frame(visit_id = rownames(X03), subj_id = sub("_.*", "", rownames(X03)),
                        wave_m = rep(c(6, 12), 200))
  s <- kits_subj(X03, item03, persons = persons)
  expect_true(all(c("pattern_id", "visit_id", "raw", "raw_max", "subj_id") %in% names(s)))
  expect_false(any(c("person_id", "pid", "booklet_id", "booklet_score",
                     "booklet_max") %in% names(s)))
  expect_identical(names(kits_subj_table(item03))[1:2], c("pattern_id", "raw"))
})

test_that("舊的欄名 person_id、pid 這一版照常認得，並說一聲（kits_message_data）", {
  wide <- data.frame(person_id = rownames(X03), X03, row.names = NULL)
  expect_message(item_old <- kits_item(wide, scale = "lan M03"), "person_id",
                 class = "kits_message_data")
  expect_identical(coef(item_old), coef(item03))
  persons <- data.frame(person_id = rownames(X03), pid = rownames(X03), wave_m = 6)
  expect_message(s <- kits_subj(X03, item03, persons = persons), "pid", class = "kits_message_data")
  expect_identical(s$subj_id, rownames(X03))
  long <- kitm03$responses |>
    dplyr::filter(domain == "lan", wave_m == 3) |>
    dplyr::rename(person_id = visit_id, pid = subj_id)
  expect_message(resp <- kits_resp(long), "visit_id, subj_id", class = "kits_message_data")
  expect_true(all(c("visit_id", "subj_id") %in% names(resp)))
})

test_that("print(kits())：下一步是分數表之後的", {
  k <- kits(X03, scale = "lan M03")
  out <- capture.output(print(k))
  expect_equal(sum(startsWith(out, "Next:")), 1)
  expect_true(any(grepl("kits_subj_table", out[startsWith(out, "Next:")])))
})
