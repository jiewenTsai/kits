# 題目表的 RP 和 DIF 的方向（0.23.0 的兩個錯）：kits_item_table() 的 D 和月齡要同一個 RP；DIF 的方向
# 在 kits_dif_combine()、kits_item_table()、plot() 用同一個 dif_direction()，標籤用真正的組和 cut。

set.seed(23)
betas <- list(a1 = c(-1, 0), a2 = 0.5, a3 = -0.5, a4 = 0, a5 = 1)
th <- seq(-3, 3, length.out = 300)
item <- suppressMessages(kits_item(kits_sim(th, betas, seed = 1)))
age <- runif(800, 3, 24)
norm <- kits_norms(-3 + 2 * log(age) + rnorm(800, 0, 0.5), age, ages = 3:24,
                   degrees = c(3, 2, 1, 1), by = data.frame(domain = rep("lan", 800)))
quiet <- function(expr) suppressWarnings(suppressMessages(expr))
milestone <- function(thr) {
  thr <- as.data.frame(thr)
  stats::setNames(thr$thr[thr$item_score == 1], thr$item_id[thr$item_score == 1])
}

# ---- RP ------------------------------------------------------------------------------------------

test_that("kits_item 給 rp = 0.8：D 是 RP80 的門檻，不是 RP50", {
  t50 <- quiet(kits_item_table(item, norm))
  t80 <- quiet(kits_item_table(item, norm, rp = 0.8))
  expect_equal(t50$D, unname(milestone(coef(item))[t50$item_id]), tolerance = 0)
  expect_equal(t80$D, unname(milestone(coef(item, rp = 0.8))[t80$item_id]), tolerance = 0)
  expect_true(all(t80$D > t50$D))
  expect_equal(attr(t80, "rp"), 0.8)
  # 表給 RP80 的門檻：和 kits_item 給的相同
  from_table <- quiet(kits_item_table(coef(item, rp = 0.8), norm, rp = 0.8))
  expect_equal(as.data.frame(from_table), as.data.frame(t80), tolerance = 0)
})

test_that("表的 RP 和 rp 不同時停下來（kits_error_input）", {
  expect_error(kits_item_table(coef(item, rp = 0.8), norm), "RP80", class = "kits_error_input")
  expect_error(kits_item_table(coef(item), norm, rp = 0.8), "RP50", class = "kits_error_input")
  expect_error(kits_item_table(list(coef(item), coef(item, rp = 0.8)), norm),
               class = "kits_error_input")
  # 沒記 RP 的表是 RP50
  plain <- data.frame(item_id = c("a1", "a2"), item_score = 1, thr = c(-1, 0))
  expect_error(kits_item_table(plain, norm, rp = 0.8), class = "kits_error_input")
  expect_identical(quiet(kits_item_table(plain, norm))$D, c(-1, 0))
  # 換尺之後 RP 跟著表走
  dl <- kits_linking(item, data.frame(item_id = c("a1", "a2", "a3", "a4", "a5"),
                                      tau = c(40, 47, 51, 55, 60)))
  d80 <- quiet(kits_convert(coef(item, rp = 0.8), dl))
  expect_error(kits_item_table(d80, norm), "RP80", class = "kits_error_input")
})

# ---- DIF 的方向 -----------------------------------------------------------------------------------

sim_dif <- function(girl, seed) {
  X <- rbind(kits_sim(th, betas, prefix = "b", seed = seed),
             kits_sim(th, modifyList(betas, girl), prefix = "g", seed = seed + 1))
  list(X = X, sex = rep(c("boy", "girl"), each = 300))
}
c1 <- sim_dif(list(a2 = -0.5, a4 = 1, a5 = 2), 1)
c2 <- sim_dif(list(a2 = -0.5, a4 = -1, a5 = 2), 3)
d1_boy <- quiet(kits_dif(c1$X, c1$sex))
d2_boy <- quiet(kits_dif(c2$X, c2$sex))
d1_girl <- quiet(kits_dif(c1$X, c1$sex, ref = "girl"))
d2_girl <- quiet(kits_dif(c2$X, c2$sex, ref = "girl"))
thr <- data.frame(item_id = d1_boy$item_id, item_score = 1, thr = seq(-2, 2, length.out = 5))

test_that("dif_direction()：大於 cut 對 focal 有利，小於 -cut 對 ref 有利；NA 留 NA", {
  x <- c(0.5, -0.5, 0.1, NA)
  expect_identical(dif_direction(x, 0.3, "girl", "boy"),
                   c("favors girls", "favors boys", "none", NA))
  expect_identical(dif_direction(x, 0.3, "boy", "girl"),
                   c("favors boys", "favors girls", "none", NA))
  expect_identical(dif_direction(x, 0.05, "F", "M"), c("favors F", "favors M", "favors F", NA))
  expect_identical(dif_direction(x, 0.6), c("none", "none", "none", NA))
})

test_that("ref = \"girl\"：符號相反，標籤仍然是對的組（三處相同）", {
  expect_equal(d1_girl$dif, -d1_boy$dif, tolerance = 1e-8)
  boy <- kits_dif_combine(d1_boy, d2_boy)
  girl <- kits_dif_combine(d1_girl, d2_girl)
  expect_identical(girl$DIF, boy$DIF)
  expect_identical(boy$DIF[match(c("a2", "a4", "a5"), boy$item_id)],
                   c("favors girls", "inconsistent", "favors boys"))
  # 題目表給一個 kits_dif()：ref 是 girl 時 a5（女生比較難）仍然 favors boys
  tab_boy <- quiet(kits_item_table(thr, norm, dif = d1_boy))
  tab_girl <- quiet(kits_item_table(thr, norm, dif = d1_girl))
  expect_identical(tab_girl$DIF, tab_boy$DIF)
  expect_identical(tab_girl$DIF[tab_girl$item_id == "a5"], "favors boys")
  # plot()：a5 在 "Favors boy" 那一組
  traces <- plotly::plotly_build(plot(d1_girl))$x$data
  a5 <- traces[[which(vapply(traces, \(t) "a5" %in% t$y, logical(1)))]]
  expect_identical(a5$name, "Favors boy (> 0.3)")
  # 兩個 kits_dif() 的 ref 不同：停下來
  expect_error(kits_dif_combine(d1_boy, d2_girl), "same `ref`", class = "kits_error_input")
})

test_that("其他的 cut：題目表的 dif_cut、kits_dif_combine() 的 cut、plot() 的 cut 判斷相同", {
  for (cut in c(0.1, 0.7)) {
    expected <- dif_direction(d1_boy$dif, cut, "girl", "boy")
    tab <- quiet(kits_item_table(thr, norm, dif = d1_boy, dif_cut = cut))
    expect_identical(tab$DIF, expected[match(tab$item_id, d1_boy$item_id)])
    only_one <- kits_dif_combine(d1_boy, d1_boy, cut = cut)
    expect_identical(only_one$DIF, expected[match(only_one$item_id, d1_boy$item_id)])
    traces <- plotly::plotly_build(plot(d1_boy, cut = cut))$x$data
    up <- unlist(lapply(traces, \(t) if (identical(t$name, sprintf("Favors girl (> %g)", cut))) t$y))
    expect_setequal(up, d1_boy$item_id[expected %in% "favors girls"])
  }
  # 0.3 和 0.7 給不同的結果（這組資料）
  expect_false(identical(quiet(kits_item_table(thr, norm, dif = d1_boy))$DIF,
                         quiet(kits_item_table(thr, norm, dif = d1_boy, dif_cut = 0.7))$DIF))
  expect_error(kits_item_table(thr, norm, dif = d1_boy, dif_cut = -1), "dif_cut",
               class = "kits_error_input")
})

test_that("kits_item_table(step = )：每題自己的通過步驟", {
  set.seed(1)
  age <- runif(800, 3, 24)
  D <- 20 + 15 * log(age) + rnorm(800, 0, 3)
  norm <- kits_norms(D, age, ages = 3:24, degrees = c(3, 2, 1, 1),
                     by = data.frame(domain = rep("lan", 800)))
  thr <- data.frame(item_id = rep(c("a1", "a5"), c(1, 4)), item_score = c(1, 1:4),
                    thr = c(40, 30, 45, 58, 70), domain = "lan")
  tab1 <- kits_item_table(thr, norm)
  expect_false("step" %in% names(tab1))
  expect_identical(tab1$D[tab1$item_id == "a5"], 30)
  tab3 <- kits_item_table(thr, norm, step = c(a5 = 3))
  expect_identical(tab3$D[tab3$item_id == "a5"], 58)
  expect_identical(tab3$step[match(c("a1", "a5"), tab3$item_id)], c(1, 3))
  expect_identical(tab3$D[tab3$item_id == "a1"], 40)
  # 通過 = 至少 step 分
  score <- sapply(list(40, c(30, 45, 58, 70)), \(b) {
    rowSums(sapply(b, \(t) rbinom(800, 1, plogis((D - t) / 3))))
  })
  colnames(score) <- c("a1", "a5")
  resp <- data.frame(visit_id = 1:800, age_m = age, wave_m = round(age), score)
  tab_data <- kits_item_table(thr, norm, data = resp, step = c(a5 = 3))
  pass <- kits:::age_pass(cbind(a1 = score[, 1] >= 1, a5 = score[, 2] >= 3), age, round(age))
  expect_identical(tab_data$age50[match(c("a1", "a5"), tab_data$item_id)], pass$p50)
  expect_error(kits_item_table(thr, norm, step = c(a9 = 2)), "not in `thr`", class = "kits_error_input")
  expect_error(kits_item_table(thr, norm, step = c(a1 = 2)), "without their pass step", class = "kits_error_input")
  expect_error(kits_item_table(thr, norm, step = c(2, 3)), "named by item", class = "kits_error_input")
})
