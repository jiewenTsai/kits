# plot() 的新方法：校準的類別曲線、性別 DIF、題目適配、通過月齡、偏移、篩波次。都是 plotly，有該有的線和點。

betas <- list(a1 = c(-1, 0), a2 = c(0, 1), a3 = 0.5, a4 = -0.5, a5 = c(-2, -1, 1), a6 = 1.5)
X <- kits_sim(rnorm(600, 0, 1.5), betas, seed = 1)
item <- kits_item(X, scale = "lan M03")
built_of <- function(p) plotly::plotly_build(p)$x
names_of <- function(p) vapply(built_of(p)$data, function(d) as.character(d$name)[1], character(1))

test_that("plot(kits_item)：每個分數一條類別機率曲線，一題一格；type = \"expected\" 一題一條線", {
  p <- plot(item)
  expect_s3_class(p, "plotly")
  expect_length(built_of(p)$data, 16)   # 3 + 3 + 2 + 2 + 4 + 2 個分數
  one <- built_of(plot(item, items = "a5"))$data
  expect_length(one, 4)
  probs <- kits:::pcm_probs(one[[2]]$x, coef(item)$beta[coef(item)$item_id == "a5"])
  expect_equal(as.vector(one[[2]]$y), probs[, 2])
  expected <- built_of(plot(item, type = "expected"))$data
  expect_length(expected, 6)
  expect_identical(names_of(plot(item, type = "expected")), paste0("a", 1:6))
  expect_error(plot(item, items = "zz"), "zz", class = "kits_error_input")
  expect_error(plot(item, type = "info"), "type", class = "kits_error_input")
  flat <- stats::setNames(rep(list(0), 12), sprintf("b%02d", 1:12))
  many <- kits_item(kits_sim(rnorm(600), flat, seed = 3), scale = "x")
  expect_message(p12 <- plot(many), "first 9", class = "kits_message_data")
  expect_length(built_of(p12)$data, 18)
})

test_that("plot(kits_item, type = \"observed\")：每組的通過率、上下短線、灰色的模型曲線", {
  grp <- rep(c("g1", "g2"), length.out = nrow(X))
  p <- plot(item, type = "observed", data = X, group = grp, items = c("a1", "a5"))
  tr <- built_of(p)$data
  expect_length(tr, 2 * (1 + 2 * 2))   # 兩格：模型 + 每組（線、短線）
  expect_identical(tr[[1]]$name, "Model")
  b5 <- coef(item)$beta[coef(item)$item_id == "a5"]
  model5 <- tr[[6]]
  P <- kits:::pcm_probs(model5$x, b5)
  expect_equal(as.vector(model5$y), 100 * rowSums(P[, -1]))   # step = 1：P(X >= 1)
  # 一段的通過率照手算：g1 的 a1，theta 分 0.5 一段，至少 10 人
  th <- kits_subj(X, item)$theta
  pass <- X[, "a1"] >= 1 & grp == "g1"
  bin <- (floor(th / 0.5) + 0.5) * 0.5
  hand <- tapply(pass[grp == "g1"], bin[grp == "g1"], mean)
  n <- tapply(pass[grp == "g1"], bin[grp == "g1"], length)
  expect_equal(as.vector(tr[[2]]$y), 100 * as.vector(hand[n >= 10]))
  # step 照題：a5 用 score >= 2，模型線跟著換
  p2 <- plot(item, type = "observed", data = X, items = "a5", step = c(a5 = 2))
  m2 <- built_of(p2)$data[[1]]
  expect_equal(as.vector(m2$y), 100 * rowSums(kits:::pcm_probs(m2$x, b5)[, 3:4]))
  expect_length(built_of(p2)$data, 3)   # 沒有 group：一條線
  expect_error(plot(item, type = "observed"), "responses", class = "kits_error_input")
  expect_error(plot(item, type = "observed", data = X, step = 0), "step", class = "kits_error_input")
  expect_error(plot(item, type = "observed", data = X, step = c(zz = 2)), "zz",
               class = "kits_error_input")
  expect_error(plot(item, type = "observed", data = X, items = "a3", step = 2), "a3",
               class = "kits_error_input")   # a3 只有一個步驟
  na_grp <- replace(grp, 1:5, NA)
  expect_message(plot(item, type = "observed", data = X, group = na_grp, items = "a1"),
                 "5 person-wave", class = "kits_message_data")
  one <- plot(item, type = "observed", data = X[, "a5", drop = FALSE], items = "a5")
  expect_length(built_of(one)$layout$annotations, 1)   # 適配算不出來：只有標題
  other <- kits_subj(X, item)
  attr(other, "scale") <- "D"
  expect_error(plot(item, type = "observed", data = X, scores = other), "scale",
               class = "kits_error_scale")
})

test_that("plot(kits_dif)：每題的 dif 和 95% 信賴區間；超過 cut 的上色", {
  theta <- seq(-3, 3, length.out = 400)
  boys <- kits_sim(theta, betas, prefix = "b", seed = 1)
  girls <- kits_sim(theta, modifyList(betas, list(a3 = 1.5)), prefix = "g", seed = 2)
  dif <- kits_dif(rbind(boys, girls), rep(c("boy", "girl"), each = 400))
  se <- attr(dif, "se")
  expect_identical(names(se), dif$item_id)
  expect_true(all(se > 0 & se < 1))
  p <- plot(dif)
  data <- built_of(p)$data
  flagged <- data[[which(grepl("Favors boy", names_of(p)))]]
  expect_identical(as.character(flagged$y), "a3")
  expect_equal(as.vector(flagged$error_x$array), 1.96 * unname(se["a3"]))
  expect_length(built_of(plot(dif, items = c("a1", "a2")))$data, 1)
  expect_error(plot(dif, cut = -1), "cut", class = "kits_error_input")
})

test_that("plot(kits_itemfit)：infit 和 outfit 每題一點，可以接受的範圍畫成灰色的帶", {
  fit <- kits_itemfit(X, item)
  b <- built_of(plot(fit))
  stat_of <- function(d) sub(":.*", "", d$name)   # 每條是一種統計量、一種 flag，例如 infit: ok
  stats <- vapply(b$data, stat_of, character(1))
  expect_identical(unique(stats), c("infit", "outfit"))
  infit_traces <- b$data[stats == "infit"]
  y <- unlist(lapply(infit_traces, \(d) d$y))
  x <- unlist(lapply(infit_traces, \(d) as.character(d$x)))
  expect_equal(y[match(fit$item_id, x)], fit$infit)
  band <- Filter(function(s) identical(s$fillcolor, "rgba(127,127,127,0.15)"), b$layout$shapes)
  expect_equal(c(band[[1]]$y0, band[[1]]$y1), c(0.6, 1.4))
  outfit_only <- built_of(plot(fit, type = "outfit"))$data
  expect_true(all(startsWith(vapply(outfit_only, \(d) d$name, character(1)), "outfit: ")))
  expect_error(plot(fit, range = c(1.4, 0.6)), "range", class = "kits_error_input")
})

test_that("plot(kits_drift)：新舊門檻和那條線；偏移的題上色；type = \"difference\" 是到線的距離", {
  new <- kits_sim(rnorm(600, 0.5, 1.5), modifyList(betas, list(a3 = 2)), prefix = "n", seed = 2)
  drift <- kits_drift(new, item)
  p <- plot(drift)
  nm <- names_of(p)
  expect_true(any(grepl("^Final line", nm)))
  drifted <- built_of(p)$data[[which(nm == "Drift (step dropped)")]]
  expect_match(drifted$hovertext[1], "^a3:1")
  diff <- built_of(plot(drift, type = "difference"))$data
  all_y <- unlist(lapply(diff, function(d) d$y))
  resid <- drift$points$resid
  expect_setequal(round(all_y[is.finite(all_y)], 10), round(resid[is.finite(resid)], 10))
  expect_error(plot(drift, items = "zz"), "zz", class = "kits_error_input")
})

test_that("plot(kits_norm_screen)：每個波次最低分的比例；拿掉的波次是橘色；cut 是虛線", {
  set.seed(5)
  age <- runif(900, 3, 30)
  wave <- round(age / 6) * 6
  D <- 20 + 15 * log(age) + rnorm(900, 0, 3)
  D[wave == 6][1:60] <- min(D[wave == 6])
  screen <- kits_norm_screen(D, wave)
  b <- built_of(plot(screen))
  expect_identical(vapply(b$data, function(d) d$name, character(1)), c("In the norms", "Left out"))
  expect_identical(as.character(b$data[[2]]$x), "6")
  lines <- Filter(function(s) identical(s$type, "line"), b$layout$shapes)
  expect_equal(lines[[1]]$y0, screen$cut)
})

test_that("plot(kits_item_table)：DIF 的顏色照方向（focal 橘、ref 藍）；男女以外的組也可以畫", {
  set.seed(1)
  age <- stats::runif(800, 3, 24)
  D <- 20 + 15 * log(age) + stats::rnorm(800, 0, 3)
  norm <- kits_norms(D, age, ages = 3:24, degrees = c(3, 2, 1, 1),
                     by = data.frame(domain = rep("lan", 800)))
  thr <- data.frame(item_id = c("a1", "a2", "a3", "a4"), item_score = 1, thr = c(40, 55, 65, 60),
                    domain = "lan")
  pass <- pass_resp(c(a1 = 4.5, a2 = 11, a3 = 20, a4 = 15), age)
  # 一個 kits_dif()（ref = "M"）：dif > 0 對 F 有利
  dif <- data.frame(item_id = c("a1", "a2", "a3", "a4"), dif = c(0, 0.5, -0.6, 0.1))
  attr(dif, "groups") <- list(ref = "M", focal = "F", n = c(M = 400, F = 400))
  tab <- kits_item_table(thr, norm, data = pass, dif = dif)
  expect_identical(tab$DIF[match(c("a2", "a3"), tab$item_id)], c("favors F", "favors M"))
  expect_true(all(tab$empirical))   # 每題都畫在觀察到的 50% 通過月齡
  expect_identical(tab$x, tab$age50)
  colour_of <- function(p, name) {
    tr <- Filter(\(d) identical(d$name, name), plotly::plotly_build(p)$x$data)
    tr[[1]]$marker$color
  }
  p <- plot(tab)
  expect_identical(colour_of(p, "Favors F"), "#D55E00")
  expect_identical(colour_of(p, "Favors M"), "#0072B2")
  expect_s3_class(plot(tab, type = "ladder", ages = c(6, 12, 18)), "plotly")
  # 男女：和以前相同的顏色（ref = "girl" 時也一樣）
  attr(dif, "groups") <- list(ref = "girl", focal = "boy", n = c(girl = 400, boy = 400))
  tab_g <- kits_item_table(thr, norm, data = pass, dif = dif)
  p <- plot(tab_g)
  expect_identical(colour_of(p, "Favors girls"), "#D55E00")
  expect_identical(colour_of(p, "Favors boys"), "#0072B2")
})
