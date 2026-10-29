# PCM 的單一核心（pcm_core.R）和教學用的 kits_pcm_prob()

betas <- list(a1 = c(-1, 0.5), a2 = 0.3, a3 = c(-2, 0, 1.5))
parms <- data.frame(item_id = rep(names(betas), lengths(betas)),
                    item_score = unlist(lapply(betas, seq_along)), beta = unlist(betas))
theta <- seq(-5, 5, by = 0.25)

test_that("kits_pcm_prob()：每題、每個 theta 的機率加起來是 1；長表的欄位", {
  pp <- kits_pcm_prob(parms, theta)
  expect_s3_class(pp, "kits_pcm_prob")
  expect_named(pp$prob, c("item_id", "theta", "x", "prob"))
  expect_named(pp$moments, c("item_id", "theta", "expected", "variance", "third"))
  expect_named(pp$test, c("theta", "expected", "information", "information_slope"))
  sums <- tapply(pp$prob$prob, list(pp$prob$item_id, pp$prob$theta), sum)
  expect_equal(as.vector(sums), rep(1, length(sums)), tolerance = 1e-14)
  expect_true(all(pp$prob$prob >= 0))
  expect_identical(nrow(pp$test), length(theta))
  expect_identical(as.data.frame(pp)$information, pp$test$information)
})

test_that("kits_pcm_prob()：E 和訊息和 kits_expected_score()、kits_information() 完全相同", {
  pp <- kits_pcm_prob(parms, theta)
  expect_identical(pp$test$expected, kits_expected_score(parms)(theta))
  expect_identical(pp$test$information, kits_information(parms)(theta))
  two <- kits_pcm_prob(parms, theta, items = c("a3", "a1"))
  expect_identical(two$items, c("a1", "a3"))
  expect_identical(two$test$information, kits_information(parms, items = c("a1", "a3"))(theta))
  # 每題的 E、Var 加起來就是測驗的（同一個核心）
  expect_equal(as.vector(tapply(pp$moments$variance, pp$moments$theta, sum)), pp$test$information,
               tolerance = 1e-14)
})

test_that("kits_pcm_prob()：兩個步驟的題，和逐項列舉相同；I' 是訊息的斜率", {
  b <- c(-0.7, 1.2)
  one <- data.frame(item_id = "q", item_score = 1:2, beta = b)
  th <- c(-3, -0.5, 0, 0.8, 4)
  pp <- kits_pcm_prob(one, th)
  num <- cbind(1, exp(th - b[1]), exp(2 * th - b[1] - b[2]))
  brute <- num / rowSums(num)
  got <- matrix(pp$prob$prob, ncol = 3, byrow = TRUE)
  expect_equal(got, brute, tolerance = 1e-14)
  E <- drop(brute %*% 0:2)
  V <- drop(brute %*% (0:2)^2) - E^2
  expect_equal(pp$moments$expected, E, tolerance = 1e-14)
  expect_equal(pp$moments$variance, V, tolerance = 1e-12)
  # dI/dtheta = 三階中央動差（數值微分）
  h <- 1e-5
  slope <- (kits_information(one)(th + h) - kits_information(one)(th - h)) / (2 * h)
  expect_equal(pp$test$information_slope, slope, tolerance = 1e-7)
})

test_that("algo_pcm()：極端的 theta 變異數不會是負的；RTMB 的 I' 等於三階中央動差的和", {
  m <- algo_pcm(c(-60, -30, 0, 30, 60), c(-1, 0, 2), third = TRUE)
  expect_true(all(m$Var >= 0))
  expect_true(all(is.finite(m$M3)))
  expect_equal(rowSums(m$P), rep(1, 5), tolerance = 1e-15)
  # I'（RTMB 對 I 微分）和 sum_i M3_i（dVar/dtheta = M3 的恆等式）相同
  parms2 <- data.frame(item_id = rep(c("a", "b"), c(3, 2)), item_score = c(1:3, 1:2),
                       beta = c(-1, 0, 2, 0.5, 1.5))
  th <- c(-4, -1, 0, 1.3, 5)
  m3 <- algo_pcm(th, c(-1, 0, 2), third = TRUE)$M3 + algo_pcm(th, c(0.5, 1.5), third = TRUE)$M3
  expect_equal(information_slope(th, parms2, c("a", "b")), m3, tolerance = 1e-12)
  # 門檻：P(X >= k) = 0.5 的地方
  thr <- step_threshold(c(-1, 0, 2), k = 2)[["thr"]]
  expect_equal(rowSums(pcm_probs(thr, c(-1, 0, 2))[, 3:4, drop = FALSE]), 0.5, tolerance = 1e-9)
})

test_that("kits_pcm_prob()：印出、圖、錯誤", {
  pp <- kits_pcm_prob(parms, theta)
  out <- capture.output(print(pp))
  expect_true(any(grepl("P(X = x | theta) of the PCM; information = Var(X | theta)", out, fixed = TRUE)))
  expect_true(any(grepl("kits_esf() is the other way round", out, fixed = TRUE)))
  expect_true(any(grepl("^Next:", out)))
  expect_s3_class(plot(pp), "plotly")
  expect_s3_class(plot(pp, items = "a3"), "plotly")
  expect_s3_class(plot(pp, type = "expected"), "plotly")
  expect_s3_class(plot(pp, type = "information"), "plotly")
  expect_error(plot(pp, items = "zz"), class = "kits_error_input")
  expect_error(kits_pcm_prob(parms, c(0, NA)), class = "kits_error_input")
  expect_error(kits_pcm_prob(parms), class = "kits_error_input")
  expect_error(kits_pcm_prob(parms, 0, items = "zz"), class = "kits_error_items")
})

test_that("kits_pcm_prob()：可以分享，不帶人次或孩子的編號", {
  resp <- kitm03$responses[kitm03$responses$domain == "lan", ]
  item <- suppressMessages(kits_item(kits_resp(resp), scale = "lan M03"))
  pp <- kits_pcm_prob(item, seq(-4, 4, by = 0.5))
  s <- rawToChar(serialize(pp, NULL, ascii = TRUE))
  hits <- function(ids) sum(vapply(ids, \(i) grepl(i, s, fixed = TRUE), logical(1)))
  expect_identical(hits(unique(resp$visit_id)), 0L)
  expect_identical(hits(unique(resp$subj_id)), 0L)
  expect_identical(pp$scale, "lan M03")
  expect_identical(pp$test$information, kits_information(item)(seq(-4, 4, by = 0.5)))
})
