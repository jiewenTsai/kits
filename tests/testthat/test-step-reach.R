# step_reach()：每個步驟上下的人數，CML 的存在條件、n_ig、kits_drift() 共用

resp_small <- data.frame(visit_id = rep(paste0("p", 1:6), each = 2),
                         item_id = rep(c("a", "b"), 6),
                         item_score = c(0, 0, 1, 0, 2, 1, 2, 1, 1, 0, 0, 1))

test_that("step_reach(): n_below and n_above per step, in the order of the steps", {
  steps <- data.frame(item_id = c("b", "a", "a", "c"), item_score = c(1, 1, 2, 1))
  out <- kits:::step_reach(resp_small, steps)
  expect_identical(out$item_id, steps$item_id)
  expect_identical(out$item_score, steps$item_score)
  brute <- function(i, k) {
    s <- resp_small$item_score[resp_small$item_id == i]
    c(sum(s < k), sum(s >= k))
  }
  expect_identical(out$n_below, c(brute("b", 1)[1], brute("a", 1)[1], brute("a", 2)[1], 0L))
  expect_identical(out$n_above, c(brute("b", 1)[2], brute("a", 1)[2], brute("a", 2)[2], 0L))
})

test_that("n_ig of the CML statistics is n_above of step_reach()", {
  st <- kits:::step_table(resp_small, c("a", "b"), c(2L, 1L))
  expect_identical(st$n, kits:::step_reach(resp_small, st)$n_above)
  expect_identical(st$n, c(4L, 2L, 3L))
})

test_that("CML existence: the message names the steps reached (missed) by no informative person", {
  X <- kits_sim(rnorm(300), list(b1 = -1, b2 = 0, b3 = c(0, 1), b4 = 1), seed = 3)
  X[, "b4"] <- as.integer(rowSums(X, na.rm = TRUE) == 5)
  expect_error(kits_item(X, scale = "x"),
               "Item step\\(s\\) with no finite CML estimate: b4:1 \\(reached by none\\)\\. ",
               class = "kits_error_items")
  # 固定的步驟不檢查
  fixed <- data.frame(item_id = "b4", item_score = 1, beta = 3)
  expect_s3_class(suppressMessages(kits_item(X, scale = "x", fixed_params = fixed)), "kits_item")
})

test_that("kits_drift(): n_below and n_above of every common step", {
  betas <- list(a1 = c(-1, 0), a2 = c(0, 1), a3 = 0.5, a4 = -0.5, a5 = 1.5)
  ref <- kits_item(kits_sim(seq(-3, 3, length.out = 300), betas, seed = 1), scale = "ref")
  X_new <- kits_sim(seq(-2, 4, length.out = 300), betas, prefix = "n", seed = 2)
  pts <- kits_drift(X_new, ref)$points
  below <- mapply(\(i, k) sum(X_new[, i] < k, na.rm = TRUE), pts$item_id, pts$item_score)
  above <- mapply(\(i, k) sum(X_new[, i] >= k, na.rm = TRUE), pts$item_id, pts$item_score)
  expect_identical(pts$n_below, unname(below))
  expect_identical(pts$n_above, unname(above))
})

test_that("shash: one set of bounds for the fit and the norm table", {
  expect_identical(kits:::SHASH_BOUNDS, c(log_sigma = 10, epsilon = 10, log_delta = 2))
  model <- list(age_mean = 0, age_sd = 1, mu_est = c(1, 2), sigma_est = c(0, 30),
                epsilon_est = c(0, -30), delta_est = c(0, 5))
  p <- kits:::shash_parameters(model, c(-1, 0, 1))
  expect_identical(names(p), c("mu", "sigma", "epsilon", "delta"))
  expect_identical(p$sigma, exp(c(-10, 0, 10)))
  expect_identical(p$epsilon, c(10, 0, -10))
  expect_identical(p$delta, exp(c(-2, 0, 2)))
})
