# kits 0.23.0：kits_item_restscore()，Kreiner (2011) 的題目–其餘分數 gamma。
# 完整資料時和 iarm::item_restscore() 相同；不完整設計時每個題目組合各算一次，小的組合跳過。

# kits 的方法（NAMESPACE 還沒登記時也能測）
rs_method <- function(generic, x, ...) {
  get(paste0(generic, ".kits_item_restscore"), envir = asNamespace("kits"))(x, ...)
}

# 參考值：iarm::item_restscore()（iarm 0.4.3，eRm 1.0.10）和 vcdExtra::GKgamma() 的數字，
# 由 fixtures/fit_ref/make_fit_ref.R 算好存下來；kits 不依賴這些套件
fit_ref <- function(name) readRDS(test_path("fixtures", "fit_ref", paste0(name, ".rds")))

betas <- list(a1 = c(-1, 0.2), a2 = c(-0.5, 0, 1), a3 = 0.5, a4 = c(-0.3, 0.6),
              a5 = c(0, 0.4, 1.1), a6 = -0.7)
set.seed(3)
X <- kits_sim(rnorm(800, 0, 1.3), betas, seed = 11)
item <- kits_item(X)
rs <- kits_item_restscore(X, item)

test_that("the table has one row per item on complete data", {
  expect_s3_class(rs, c("kits_item_restscore", "kits_table", "data.frame"))
  expect_identical(names(rs), c("item_id", "pattern_id", "n", "gamma_obs", "gamma_exp", "se", "z",
                                "p", "p_adj"))
  expect_identical(rs$item_id, colnames(X))
  expect_identical(unique(rs$pattern_id), "bk1")
  expect_true(all(rs$n == 800))
  expect_equal(rs$z, (rs$gamma_obs - rs$gamma_exp) / rs$se)
  expect_equal(rs$p_adj, stats::p.adjust(rs$p, "BH"))
  expect_true(all(rs$se > 0 & rs$p >= 0 & rs$p <= 1))
})

test_that("complete PCM data: the same numbers as iarm::item_restscore(eRm::PCM(X))", {
  ref <- fit_ref("restscore_pcm")
  expect_identical(unname(ref$data), unname(X))   # the same simulated data
  # observed gamma and its SE use only the data: identical
  expect_equal(rs$gamma_obs, unname(ref$ref$observed), tolerance = 1e-12)
  expect_equal(rs$se, unname(ref$ref$se), tolerance = 1e-12)
  # expected gamma and p use the item parameters: kits' CML and eRm's agree to eRm's convergence
  expect_lt(max(abs(rs$gamma_exp - ref$ref$expected)), 5e-6)
  expect_lt(max(abs(rs$p - ref$ref$p)), 1e-4)
  expect_lt(max(abs(rs$p_adj - ref$ref$p_adj)), 1e-4)
  # with eRm's item parameters, everything agrees to rounding
  same <- kits_item_restscore(X, ref$parms)
  expect_equal(same$gamma_exp, unname(ref$ref$expected), tolerance = 1e-10)
  expect_equal(same$p, unname(ref$ref$p), tolerance = 1e-10)
  expect_equal(same$p_adj, unname(ref$ref$p_adj), tolerance = 1e-10)
})

test_that("dichotomous items: the same numbers as iarm on eRm::RM()", {
  b2 <- list(d1 = -1.2, d2 = -0.5, d3 = 0, d4 = 0.3, d5 = 0.8, d6 = 1.4)
  set.seed(4)
  X2 <- kits_sim(rnorm(600), b2, seed = 5)
  item2 <- kits_item(X2)
  rs2 <- kits_item_restscore(X2, item2)
  expect_identical(nrow(rs2), 6L)
  expect_true(all(abs(rs2$gamma_obs) <= 1 & abs(rs2$gamma_exp) <= 1))
  ref <- fit_ref("restscore_rm")
  expect_identical(unname(ref$data), unname(X2))
  expect_equal(rs2$gamma_obs, unname(ref$ref$observed), tolerance = 1e-12)
  expect_equal(rs2$se, unname(ref$ref$se), tolerance = 1e-12)
  expect_lt(max(abs(rs2$gamma_exp - ref$ref$expected)), 5e-6)
  expect_lt(max(abs(rs2$p - ref$ref$p)), 1e-4)
  same <- kits_item_restscore(X2, ref$parms)
  expect_equal(same$gamma_exp, unname(ref$ref$expected), tolerance = 1e-10)
})

test_that("the expected table follows the PCM conditional distribution (brute force)", {
  b <- list(c(-0.4, 0.3), 0.2, c(-1, 0.5, 0.9))
  grid <- expand.grid(map(b, \(v) 0:length(v)))
  cum_beta <- sapply(seq_along(b), \(i) c(0, cumsum(b[[i]]))[grid[[i]] + 1])
  weight <- exp(-rowSums(cum_beta))
  total <- rowSums(grid)
  n_total <- c(3, 10, 20, 25, 18, 9, 4)   # persons with total 0..6
  ex <- expected_rest_table(algo_esf(stats::setNames(b, c("i1", "i2", "i3")))$prob[[1]], n_total)
  brute <- matrix(0, 5, 3)
  for (k in seq_len(nrow(grid))) {
    t <- total[k]
    share <- weight[k] / sum(weight[total == t])
    brute[t - grid[k, 1] + 1, grid[k, 1] + 1] <- brute[t - grid[k, 1] + 1, grid[k, 1] + 1] +
      n_total[t + 1] * share
  }
  expect_equal(ex, brute, tolerance = 1e-12)
  expect_equal(sum(ex), sum(n_total))
})

test_that("gamma and its SE are those of vcdExtra::GKgamma()", {
  ref <- fit_ref("gkgamma")
  for (k in seq_along(ref$tables)) {
    got <- gk_gamma(ref$tables[[k]])
    expect_equal(got$gamma, ref$gamma[k], tolerance = 1e-12)
    expect_equal(got$se, ref$sigma[k], tolerance = 1e-12)
  }
  expect_true(is.na(gk_gamma(matrix(c(5, 0, 3, 0), 2))$gamma))   # one row only: no pairs
})

test_that("incomplete designs: one row per item and item pattern; small patterns skipped", {
  b3 <- list(a1 = c(-1, 0), a2 = c(-0.5, 0.5), a3 = 0.3, a4 = c(-0.2, 0.8), a5 = -0.4, a6 = 0.6)
  set.seed(7)
  X3 <- kits_sim(rnorm(1250, 0, 1.3), b3, seed = 8)
  X3[1:600, "a6"] <- NA         # pattern a1..a5
  X3[601:1180, "a1"] <- NA      # pattern a2..a6
  X3[1181:1240, "a3"] <- NA     # 60 persons: skipped
  item3 <- kits_item(X3)
  rs3 <- kits_item_restscore(X3, item3)
  expect_identical(nrow(rs3), 10L)
  expect_setequal(rs3$n, c(600L, 580L))
  expect_identical(as.vector(table(rs3$item_id)), c(1L, 2L, 2L, 2L, 2L, 1L))
  skipped <- attr(rs3, "skipped")
  expect_identical(sort(skipped$n), c(10L, 60L))   # 60 without a3, 10 with all items
  d <- attr(rs3, "diagnostics")
  expect_identical(d$status[d$check == "item patterns skipped"], "note")
  # the pattern ids are those of kits_subj()
  subj <- kits_subj(X3, item3)
  ids <- unique(subj$pattern_id[subj$visit_id %in% rownames(X3)[1:600]])
  expect_true(all(rs3$pattern_id[rs3$n == 600] == ids))
  # with min_n = 0 every pattern with 2 or more items is tested
  expect_identical(nrow(kits_item_restscore(X3, item3, min_n = 0)), 10L + 5L + 6L)
  # each pattern is the complete-data test on its own items, with the same item parameters
  sub <- X3[1:600, paste0("a", 1:5)]
  ref <- kits_item_restscore(sub, item3)
  got <- rs3[rs3$n == 600, ]
  expect_equal(got$gamma_obs, ref$gamma_obs, tolerance = 1e-12)
  expect_equal(got$gamma_exp, ref$gamma_exp, tolerance = 1e-12)
  expect_equal(got$se, ref$se, tolerance = 1e-12)
  expect_equal(got$z, ref$z, tolerance = 1e-12)
})

test_that("a misfitting item is found", {
  set.seed(1)
  Xn <- kits_sim(seq(-3, 3, length.out = 500), betas, seed = 1)
  Xn[, "a3"] <- sample(0:1, 500, replace = TRUE)
  rsn <- kits_item_restscore(Xn, kits_item(Xn))
  expect_lt(rsn$p_adj[rsn$item_id == "a3"], 1e-6)
  expect_lt(rsn$gamma_obs[rsn$item_id == "a3"], rsn$gamma_exp[rsn$item_id == "a3"])
})

test_that("bad inputs stop with kits_error_input", {
  expect_error(kits_item_restscore(X, "not an item"), class = "kits_error_input")
  expect_error(kits_item_restscore(X, item, min_n = -1), class = "kits_error_input")
  expect_error(kits_item_restscore(X, item, min_n = "a"), class = "kits_error_input")
  expect_error(kits_item_restscore(X, item, min_n = c(10, 20)), class = "kits_error_input")
  expect_error(kits_item_restscore(X, item, alpha = 0.01), class = "kits_error_input")
  expect_error(kits_item_restscore(X[, 1], item), class = "kits_error_input")
  expect_error(kits_item_restscore(X, item, min_n = 1000), class = "kits_error_data")
  extra <- cbind(X, z1 = 0L)
  expect_error(kits_item_restscore(extra, item), class = "kits_error_items")
})

test_that("print shows the counts, the significant and the largest rows, and Next", {
  out <- utils::capture.output(rs_method("print", rs))
  expect_match(out[1], "6 rows")
  expect_true(any(grepl("p_adj < 0.05", out)))
  expect_true(any(grepl("Largest", out)))
  expect_true(any(grepl("large N", out)))
  expect_true(any(grepl("^Next:", out)))
})

test_that("plot() draws observed against expected gamma with the identity line", {
  p <- rs_method("plot", rs)
  expect_s3_class(p, "plotly")
  traces <- p$x$data
  expect_identical(traces[[1]]$name, "observed = expected")
  expect_equal(sum(map_int(traces[-1], \(tr) length(tr$x))), nrow(rs))
  expect_error(rs_method("plot", rs, items = "zz"), class = "kits_error_input")
})
