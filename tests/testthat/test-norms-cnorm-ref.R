# 常模和 cNORM 比：cNORM 3.7.0 算的值存在 inst/validation/reference/norms_cnorm.rds
# （reference/make_norms_cnorm.R 產生；資料是模擬的）。cNORM 沒裝也能比。
test_that("kits_norms、kits_norm_z 和存下的 cNORM 值相同（負對數概似、百分位、同齡 z 在 1e-4 內）", {
  dir <- system.file("validation", package = "kits")
  skip_if(!nzchar(dir) || !file.exists(file.path(dir, "reference", "norms_cnorm.rds")))
  source(file.path(dir, "norms_sim.R"), local = TRUE)
  ns <- norm_sim()
  ref <- readRDS(file.path(dir, "reference", "norms_cnorm.rds"))
  nr <- suppressWarnings(kits_norms(ns$D, ns$age, ages = ns$ages, degrees = ns$degrees, pad = 0,
                                    delta_starts = ns$delta_starts))
  expect_lt(abs(nr$neg_loglik[1] - ref$neg_loglik), 1e-4)
  expect_lt(max(abs(as.matrix(as.data.frame(nr)[colnames(ref$pct)]) - ref$pct)), 1e-4)
  inside <- abs(ref$z) < 3   # cNORM 預設把常模分數截在 z = ±3
  expect_lt(max(abs(kits_norm_z(ns$D, ns$age_row, nr)[inside] - ref$z[inside])), 1e-4)
})
