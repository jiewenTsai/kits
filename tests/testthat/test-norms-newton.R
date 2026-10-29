# 常模配適的牛頓步（shash_newton()）：nlminb 差一點收斂時才用；只動不在界上的係數，-log L 不變大

test_that("shash_newton()：走到最佳解、留在界內、在界上的係數不動", {
  target <- c(1, 2, 3)
  A <- diag(c(1e4, 1, 100))
  obj <- list(fn = \(p) drop(t(p - target) %*% A %*% (p - target)),
              gr = \(p) drop(2 * A %*% (p - target)),
              he = \(p) 2 * A)
  lower <- c(-Inf, -5, -5)
  upper <- c(Inf, 5, 2.5)
  start <- list(par = c(1.001, 1.9, 2.5), objective = obj$fn(c(1.001, 1.9, 2.5)))
  out <- kits:::shash_newton(start, obj, lower, upper)
  expect_equal(out$par[1:2], target[1:2])   # 自由的係數到最佳解
  expect_equal(out$par[3], 2.5)              # 在上界的係數不動（最佳解 3 在界外）
  expect_lte(out$objective, start$objective)
  # Hessian 不正定：不走
  bad <- modifyList(obj, list(he = \(p) -2 * A))
  expect_identical(kits:::shash_newton(start, bad, lower, upper)$par, start$par)
})
