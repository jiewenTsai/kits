# 診斷、summary()、kits_check()、來源紀錄（kits_info）、有類別的錯誤和警告
set.seed(3)

betas <- list(a1 = c(-1, 0), a2 = c(0, 1), a3 = 0.5, a4 = -0.5, a5 = c(-2, -1, 1))
X <- kits_sim(rnorm(400, 0, 1.5), betas, seed = 1)
item <- kits_item(X)
s <- kits_subj(X, item)
# kitsPreview 載入時會蓋掉一些 print 方法，所以印出來的檢查直接叫 kits 的方法
pr <- function(x, ...) {
  fun <- get(paste0("print.", class(x)[1]), envir = asNamespace("kits"))
  fun(x, ...)
}

test_that("診斷表的欄位；kits_item 的檢查都在", {
  dg <- item$diagnostics
  expect_s3_class(dg, "kits_diagnostics")
  expect_equal(names(dg), c("check", "value", "status", "message"))
  expect_true(all(dg$status %in% c("ok", "note", "warn")))
  expect_equal(dg$check, c("optimizer convergence", "max |gradient|", "Hessian positive definite", "extreme estimates",
                           "smallest category count",
                           "fixed items left out", "persons", "item patterns",
                           "smallest main item pattern", "item fit"))
  expect_true(all(dg$status == "ok"))
  expect_lt(dg$value[dg$check == "max |gradient|"], 1e-4)
  expect_output(pr(item), "Checks: OK")
  # 很少人的分數：note，說出題號和分數
  Y <- X
  Y[Y[, "a5"] == 3, "a5"][-(1:2)] <- 2L
  dy <- kits_item(Y)$diagnostics
  expect_equal(dy$status[dy$check == "smallest category count"], "note")
  expect_match(dy$message[dy$check == "smallest category count"], "a5:3 \\(2\\)")
  # 固定的題沒放進校準：note
  Z <- X
  Z[, "a1"] <- pmax(Z[, "a1"], 1L)
  dz <- suppressMessages(kits_item(Z, fixed_params = item))$diagnostics
  expect_equal(dz$status[dz$check == "fixed items left out"], "note")
})

test_that("kits_subj 的檢查：0 分、滿分、不適配的比例", {
  dg <- attr(s, "diagnostics")
  expect_equal(dg$check[1:3], c("share of zero scores", "share of perfect scores", "share of misfit"))
  expect_equal(dg$value[1], 100 * mean(s$flag == "zero"))
  hard <- kits_sim(rnorm(300, -4), betas, seed = 2)
  dh <- attr(kits_subj(hard, item), "diagnostics")
  expect_equal(dh$status[dh$check == "share of zero scores"], "note")
  # 每個波次都一樣的孩子：信度低於 0（警告）
  flat <- kits_sim(rep(-3, 400), betas, seed = 3)   # 0.25.0 的 SE 比較小：在 0 時信度剛好在 0 附近
  df <- attr(kits_subj(flat, item, persons = data.frame(visit_id = rownames(flat), wave_m = 12)),
             "diagnostics")
  expect_equal(df$status[df$check == "reliability, wave 12"], "warn")
  expect_match(df$message[df$check == "reliability, wave 12"], "cannot separate children")
})

test_that("連結的檢查：留下的點、剔除的比例、相關、連結誤差、各領域的截距", {
  from <- seq(-2, 3, length.out = 12)
  to <- 1 + 2 * from
  to[5] <- to[5] + 4
  link <- link_of_points(from, to, id = paste0("p", 1:12))
  dg <- link$diagnostics
  expect_equal(dg$check, c("points kept", "share dropped", "correlation of kept points"))   # n_boot = 0: no linking error row
  expect_equal(dg$value[1:2], c(11, 100 / 12))
  expect_equal(dg$status, c("ok", "ok", "ok"))
  expect_output(pr(link), "Checks: OK")
  expect_false(any(grepl("linking error", capture.output(pr(link)))))
  boot <- link_of_points(from, to, id = paste0("p", 1:12), n_boot = 50)
  expect_equal(boot$diagnostics$check[4], "linking error")
  expect_output(pr(boot), "linking error")
  expect_warning(one <- link_of_points(from, to, id = paste0("p", 1:12), item = rep("i", 12),
                                       n_boot = 50), "one item")
  expect_equal(one$diagnostics$status[one$diagnostics$check == "linking error"], "note")
  expect_no_warning(link_of_points(from, to, id = paste0("p", 1:12), item = rep("i", 12)))
  few <- link_of_points(c(-1, 0, 1, 2), c(0, 1, 2, 3.5), letters[1:4], n_boot = 50)
  expect_equal(few$diagnostics$status[1], "warn")
  expect_equal(few$diagnostics$status[4], "ok")
  noisy <- link_of_points(1:8, c(1, 5, 2, 7, 3, 4, 8, 6), letters[1:8])
  expect_equal(noisy$diagnostics$status[noisy$diagnostics$check == "correlation of kept points"], "warn")
  expect_output(pr(noisy), "Checks: 1 warning")
  steps <- data.frame(item_id = paste0("a", 1:8), item_score = 1, thr = c(-1, 0, 1, 2, -1.5, 0.5, 1, 2.5),
                      domain = rep(c("lan", "cog"), each = 4))
  anchors <- data.frame(item_id = paste0("a", 1:8), tau = 64 + 7 * steps$thr + ifelse(1:8 > 4, 3, 0) + c(0.1, -0.1, 0.2, 0, 0.1, 0, -0.2, 0.1))
  dl <- kits_linking(steps, anchors, n_boot = 200)
  row <- dl$diagnostics[dl$diagnostics$check == "max |domain intercept - intercept|", ]
  expect_equal(row$value, max(abs(dl$domains$intercept - dl$intercept)))
  expect_equal(row$status, "note")
})

test_that("常模、篩波次、交互作用模型、偏移的檢查", {
  set.seed(1)
  age <- runif(1500, 3, 36)
  D <- 20 + 15 * log(age) + rnorm(1500, 0, 3)
  norm <- suppressWarnings(kits_norms(D, age, ages = 6:30, degrees = c(3, 2, 1, 1)))
  dg <- attr(norm, "diagnostics")
  expect_equal(dg$check, c("optimizer convergence", "max |gradient|", "coefficients at their limits", "max delta",
                           "same-age z: max |mean| by age decile", "same-age z: SD furthest from 1 by age decile"))
  expect_equal(dg$status, rep("ok", 6))
  expect_lte(dg$value[2], 1e-4)
  expect_output(pr(norm), "Checks:")
  # 中位數不對的常模：同齡 z 的平均偏掉
  bent <- suppressWarnings(kits_norms(D + ifelse(age > 20, 6, 0), age, ages = 6:30, degrees = c(1, 1, 1, 1)))
  expect_true("warn" %in% attr(bent, "diagnostics")$status[5:6])
  D2 <- c(10, 10, 10, 12, 15, 20, 22, 25, 30, 31)
  scr <- kits_norm_screen(D2, rep(c(6, 12), each = 5))
  expect_equal(scr$diagnostics$status, "note")
  expect_match(scr$diagnostics$message, "6, 12")
  new <- kits_sim(seq(-2, 4, length.out = 400), modifyList(betas, list(a3 = 2.5)), prefix = "n", seed = 2)
  drift <- kits_drift(new, item)
  expect_s3_class(drift, "kits_drift")
  expect_equal(drift$diagnostics$status[1], "note")
  expect_match(drift$diagnostics$message[1], "a3")
  expect_output(pr(drift), "drifted: a3")
  # 印出來的建議不假設資料叫 X_new：矩陣是 <new data>，data.frame 是它的名字
  expect_output(pr(drift), "kits_item(<new data>, fixed_params = x$fix)", fixed = TRUE)
  d_new <- data.frame(visit_id = rownames(new), new)
  drift_d <- kits_drift(d_new, item)
  expect_equal(kits_info(drift_d)$input$name, "d_new")
  expect_output(pr(drift_d), "kits_item(d_new, fixed_params = x$fix)", fixed = TRUE)
  expect_equal(drift_d$items, drift$items)
})

test_that("kits_itemfit 用指紋看參數是不是這份資料校準的", {
  same <- attr(kits_itemfit(X, item), "diagnostics")
  expect_equal(same$status[1], "ok")
  other <- attr(kits_itemfit(X[1:300, ], item), "diagnostics")
  expect_equal(other$status[1], "note")
  expect_match(other$message[1], "other data")
  none <- attr(kits_itemfit(X, as.data.frame(coef(item))[c("item_id", "item_score", "beta")]), "diagnostics")
  expect_match(none$message[1], "no record")
})

test_that("kits_check：一個或多個結果；只留 note 和 warn，warn 在前；all = TRUE 全部", {
  noisy <- link_of_points(1:8, c(1, 5, 2, 7, 3, 4, 8, 6), letters[1:8])
  out <- kits_check(item, s, noisy)
  expect_s3_class(out, "kits_table")
  expect_equal(names(out), c("object", "check", "value", "status", "message"))
  expect_true(all(out$status != "ok"))
  expect_equal(out$status[1], "warn")
  expect_equal(out$object[1], "noisy")
  expect_false(is.unsorted(match(out$status, c("warn", "note"))))
  all_rows <- kits_check(calibration = item, all = TRUE)
  expect_equal(unique(all_rows$object), "calibration")
  expect_equal(nrow(all_rows), nrow(item$diagnostics))
  expect_output(print(kits_check(item)), "Checks: OK \\(no notes or warnings\\)")
  expect_output(print(out), "correlation of kept points")
  expect_error(kits_check(), class = "kits_error_input")
  # an object that is not a kits result: a row (warn), not an error, so the other objects are still checked
  bad <- kits_check(item, "hello", X)
  expect_equal(bad$object[bad$status == "warn"], c("object 2", "X"))   # a value is labelled by its position, never by its data
  expect_match(bad$message[bad$object == "X"], "not a result of a kits function \\(a matrix\\)")
  expect_false(any(grepl("map2", bad$message)))
})

test_that("summary()：簡短的印出加上整張診斷表", {
  out <- capture.output(print(summary(s)))
  expect_true(any(grepl("<kits_subj> 400 person-waves", out)))
  expect_true(any(grepl("Diagnostics:", out)))
  expect_true(any(grepl("share of misfit", out)))
  expect_true(any(grepl("Diagnostics:", capture.output(print(summary(item))))))
  short <- link_of_points(1:5, c(1, 2, 3, 4, 6), letters[1:5])
  expect_true(any(grepl("points kept", capture.output(print(summary(short))))))
  expect_true(any(grepl("<kits_itemfit>", capture.output(print(summary(kits_itemfit(X, item)))))))
})

test_that("同樣的輸入給完全相同的物件；kits_info 記下呼叫、版本、設定、資料指紋", {
  expect_identical(kits_item(X), kits_item(X))
  expect_identical(kits_subj(X, item), kits_subj(X, item))
  info <- kits_info(item)
  expect_s3_class(info, "kits_info")
  expect_equal(info$call, "kits_item(X)")
  expect_equal(info$version, as.character(utils::packageVersion("kits")))
  expect_equal(info$data$n_rows, 400)
  expect_equal(info$data$items, names(betas))
  expect_match(info$data$scores_md5, "^[0-9a-f]{32}$")
  # 矩陣和長格式是同一份資料：指紋相同
  long <- data.frame(visit_id = rep(rownames(X), ncol(X)),
                     item_id = rep(colnames(X), each = nrow(X)), item_score = as.vector(X))
  expect_equal(kits_info(kits_item(long[sample(nrow(long)), ]))$data, info$data)
  Y <- X
  Y[1, 1] <- 1L - Y[1, 1]
  expect_false(identical(kits_info(kits_item(Y))$data$scores_md5, info$data$scores_md5))
  expect_output(print(info), "Call: kits_item\\(X\\)")
  expect_output(print(kits_info(kits_linking(item, data.frame(item_id = c("a1", "a2", "a3", "a4"), tau = c(40, 50, 45, 55)), n_boot = 0))),
                "kits_linking\\(item, ")
  expect_equal(kits_info(kits_linking(item, data.frame(item_id = c("a1", "a2", "a3", "a4"), tau = c(40, 50, 45, 55)), n_boot = 200))$call,
               "kits_linking(item, data.frame(item_id = c(<data>, <data>, <data>, <data>), tau = c(<data>, <data>, <data>, <data>)), n_boot = 200)")
  expect_error(kits_info(1), class = "kits_error_input")
})

test_that("錯誤、警告、訊息都有類別", {
  e <- tryCatch(kits_item(X, scale = "a", method = "MML"), error = identity)
  expect_s3_class(e, c("kits_error_input", "kits_error", "error", "condition"))
  expect_null(conditionCall(e))
  gap <- X
  gap[gap[, "a3"] == 1, "a3"] <- 2L
  gap[gap[, "a1"] == 1, "a1"] <- 2L
  expect_error(kits_item(gap), class = "kits_error_data")
  expect_error(kits_subj(cbind(X, zz = 0L), item), class = "kits_error_items")
  z <- X
  z[] <- 0L
  expect_error(suppressMessages(kits_item(z, fixed_params = item)), class = "kits_error_fit")
  expect_warning(link_of_points(1:5, -(1:5), letters[1:5]), class = "kits_warning_scale")
  set.seed(1)
  age <- runif(1500, 3, 36)
  D <- 20 + 15 * log(age) + rnorm(1500, 0, 3)
  expect_warning(kits_norms(D, age, ages = 6:24, degrees = c(3, 2, 1, 1), delta_max = 0), class = "kits_warning_norms")
  Y <- X
  Y[, "a1"] <- pmax(Y[, "a1"], 1L)
  expect_message(kits_item(Y, fixed_params = item), class = "kits_message_fit")
  expect_message(kits_subj(rbind(X, zz = NA), item), class = "kits_message_data")
  caught <- tryCatch(kits_linking(item, item), kits_error_scale = function(e) "caught")
  expect_equal(caught, "caught")
})
