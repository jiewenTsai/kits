# kits_dif() 的 LR 檢定：Andersen (1973)。
# eRm::LRtest() 的參考值存在 fixtures/fit_ref/lrtest_ref.rds（eRm 1.0.10 算的；eRm 不是 kits 的相依套件，
# 這裡不呼叫 eRm）。

ref <- readRDS(test_path("fixtures", "fit_ref", "lrtest_ref.rds"))

# 結果的各部分（kits_dif() 的屬性）
lr <- function(x) attr(x, "lrtest")
part <- function(x, what) as.data.frame(attr(x, what))
expect_lr <- function(x, r, tol = 1e-6) {
  expect_equal(lr(x)$LR, r$LR, tolerance = tol)
  expect_identical(lr(x)$df, as.integer(r$df))
  expect_equal(lr(x)$p, r$p, tolerance = tol)
}

# ---- 和 eRm 相同（完整資料，沒有觀察不到的步驟） ----------------------------------------------------

test_that("依原始分的中位數分組：LR、df、p 和 eRm::LRtest(splitcr = \"median\") 相同", {
  x <- kits_dif(ref$X, "raw")
  expect_s3_class(x, "kits_dif")
  expect_lr(x, ref$median)
  expect_identical(lr(x)$groups$group, c("low", "high"))
  expect_identical(lr(x)$groups$n, unname(ref$median$n))
  expect_identical(nrow(part(x, "unobserved")), 0L)
  # 有 kits_item 時用它的概似，和自己校準的結果相同
  item <- kits_item(ref$X, scale = "lr")
  expect_equal(lr(kits_dif(ref$X, "raw", item))$LR, lr(x)$LR, tolerance = 1e-8)
  expect_equal(attr(attr(x, "full"), "logLik"), as.numeric(logLik(item)), tolerance = 1e-8)
})

test_that("同分歸低分組（和 eRm 相同）：raw <= 中位數的是 low", {
  rv <- rowSums(ref$X)
  x <- kits_dif(ref$X, "raw")
  expect_true(any(rv == stats::median(rv)))   # 真的有同分
  expect_identical(lr(x)$groups$n, c(sum(rv <= stats::median(rv)), sum(rv > stats::median(rv))))
})

test_that("使用者給的分組向量：LR、df、p 和 eRm::LRtest(splitcr = 向量) 相同", {
  x <- kits_dif(ref$X, ref$split3)
  expect_lr(x, ref$user)
  expect_identical(lr(x)$groups$group, c("1", "2", "3"))
  expect_identical(lr(x)$groups$n, unname(ref$user$n))
  # 因子：組照因子的水準排；照 visit_id 命名的向量也可以
  f <- factor(ref$split3, levels = c(3, 1, 2))
  expect_identical(lr(kits_dif(ref$X, f))$groups$group, c("3", "1", "2"))
  named <- stats::setNames(rev(ref$split3), rev(rownames(ref$X)))
  expect_equal(lr(kits_dif(ref$X, named))$LR, lr(x)$LR, tolerance = 1e-10)
})

test_that("每組的步驟參數：和全部一起的 beta 同一個原點；差和 Wald z", {
  x <- kits_dif(ref$X, "raw")
  s <- as.data.frame(part(x, "steps"))
  # 完整資料：每組的 beta 和是 0，和全部一起的相同，所以不用平移
  expect_equal(as.vector(tapply(s$beta, s$group, sum)), c(0, 0), tolerance = 1e-8)
  d <- as.data.frame(part(x, "diff"))
  low <- s[s$group == "low", ]
  high <- s[s$group == "high", ]
  expect_equal(d$diff, high$beta - low$beta)
  expect_equal(d$z, (high$beta - low$beta) / sqrt(low$se^2 + high$se^2))
  expect_identical(unique(paste(d$group_1, d$group_2)), "low high")
})

# ---- 觀察不到的步驟 -------------------------------------------------------------------------------

test_that("一組觀察不到的步驟：這組只估觀察得到的步驟（和 eRm 的 PCM() 各組的概似相同）", {
  # 每題的 DIF：a5 兩組的分數不同，不比（警告）
  quiet_dif <- function(e) suppressWarnings(e, classes = "kits_warning_data")
  expect_message(x <- quiet_dif(kits_dif(ref$Xu, "raw")), class = "kits_message_data")
  expect_lr(x, ref$steps)
  expect_equal(lr(x)$groups$logLik, unname(ref$steps$loglik[c("low", "high")]), tolerance = 1e-6)
  expect_equal(attr(attr(x, "full"), "logLik"), unname(ref$steps$loglik["all"]), tolerance = 1e-6)
  u <- as.data.frame(part(x, "unobserved"))
  expect_identical(u$group, "low")
  expect_identical(u$item_id, "a5")
  expect_identical(u$item_score, 2L)
  expect_false(u$left_out)
  st <- part(x, "steps")
  expect_false(any(st$group == "low" & st$item_id == "a5" & st$item_score == 2))
  expect_identical(lr(x)$groups$n_par, c(7L, 8L))
  dg <- attr(x, "diagnostics")
  expect_true("steps not observed in a group" %in% dg$check[dg$status == "note"])
})

test_that("unobserved = \"drop\"：每一組都不用這題，和 eRm::LRtest() 相同", {
  drop_dif <- function() kits_dif(ref$Xu, "raw", unobserved = "drop")
  expect_warning(expect_warning(x <- drop_dif(), "joint calibration"), "different scores")
  expect_lr(x, ref$drop)
  expect_identical(lr(x)$groups$n, unname(ref$drop$n))
  expect_true(all(part(x, "unobserved")$left_out))
  expect_identical(unique(part(x, "unobserved")$item_id), "a5")
  expect_false("a5" %in% attr(x, "full")$item_id)
})

test_that("最低分觀察不到（沒有 0 分）：分數往下平移，估其他的步驟", {
  X <- kits_sim(seq(-2, 2, length.out = 600), list(a1 = c(-1, 0.5), a2 = c(0, 1), a3 = 0.3, a4 = -0.4),
                seed = 4)
  g <- rep(c("A", "B"), each = 300)
  X[g == "B" & X[, "a1"] == 0, "a1"] <- 1   # B 組沒有人 a1 是 0 分
  expect_message(x <- kits_dif(X, g), class = "kits_message_data")
  u <- as.data.frame(part(x, "unobserved"))
  expect_identical(paste(u$group, u$item_id, u$item_score), "B a1 1")
  expect_match(u$reason, "reached by every")
  # B 組的概似 = a1 減 1 分之後的校準
  Xb <- X[g == "B", ]
  Xb[, "a1"] <- Xb[, "a1"] - 1
  expect_equal(lr(x)$groups$logLik[2], -cml_fit(resp_long(Xb))$est$nll, tolerance = 1e-8)
  expect_true(any(part(x, "steps")$group == "B" & part(x, "steps")$item_id == "a1" & part(x, "steps")$item_score == 2))
})

test_that("一組裡分數有缺（0、1、3）：這組不用這題，全部一起的校準也不用這些作答（警告）", {
  X <- kits_sim(seq(-2, 2, length.out = 600), list(a1 = c(-1, 0, 1), a2 = c(0, 1), a3 = 0.3, a4 = -0.4),
                seed = 5)
  g <- rep(c("A", "B"), each = 300)
  X[g == "A" & X[, "a1"] == 2, "a1"] <- 1   # A 組 a1 只有 0、1、3 分
  expect_true(all(c(0, 1, 3) %in% X[g == "A", "a1"]))
  expect_warning(x <- suppressMessages(kits_dif(X, g)), class = "kits_warning_data")
  u <- as.data.frame(part(x, "unobserved"))
  expect_true(all(u$left_out[u$item_id == "a1"]))
  expect_match(u$reason[u$item_id == "a1"][1], "gap")
  # 全部一起的校準：A 組的 a1 不用
  joint <- resp_long(X)
  joint <- joint[!(joint$item_id == "a1" & joint$visit_id %in% rownames(X)[g == "A"]), ]
  expect_equal(attr(attr(x, "full"), "logLik"), -cml_fit(joint)$est$nll, tolerance = 1e-8)
  expect_true(any(attr(x, "diagnostics")$check == "items left out" & attr(x, "diagnostics")$status == "note"))
})

# ---- 不完整的設計 ---------------------------------------------------------------------------------

betas7 <- list(a1 = c(-1, 0), a2 = -0.5, a3 = c(0, 0.5), a4 = 0.2, a5 = c(-0.3, 0.6), a6 = 0.8,
               a7 = -0.2)
X7 <- kits_sim(rnorm(900), betas7, seed = 6)
book <- rep(c("A", "B", "C"), 300)
X7[book == "A", c("a5", "a6", "a7")] <- NA
X7[book == "B", c("a1", "a2", "a7")] <- NA
X7[book == "C", c("a1", "a2", "a3")] <- NA

test_that("不完整的設計（好幾個作答組合）：在每個作答組合裡依中位數分組", {
  x <- suppressMessages(kits_dif(X7, "raw"))
  rv <- rowSums(X7, na.rm = TRUE)
  low <- unlist(tapply(seq_along(rv), book, \(i) rv[i] <= stats::median(rv[i])))
  expect_identical(lr(x)$groups$n, c(sum(low), sum(!low)))
  expect_equal(lr(x)$LR, 2 * (sum(lr(x)$groups$logLik) - attr(attr(x, "full"), "logLik")))
  expect_identical(lr(x)$df, sum(lr(x)$groups$n_par) - attr(attr(x, "full"), "n_par"))
  expect_gt(lr(x)$p, 0.001)   # 資料照 PCM 模擬
  expect_identical(attr(attr(x, "full"), "n_par"), length(unlist(betas7)) - 1L)
})

test_that("一組答的題分成沒有共同題的兩片：各自校準，自由參數是每片的步驟數減 1", {
  g2 <- ifelse(book == "B", "mid", "ends")   # ends 組：A 和 C
  X2 <- X7
  X2[book == "C", "a4"] <- NA   # 拿掉 a4：A 和 C 沒有共同題，B 把它們連起來
  X2[book == "A", "a4"] <- NA
  x <- suppressMessages(kits_dif(X2, g2))
  ends <- lr(x)$groups[lr(x)$groups$group == "ends", ]
  steps_a <- 2 + 1 + 2   # a1, a2, a3
  steps_c <- 2 + 1 + 1   # a5, a6, a7
  expect_identical(ends$n_steps, as.integer(steps_a + steps_c))
  expect_identical(ends$n_par, as.integer(steps_a - 1 + steps_c - 1))
  r <- suppressMessages(resp_long(X2[g2 == "ends", ]))   # a4 沒有作答，拿掉（訊息）
  ll <- -cml_fit(r[r$item_id %in% c("a1", "a2", "a3"), ])$est$nll -
    cml_fit(r[r$item_id %in% c("a5", "a6", "a7"), ])$est$nll
  expect_equal(ends$logLik, ll, tolerance = 1e-8)
})

# ---- 依波次、依孩子的欄 ---------------------------------------------------------------------------

lan <- kitm03$responses[kitm03$responses$domain == "lan", ]
resp <- kits_resp(lan)

test_that("split = \"wave\"：kitm03 lan，每個波次一組，各組只用觀察得到的步驟", {
  item <- kits_item(resp, scale = "lan M03")
  x <- suppressWarnings(suppressMessages(kits_dif(resp, "wave", item)))
  expect_identical(lr(x)$groups$group, as.character(sort(unique(resp$wave_m))))
  expect_identical(sum(lr(x)$groups$n), nrow(resp))
  expect_gt(lr(x)$df, 0)
  expect_equal(lr(x)$LR, 2 * (sum(lr(x)$groups$logLik) - attr(attr(x, "full"), "logLik")))
  expect_identical(lr(x)$df, sum(lr(x)$groups$n_par) - attr(attr(x, "full"), "n_par"))
  expect_gt(nrow(part(x, "unobserved")), 0)
  expect_identical(lr(x)$split, "wave_m")
  # 波次 3 沒有人 lana04 拿到最高分
  expect_true(any(part(x, "unobserved")$group == "3" & part(x, "unobserved")$item_id == "lana04"))
  # 訊息（沒有估的步驟）和警告（分數有缺，lanc09 在 18 個月）
  expect_message(suppressWarnings(kits_dif(resp, "wave", item)), class = "kits_message_data")
  expect_warning(suppressMessages(kits_dif(resp, "wave", item)), class = "kits_warning_data")
  # 矩陣加上 persons 給一樣的結果
  X <- wide_matrix(resp, "resp")
  persons <- data.frame(visit_id = resp$visit_id, wave_m = resp$wave_m)
  y <- suppressWarnings(suppressMessages(kits_dif(X, "wave", persons = persons)))
  expect_equal(lr(y)$LR, lr(x)$LR, tolerance = 1e-8)
  expect_identical(lr(y)$df, lr(x)$df)
})

test_that("孩子的欄（sex）、沒有組（NA）的人次不用（訊息）", {
  x <- suppressMessages(kits_dif(resp, "sex"))
  expect_identical(lr(x)$groups$group, c("boy", "girl"))
  expect_identical(lr(x)$split, "sex")
  g <- ref$split3
  g[1:10] <- NA
  expect_message(y <- kits_dif(ref$X, g), "10 person-wave", class = "kits_message_data")
  expect_identical(sum(lr(y)$groups$n), 990L)
  dg <- attr(y, "diagnostics")
  expect_identical(dg$value[dg$check == "person-waves without a group"], 10)
})

# ---- 輸入錯誤 -------------------------------------------------------------------------------------

test_that("錯誤的輸入：kits_error_input；只有一組：kits_error_data", {
  X <- ref$X
  expect_error(kits_dif(X, "wave"), "kits_resp", class = "kits_error_input")
  expect_error(kits_dif(resp, "nope"), "no column", class = "kits_error_input")
  expect_error(kits_dif(X, 1:3), class = "kits_error_input")
  expect_error(kits_dif(X, NA_character_), class = "kits_error_input")
  expect_error(kits_dif(X, list(1)), class = "kits_error_input")
  expect_error(kits_dif(X, "raw", unobserved = "some"), class = "kits_error_input")
  expect_error(kits_dif(X, "raw", splitcr = "median"), "splitcr", class = "kits_error_input")
  item <- kits_item(X, scale = "lr")
  other <- kits_item(X[1:500, ], scale = "lr")
  expect_error(kits_dif(X, "raw", item = coef(item)), class = "kits_error_input")
  expect_error(kits_dif(X, "raw", item = other), "other data", class = "kits_error_input")
  persons <- data.frame(visit_id = rownames(X)[-1], wave_m = 1)
  expect_error(kits_dif(X, "wave", persons = persons), "no row", class = "kits_error_input")
  expect_error(kits_dif(X, rep("a", nrow(X))), "two groups", class = "kits_error_data")
  expect_error(kits_dif(X, rep(1:2, 500), ref = "3"), "groups are 1, 2", class = "kits_error_input")
})

test_that("persons 有重複的 visit_id：kits_error_data（和 kits_subj()、kits_resp() 相同）", {
  X <- ref$X
  persons <- data.frame(visit_id = rownames(X)[c(1, seq_len(nrow(X)))], wave_m = 1)
  expect_error(kits_dif(X, "wave", persons = persons), "repeated visit_id", class = "kits_error_data")
})

test_that("item 有固定參數：全部一起的模型重新校準（訊息），結果相同", {
  X <- ref$X
  item <- kits_item(X, scale = "lr")
  fixed <- suppressMessages(kits_item(X, fixed_params = coef(item)[1:2, ]))
  expect_message(x <- kits_dif(X, "raw", fixed), class = "kits_message_fit")
  expect_equal(lr(x)$LR, lr(kits_dif(X, "raw", item))$LR, tolerance = 1e-8)
})

# ---- print、summary、plot ------------------------------------------------------------------------

test_that("print：組、LR 檢定和最大的步驟差，兩組時再印每題的 DIF；summary", {
  x <- kits_dif(ref$X, "raw")
  out <- capture.output(print(x))
  expect_match(out[1], "^<kits_dif> groups by raw score")
  expect_identical(out[2], "Groups: low (n = 561), high (n = 439)")
  expect_identical(out[3], "LR test (Andersen): LR = 4.82, df = 8, p = 0.777")
  expect_match(out[4], "^Largest \\|step difference\\|: .* logits")
  expect_true(any(grepl("^Item DIF in logits: low minus high", out)))
  expect_true(any(grepl("^Overall DIF test: chi-square", out)))
  expect_true(any(grepl("^Checks: OK", out)))
  big <- max(abs(part(x, "diff")$diff))
  dg <- attr(x, "diagnostics")
  expect_equal(dg$value[dg$check == "largest |step difference|"], big)
  s <- capture.output(print(summary(x, n = 3)))
  expect_true(any(grepl("Largest step differences \\(3 of 9", s)))
  expect_true(any(grepl("^Diagnostics:", s)))
  expect_s3_class(kits_check(x, all = TRUE), "kits_check")
  # 三組：沒有每題的 DIF（0 列，沒有 overall、groups、se、left_out）
  x3 <- kits_dif(ref$X, ref$split3)
  expect_identical(nrow(x3), 0L)
  expect_null(attr(x3, "overall"))
  expect_null(attr(x3, "groups"))
  out3 <- capture.output(print(x3))
  expect_false(any(grepl("Item DIF", out3)))
  expect_error(confint(x3), "two groups", class = "kits_error_input")
})

test_that("兩組時每題的 DIF 和 LR 檢定在同一個結果裡；ref 沒給時 boy 或第一組", {
  X <- ref$X
  g <- rep(c("M", "F"), 500)
  x <- kits_dif(X, g)
  expect_identical(attr(x, "groups")$ref, "F")
  expect_identical(lr(x)$groups$group, c("F", "M"))
  expect_equal(kits_dif(X, g, ref = "M")$dif, -x$dif, tolerance = 1e-12)
  sex <- rep(c("girl", "boy"), 500)
  expect_identical(attr(kits_dif(X, sex), "groups")$ref, "boy")
  # 一組不能校準的題（只有 0 分）：這組不用它，每題的 DIF 不比它（舊版停下來）
  X[g == "F", "a1"] <- 0
  expect_message(y <- kits_dif(X, g), "a1", class = "kits_message_data")
  expect_false("a1" %in% y$item_id)
  expect_true("a1" %in% attr(y, "left_out"))
})

test_that("plot：兩組是組對組（identity 線、+/- 1.96 SE），多組是每組對全部", {
  x <- kits_dif(ref$X, "raw")
  p <- plot(x, type = "steps")
  expect_s3_class(p, "plotly")
  types <- vapply(p$x$data, \(tr) if (is.null(tr$mode)) "" else tr$mode, "")
  expect_true("lines" %in% types)
  pts <- Filter(\(tr) identical(tr$mode, "markers"), p$x$data)
  expect_identical(sum(vapply(pts, \(tr) length(tr$x), 1L)), nrow(part(x, "diff")))
  expect_true(all(vapply(pts, \(tr) !is.null(tr$error_x) && !is.null(tr$error_y), TRUE)))
  p3 <- plot(kits_dif(ref$X, ref$split3))
  expect_s3_class(p3, "plotly")
  expect_true(all(c("xaxis3", "yaxis3") %in% names(p3$x$layout)))
  expect_s3_class(plot(x, type = "steps", items = c("a1", "a2")), "plotly")
  expect_error(plot(x, type = "steps", items = "zz"), class = "kits_error_input")
  expect_error(plot(kits_dif(ref$X, ref$split3), type = "items"), "two groups", class = "kits_error_input")
})

test_that("plot(type = \"steps\", items = )：三組、有一組在選的題裡沒有步驟時也畫得出來", {
  X <- ref$X
  g <- rep(c("A", "B", "C"), length.out = nrow(X))
  X[g == "A", "a5"] <- 0   # A 組 a5 只有 0 分：沒有步驟
  x <- suppressMessages(kits_dif(X, g))
  expect_false(any(part(x, "steps")$group == "A" & part(x, "steps")$item_id == "a5"))
  p <- plot(x, items = "a5")
  expect_s3_class(p, "plotly")
  expect_true(all(c("xaxis3", "yaxis3") %in% names(p$x$layout)))
})
