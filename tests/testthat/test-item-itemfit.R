# kits 0.23.0：kits_item() 算題目適配（item$itemfit），和 kits_itemfit(data, item) 同一個表；
# 印出一行、診斷一列（note，不是 warn）、plot(item, type = "fit")、算不出來時 NULL、可以分享、kits() 也有

betas <- list(a1 = c(-1, 0), a2 = c(0, 1), a3 = 0.5, a4 = -0.5, a5 = c(-2, -1, 1), a6 = -1.5, a7 = 1.5,
              a8 = c(-0.5, 0.5))
X <- kits_sim(seq(-3, 3, length.out = 400), betas, seed = 2)
item <- kits_item(X, scale = "fit demo")

# 一題是純雜訊：它的 infit、outfit 在範圍外（其他題在範圍內）
noisy_X <- X
set.seed(1)
noisy_X[, "a3"] <- sample(0:1, 400, replace = TRUE)
noisy <- kits_item(noisy_X, scale = "noisy")

# info 的呼叫是校準的（kits_item(...)），不是 kits_itemfit(...)；其他都要完全相同
without_call <- function(fit) {
  info <- attr(fit, "info")
  info$call <- NULL
  attr(fit, "info") <- info
  fit
}

expect_same_itemfit <- function(item, data) {
  fit <- kits_itemfit(data, item)
  expect_identical(without_call(item$itemfit), without_call(fit))
  expect_identical(class(item$itemfit), class(fit))
  expect_identical(names(item$itemfit),
                   c("item_id", "n", "pvalue", "rir", "infit", "outfit", "flag"))
  expect_identical(attr(item$itemfit, "info")$call, item$info$call)
}

# 只有欄（不比 info：每個校準的呼叫和量尺不同）
columns <- function(fit) stats::setNames(lapply(names(fit), \(nm) fit[[nm]]), names(fit))

# kits 的 print 方法（kitsPreview 載入時會蓋掉一些 print 方法）
pr <- function(x, ...) {
  fun <- get(paste0("print.", class(x)[1]), envir = asNamespace("kits"))
  fun(x, ...)
}
printed <- function(x) utils::capture.output(pr(x))

# 印出的 Item fit 那幾行（到 until 之前），接成一行
fit_text <- function(out, until = "^Next:") {
  lines <- out[seq(grep("^Item fit", out), grep(until, out) - 1)]
  gsub(" +", " ", paste(lines, collapse = " "))
}

test_that("item$itemfit 和 kits_itemfit(data, item) 相同：矩陣、寬的、長的 data.frame、kits_resp", {
  expect_same_itemfit(item, X)
  wide <- data.frame(visit_id = rownames(X), X)
  item_wide <- kits_item(wide, scale = "wide")
  expect_same_itemfit(item_wide, wide)
  long <- data.frame(visit_id = rep(rownames(X), ncol(X)), item_id = rep(colnames(X), each = nrow(X)),
                     item_score = as.vector(X))
  item_long <- kits_item(long, scale = "long")
  expect_same_itemfit(item_long, long)
  resp <- kitm03$responses |>
    dplyr::filter(domain == "lan", wave_m <= 6) |>
    kits_resp() |>
    suppressMessages()
  item_resp <- kits_item(resp, scale = "lan M03")
  expect_same_itemfit(item_resp, resp)
  # 同樣的資料：四種輸入的數字相同
  expect_identical(columns(item_wide$itemfit), columns(item$itemfit))
  expect_identical(columns(item_long$itemfit), columns(item$itemfit))
  # 元素的名字：fit 還是最佳化的結果
  expect_s3_class(item$fit, "kits_fit")
  expect_identical(names(item), c("thresholds", "fit", "stats", "n_persons", "left_out", "design", "scale",
                                  "itemfit", "diagnostics", "info"))
})

test_that("固定參數（含沒放進校準的題）時也算，和 kits_itemfit() 相同", {
  betas_new <- c(betas, list(a9 = c(0.5, 1.5)))
  X_new <- kits_sim(seq(-2, 4, length.out = 400), betas_new, prefix = "n", seed = 2)
  item_new <- kits_item(X_new, fixed_params = coef(item))
  expect_same_itemfit(item_new, X_new)
  # a4 在新的資料裡每個人都通過：不放進校準，題目適配照樣算（照固定的參數）
  X_left <- X_new
  X_left[, "a4"] <- 1
  item_left <- suppressMessages(kits_item(X_left, fixed_params = coef(item)))
  expect_identical(item_left$left_out, "a4")
  expect_same_itemfit(item_left, X_left)
  expect_true("a4" %in% item_left$itemfit$item_id)
})

test_that("印出：範圍內和範圍外各一行，在 Next: 之前，每行不超過 100 字", {
  out <- printed(item)
  line <- "Item fit (conditional): infit and outfit within 0.6-1.4 for all 8 items."
  expect_true(line %in% out)
  expect_identical(which(out == line) + 1L, which(startsWith(out, "Next:")))
  out <- printed(noisy)
  expect_identical(fit_text(out), paste0("Item fit (conditional): 1 of 8 items outside 0.6-1.4 ",
                                         "(1 underfit: a3); see x$itemfit, plot(x, type = \"fit\")."))
  # 很多題在範圍外（7 題反向計分）：最多 5 個題號，再加 ...；太長時斷行，每行最多 100 字
  b <- stats::setNames(as.list(seq(-2, 2, length.out = 24)), sprintf("motor_item_%02d", 1:24))
  many_X <- kits_sim(seq(-3, 3, length.out = 400), b, seed = 3)
  for (j in 1:7) many_X[, j] <- 1 - many_X[, j]
  many <- kits_item(many_X, scale = "many")
  out <- printed(many)
  fit_lines <- out[seq(grep("^Item fit", out), grep("^Next:", out) - 1)]
  expect_true(all(nchar(out) <= 100))
  expect_gt(length(fit_lines), 1)
  head <- "^Item fit \\(conditional\\): 7 of 24 items outside 0.6-1.4 \\(7 underfit: motor_item_01, "
  expect_match(fit_lines[1], head)
  all_text <- paste(fit_lines, collapse = " ")
  expect_identical(lengths(regmatches(all_text, gregexpr("motor_item_", all_text))), 5L)
  expect_match(all_text, "motor_item_05, ...); see x$itemfit, plot(x, type = \"fit\").", fixed = TRUE)
  # 範圍是 plot() 的預設（只寫一次）
  expect_identical(kits:::itemfit_range(), eval(formals(kits:::plot.kits_itemfit)$range))
})

test_that("診斷：item fit 一列；範圍外是 note（說出題號），範圍內是 ok；從來不是 warn", {
  dg <- item$diagnostics
  expect_identical(dg$check[nrow(dg)], "item fit")
  expect_identical(dg$status[dg$check == "item fit"], "ok")
  expect_identical(dg$value[dg$check == "item fit"], 0)
  dn <- noisy$diagnostics
  row <- dn[dn$check == "item fit", ]
  expect_identical(row$status, "note")
  expect_identical(row$value, 1)
  msg <- "1 of 8 items outside 0.6-1.4 (conditional infit or outfit): 1 underfit (a3)"
  expect_match(row$message, msg, fixed = TRUE)
  chk <- kits_check(noisy)
  expect_true("item fit" %in% chk$check)
  expect_identical(chk$status[chk$check == "item fit"], "note")
  expect_false("item fit" %in% kits_check(item)$check)
  expect_output(print(summary(noisy)), "item fit")
  expect_false(any(dn$status == "warn"))
})

test_that("plot(item, type = \"fit\") 和 plot(kits_itemfit()) 同一張圖；items = 選題號；預設的 type 不變", {
  data_of <- function(p) plotly::plotly_build(p)$x$data
  layout_of <- function(p) plotly::plotly_build(p)$x$layout[c("title", "shapes", "xaxis", "yaxis")]
  p_item <- plot(noisy, type = "fit")
  p_fit <- plot(kits_itemfit(noisy_X, noisy))
  expect_s3_class(p_item, "plotly")
  expect_equal(data_of(p_item), data_of(p_fit))
  expect_equal(layout_of(p_item), layout_of(p_fit))
  p_two <- plot(noisy, type = "fit", items = c("a1", "a3"))
  expect_equal(data_of(p_two), data_of(plot(kits_itemfit(noisy_X, noisy), items = c("a1", "a3"))))
  expect_identical(sort(unlist(lapply(data_of(p_two), \(d) as.character(d$x)))),
                   c("a1", "a1", "a3", "a3"))
  expect_error(plot(noisy, type = "fit", items = "zz"), class = "kits_error_input")
  expect_identical(eval(formals(kits:::plot.kits_item)$type)[1], "category")
})

test_that("算不出來時 item$itemfit 是 NULL、診斷有 note，校準照常完成", {
  Y <- cbind(a1 = rep(0:1, 50), a2 = rep(0:1, 50))   # 每個人都是 0 分或滿分
  rownames(Y) <- paste0("p", 1:100)
  fixed <- data.frame(item_id = c("a1", "a2"), item_score = 1, beta = c(-0.5, 0.5))
  expect_no_error(it <- kits_item(Y, fixed_params = fixed, scale = "null case"))
  expect_true("itemfit" %in% names(it))
  expect_null(it$itemfit)
  row <- it$diagnostics[it$diagnostics$check == "item fit", ]
  expect_identical(row$status, "note")
  expect_match(row$message, "not computed .*zero or a perfect score")
  expect_true("Item fit: not computed; see summary(x)." %in% printed(it))
  expect_error(plot(it, type = "fit"), "no item fit", class = "kits_error_input")
  expect_error(kits_itemfit(Y, it), class = "kits_error_fit")   # kits_itemfit() 照舊停下來
  # 舊的物件（沒有 itemfit 這個元素）：不印那一行
  old <- item
  old$itemfit <- NULL
  expect_false(any(startsWith(printed(old), "Item fit")))
})

test_that("可以分享：item$itemfit 只有一題一列的彙總，沒有一個人次一列的欄或屬性", {
  n_visits <- nrow(X)
  fit <- item$itemfit
  expect_equal(nrow(fit), length(betas))
  expect_true(all(lengths(as.list(fit)) != n_visits))
  leaf_lengths <- function(v) if (is.list(v)) unlist(lapply(v, leaf_lengths)) else length(v)
  for (a in names(attributes(fit))) {
    v <- attr(fit, a)
    if (is.data.frame(v)) expect_false(nrow(v) == n_visits, info = a)
    expect_false(any(leaf_lengths(v) == n_visits), info = a)
  }
  expect_false(any(rownames(X) %in% unlist(attributes(fit), use.names = FALSE)))
  expect_false(any(rownames(X) %in% unlist(as.list(fit), use.names = FALSE)))
})

test_that("kits()：x$item$itemfit、印出（指向 x$item）、診斷", {
  k <- kits(noisy_X, scale = "noisy")
  expect_identical(without_call(k$item$itemfit), without_call(kits_itemfit(noisy_X, k$item)))
  expect_identical(columns(k$item$itemfit), columns(noisy$itemfit))
  out <- utils::capture.output(print(k))
  expect_identical(fit_text(out, "^Abilities"),
                   paste0("Item fit (conditional): 1 of 8 items outside 0.6-1.4 (1 underfit: a3); ",
                          "see x$item$itemfit, plot(x$item, type = \"fit\")."))
  expect_true(all(nchar(out) <= 100))
  expect_lt(grep("^Item fit", out), grep("^Abilities", out))
  expect_identical(sum(startsWith(out, "Next:")), 1L)
  expect_identical(k$diagnostics$status[k$diagnostics$check == "item: item fit"], "note")
  expect_false(any(startsWith(utils::capture.output(print(kits(X, scale = "clean"))), "Checks: 1 warning")))
})

# 一張改過 infit、outfit 的表（flag 照新的值重算），給印出和診斷的測試
with_fit <- function(item, infit, outfit) {
  fit <- item$itemfit
  fit$infit <- infit
  fit$outfit <- outfit
  fit$flag <- kits:::itemfit_flag(infit, outfit)
  item$itemfit <- fit
  item
}

test_that("flag：underfit 優先、不然 overfit、其他 ok；剛好在 0.6、1.4 是 ok；兩個都 NA 時 NA", {
  infit <- c(1.4, 0.6, 1.41, 0.59, 1.5, 0.5, 1.0, NA, NA, 1.0, NaN)
  outfit <- c(1.4, 0.6, 1.0, 1.0, 0.5, 1.5, 0.59, NA, 1.5, NA, NaN)
  want <- c("ok", "ok", "underfit", "overfit", "underfit", "underfit", "overfit", NA, "underfit",
            "ok", NA)
  expect_identical(kits:::itemfit_flag(infit, outfit), want)
  # 範圍是 plot() 的預設；給別的範圍時照那個範圍
  expect_identical(kits:::itemfit_flag(1.45, 1, c(0.5, 1.5)), "ok")
  # 表裡的欄：在 outfit 後面，照同一個規則
  fit <- noisy$itemfit
  expect_identical(names(fit)[6:7], c("outfit", "flag"))
  expect_type(fit$flag, "character")
  expect_identical(fit$flag, kits:::itemfit_flag(fit$infit, fit$outfit))
  expect_identical(fit$flag[fit$item_id == "a3"], "underfit")
  expect_identical(unique(item$itemfit$flag), "ok")
  # kits_itemfit()、item$itemfit、kits() 都有，而且相同（除了呼叫）
  expect_identical(without_call(noisy$itemfit), without_call(kits_itemfit(noisy_X, noisy)))
  k <- kits(noisy_X, scale = "noisy")
  expect_identical(k$item$itemfit$flag, fit$flag)
  expect_match(utils::capture.output(kits:::print.kits_table(fit))[1], "infit +outfit +flag$")
})

test_that("印出和診斷：範圍外的題分成 underfit 和 overfit，underfit 的題號在前", {
  both <- with_fit(item, infit = c(0.5, 1.5, 1, 0.55, 1, 1, 1, 1.45),
                   outfit = c(1, 1, 1, 1, 1, 1.6, 0.5, 0.5))
  # a2、a6、a8 是 underfit（a8 兩個都有：underfit 優先）；a1、a4、a7 是 overfit
  expect_identical(both$itemfit$flag, c("overfit", "underfit", "ok", "overfit", "ok", "underfit",
                                        "overfit", "underfit"))
  line <- paste0("Item fit (conditional): 6 of 8 items outside 0.6-1.4 ",
                 "(3 underfit, 3 overfit: a2, a6, a8, a1, a4, ...); ",
                 "see x$itemfit, plot(x, type = \"fit\").")
  out <- printed(both)
  expect_identical(fit_text(out), line)
  expect_true(all(nchar(out) <= 100))
  # 只有 overfit
  low <- with_fit(item, infit = c(0.5, rep(1, 7)), outfit = rep(1, 8))
  expect_identical(fit_text(printed(low)),
                   paste0("Item fit (conditional): 1 of 8 items outside 0.6-1.4 (1 overfit: a1); ",
                          "see x$itemfit, plot(x, type = \"fit\")."))
  # 全部在範圍內的句子不變
  line <- "Item fit (conditional): infit and outfit within 0.6-1.4 for all 8 items."
  expect_true(line %in% printed(item))
  # 診斷：還是 note，說出兩種的題數和題號
  row <- kits:::itemfit_row(list(table = both$itemfit))
  expect_identical(row$status, "note")
  expect_identical(row$value, 6)
  msg <- paste0("6 of 8 items outside 0.6-1.4 (conditional infit or outfit): ",
                "3 underfit (a2, a6, a8), 3 overfit (a1, a4, a7); see x$itemfit")
  expect_match(row$message, msg, fixed = TRUE)
  row <- kits:::itemfit_row(list(table = low$itemfit))
  expect_match(row$message, "(conditional infit or outfit): 1 overfit (a1); see", fixed = TRUE)
  # 舊的物件（表裡沒有 flag 欄）：照 infit、outfit 算，印出相同
  old <- both
  old$itemfit$flag <- NULL
  expect_identical(printed(old), out)
})

test_that("plot()：點的顏色是 flag（underfit 橘、overfit 藍、ok 灰），照圖的範圍", {
  traces <- plotly::plotly_build(plot(noisy, type = "fit"))$x$data
  names_of <- vapply(traces, \(d) d$name, character(1))
  expect_identical(names_of, c("infit: underfit", "infit: ok", "outfit: underfit", "outfit: ok"))
  expect_identical(vapply(traces, \(d) d$marker$color, character(1)),
                   c("#D55E00", "#7F7F7F", "#D55E00", "#7F7F7F"))
  expect_identical(vapply(traces, \(d) d$marker$symbol, character(1)),
                   c("circle", "circle", "square", "square"))
  expect_identical(as.character(traces[[1]]$x), "a3")
  wide <- plotly::plotly_build(plot(noisy$itemfit, range = c(0.1, 10)))$x$data
  expect_identical(vapply(wide, \(d) d$name, character(1)), c("infit: ok", "outfit: ok"))
})

test_that("type：conditional（預設，item$itemfit 也是）和 wle；記在 attr 和 info；印出第一行、圖的標題", {
  wle <- kits_itemfit(X, item, type = "wle")
  cond <- kits_itemfit(X, item)
  expect_identical(attr(cond, "type"), "conditional")
  expect_identical(attr(wle, "type"), "wle")
  expect_identical(attr(item$itemfit, "type"), "conditional")
  expect_identical(kits_info(cond)$settings$type, "conditional")
  expect_identical(names(wle), names(cond))
  expect_identical(columns(cond)[c("item_id", "pvalue", "rir")], columns(wle)[c("item_id", "pvalue", "rir")])
  expect_false(isTRUE(all.equal(cond$infit, wle$infit)))
  expect_error(kits_itemfit(X, item, type = "eap"), "type", class = "kits_error_input")
  out <- utils::capture.output(kits:::print.kits_itemfit(cond))
  expect_identical(out[1], paste0("Item fit (conditional): infit and outfit given the raw score ",
                                  "within each item pattern"))
  expect_match(utils::capture.output(kits:::print.kits_itemfit(wle))[1], "^Item fit \\(wle\\)")
  expect_identical(attr(cond[1:3, ], "type"), "conditional")   # 選列時留著
  title <- plotly::plotly_build(plot(cond))$x$layout$title$text
  expect_match(title, "Item fit (conditional): mean squares", fixed = TRUE)
})

test_that("conditional：給定原始分，各題的條件期望值加起來是原始分；一題的題本、0 分、滿分不算", {
  b <- list(i1 = c(-9, 0, 8), i2 = c(-3, 7), i3 = 5, i4 = c(-8, -7, -6))
  tab <- kits:::cond_table(kits:::set_key(names(b)), b)
  expect_equal(as.vector(tapply(tab$E, tab$raw, sum)), 0:9, tolerance = 1e-12)
  expect_true(all(tab$W >= 0))
  expect_false(any(tab$free[tab$raw %in% c(0, 9)]))
  expect_true(all(tab$free[tab$raw %in% 1:8 & tab$item_id == "i1"]))
  # 原始分 1 時，3 分的題只能是 0 或 1
  expect_equal(sum(tab$free[tab$raw == 1]), 4)
  one <- kits:::cond_table("i1", b)
  expect_false(any(one$free))
  # n：只算條件變異數不是 0 的作答（0 分和滿分的人次不算）
  raw <- rowSums(X)
  inner <- raw > 0 & raw < sum(lengths(betas))
  expect_identical(item$itemfit$n, rep(sum(inner), ncol(X)))
})

test_that("conditional：短題本時沒有 WLE 的向下偏誤（PCM 成立時平均約 1）", {
  b <- stats::setNames(as.list(seq(-1.5, 1.5, length.out = 4)), paste0("s", 1:4))
  S <- kits_sim(stats::rnorm(3000, 0, 1.5), b, seed = 5)
  short <- kits_item(S)
  wle <- kits_itemfit(S, short, type = "wle")
  expect_lt(abs(mean(short$itemfit$infit) - 1), 0.03)
  expect_lt(mean(wle$infit), 0.9)
})

test_that("conditional：完整的資料時和 iarm::out_infit() 相同（存下來的參考值）", {
  # fixtures/fit_ref/make_fit_ref.R：eRm 1.0.10 的 PCM、iarm 0.4.3 的 out_infit()；kits 不依賴它們
  ref <- readRDS(test_path("fixtures", "fit_ref", "cond_fit.rds"))
  # 同樣的參數（eRm 的門檻）：完全相同
  fit <- kits_itemfit(ref$data, ref$parms)
  expect_lt(max(abs(fit$infit - ref$infit)), 1e-10)
  expect_lt(max(abs(fit$outfit - ref$outfit)), 1e-10)
  # kits 自己的 CML：只差 eRm 的收斂精度
  own <- kits_item(ref$data)$itemfit
  expect_lt(max(abs(own$infit - ref$infit)), 1e-4)
  expect_lt(max(abs(own$outfit - ref$outfit)), 1e-4)
})
