# 0.16.0：三份審查（專家、學生、Shiny 計畫的隱私）找到的問題，每一個一個測試。
# 每個測試在 0.15.0 會失敗，修正後通過；重現的程式來自審查（review_20261027/）。

betas7 <- list(a1 = -2, a2 = -1, a3 = c(-0.5, 0.5), a4 = 0, a5 = c(0.5, 1.5), a6 = 2, a7 = 3)
young <- kits_sim(rnorm(500, 0), betas7[1:5], prefix = "y", seed = 1)
old <- kits_sim(rnorm(500, 2), betas7[3:7], prefix = "o", seed = 2)
item_y <- kits_item(young, scale = "young")
item_o <- kits_item(old, scale = "old")
link_oy <- kits_linking(item_o, item_y)
anchors <- data.frame(item_id = paste0("a", 1:5), tau = c(40, 47, 51, 55, 60))

# ---- A1 隱私 -------------------------------------------------------------------------------------------------

test_that("A1 kits_info()$call 不帶資料：do.call()、map() 給的值寫 <data>；名字和設定照寫", {
  via_do_call <- kits_info(do.call(kits_item, list(young, scale = "young")))$call
  expect_equal(via_do_call, "kits_item(<data>, scale = \"young\")")
  expect_false(grepl("[0-9]{2}", via_do_call))
  via_map <- purrr::map(list(young), \(x) kits_item(x, scale = "s"))[[1]]
  expect_false(grepl("structure|c\\(", kits_info(via_map)$call))
  expect_equal(kits_info(item_y)$call, "kits_item(young, scale = \"young\")")
  expect_equal(kits_info(kits_linking(item_o, item_y, n_boot = 20))$call, "kits_linking(item_o, item_y, n_boot = 20)")
  # a single child's score and age given by value are data
  norm <- data.frame(age_m = c(6, 12), mu = c(40, 60), sigma = 5, epsilon = 0, delta = 1)
  s1 <- kits_subj(young, item_y)
  expect_equal(kits_info(do.call(kits_subj, list(young, item_y)))$call, "kits_subj(<data>, <data>)")
  # persons given by value
  p <- data.frame(visit_id = rownames(young), wave_m = 12)
  expect_match(kits_info(do.call(kits_subj, list(young, item_y, persons = p)))$call, "persons = <data>", fixed = TRUE)
  # kits_check() labels a value by its position
  expect_equal(unique(do.call(kits_check, list(item_y, all = TRUE))$object), "object 1")
  expect_true(is.numeric(norm$mu))
  expect_true(is.data.frame(s1))
})

test_that("A1 預設的量尺標籤不帶資料：名字，否則 \"scale 1\"（有訊息）", {
  expect_equal(suppressMessages(kits_item(young))$scale, "young")
  expect_message(lab <- do.call(kits_item, list(young))$scale, "scale 1", class = "kits_message_scale")
  expect_equal(lab, "scale 1")
  expect_message(kits(young[1:300, ]), "scale 1", class = "kits_message_scale")
  # two unlabelled calibrations cannot be linked by mistake
  a <- suppressMessages(kits_item(young[1:250, ]))
  b <- suppressMessages(kits_item(young[251:500, ]))
  expect_error(kits_linking(a, b), "scale 1", class = "kits_error_scale")
})

test_that("A1 訊息和錯誤預設不給 visit_id，只給人數；kits.show_ids = TRUE 時才列出", {
  Y <- young
  rownames(Y)[2] <- rownames(Y)[1]
  e <- tryCatch(kits_item(Y), error = conditionMessage)
  expect_match(e, "1 repeated person id")
  expect_false(grepl(rownames(young)[1], e, fixed = TRUE))
  p <- data.frame(visit_id = rownames(young)[-(1:3)], wave_m = 12)
  e <- tryCatch(kits_subj(young, item_y, persons = p), error = conditionMessage)
  expect_match(e, "^3 person-wave\\(s\\) of the data have no row in `persons`")
  expect_false(grepl("y000", e))
  e <- tryCatch(kits_dif(young, stats::setNames(rep("boy", 499), rownames(young)[-1])), error = conditionMessage)
  expect_false(grepl("y000", e))
  op <- options(kits.show_ids = TRUE)
  on.exit(options(op))
  expect_match(tryCatch(kits_item(Y), error = conditionMessage), rownames(young)[1], fixed = TRUE)
})

test_that("A1 plot(norm, scores = )：預設分組畫（沒有個人的點，少於 5 人的格不畫）；persons = \"points\" 才畫每個人", {
  set.seed(1)
  age <- runif(600, 3, 36)
  D <- 20 + 15 * log(age) + rnorm(600, 0, 3)
  norm <- kits_norms(D, age, ages = seq(3, 36, by = 3), degrees = c(3, 2, 1, 1))
  children <- data.frame(theta = D, age_m = age)
  built <- plotly::plotly_build(plot(norm, scores = children))$x$data
  types <- vapply(built, \(t) if (is.null(t$type)) "" else t$type, "")
  expect_true("heatmap" %in% types)
  expect_false(any(vapply(built, \(t) identical(t$mode, "markers"), TRUE)))
  heat <- built[[which(types == "heatmap")]]
  z <- unlist(heat$z)
  expect_true(all(is.na(z) | z >= 5))
  expect_false(any(D %in% unlist(heat$y)))
  pts <- plotly::plotly_build(plot(norm, scores = children, persons = "points"))$x$data
  expect_true(any(vapply(pts, \(t) identical(t$mode, "markers") && length(t$x) == 600, TRUE)))
  expect_error(plot(norm, scores = children, persons = "dots"), class = "kits_error_input")
})

# ---- A2、A3 疊表 --------------------------------------------------------------------------------------------

test_that("A2 bind_rows()：量尺不同停下來（kits_error_scale），和 rbind() 一樣；量尺相同時留著類別和量尺", {
  sy <- kits_subj(young, item_y)
  so <- kits_subj(old, item_o)
  expect_error(dplyr::bind_rows(sy, so), "young and old", class = "kits_error_scale")
  expect_error(rbind(sy, so), "young and old", class = "kits_error_scale")
  both <- quiet_records(dplyr::bind_rows(sy, kits_convert(so, link_oy)))
  expect_s3_class(both, "kits_subj")
  expect_equal(attr(both, "scale"), "young")
  expect_equal(nrow(both), 1000)
  expect_equal(both$theta, rbind(sy, kits_convert(so, link_oy))$theta)
  # the wrong conversion of the student's example now stops
  d_link <- kits_linking(item_y, anchors)
  expect_error(kits_convert(dplyr::bind_rows(sy, so), d_link), class = "kits_error_scale")
  # thresholds tables too
  expect_error(dplyr::bind_rows(coef(item_y), coef(item_o)), class = "kits_error_scale")
  expect_equal(attr(dplyr::bind_rows(coef(item_y), coef(item_y)), "scale"), "young")
})

test_that("A2 bind_rows()、rbind() 混進一般的 data.frame 或 tibble：警告（kits_warning_scale）", {
  sy <- kits_subj(young, item_y)
  expect_warning(out <- quiet_records(dplyr::bind_rows(sy, as.data.frame(sy))),
                 class = "kits_warning_scale")
  expect_s3_class(out, "kits_subj")
  expect_equal(attr(out, "scale"), "young")
  expect_warning(quiet_records(dplyr::bind_rows(sy, dplyr::as_tibble(as.data.frame(sy)))),
                 class = "kits_warning_scale")
  expect_warning(dplyr::bind_rows(as.data.frame(sy), sy), class = "kits_warning_scale")
  expect_warning(rbind(sy, as.data.frame(sy)), class = "kits_warning_scale")
})

test_that("A3 疊起來的表：summary() 和 kits_check() 都可以用；kits_subj 重算診斷，其他的表留著各表的診斷", {
  w <- rep(c(6, 12), length.out = 500)
  sy <- kits_subj(young, item_y, persons = data.frame(visit_id = rownames(young), wave_m = w))
  so <- kits_convert(kits_subj(old, item_o,
                               persons = data.frame(visit_id = rownames(old),
                                                    wave_m = w + 30)), link_oy)
  for (both in list(rbind(sy, so), quiet_records(dplyr::bind_rows(sy, so)))) {
    out <- capture.output(print(summary(both)))
    expect_false(any(grepl("No diagnostics recorded", out)))
    expect_true(any(grepl("reliability, wave 36", out)))
    chk <- kits_check(both, all = TRUE)
    expect_true(any(grepl("share of perfect scores", chk$check)))
  }
  tab_y <- subj_table_of(item_y, young, wave_m = w)
  tab_o <- kits_convert(subj_table_of(item_o, old, wave_m = w + 30), link_oy)
  raw_table <- rbind(tab_y, tab_o)
  chk <- kits_check(tab_y, tab_o, raw_table, all = TRUE)
  expect_true(all(c("tab_y", "tab_o", "raw_table") %in% chk$object))
  expect_true(any(grepl("^table 2: ", chk$check[chk$object == "raw_table"])))
  # after bind_rows() the other tables have none: a note, not an error (no purrr backtrace)
  chk <- kits_check(dplyr::bind_rows(tab_y, tab_o))
  expect_equal(chk$status, "note")
  expect_match(chk$message, "check the tables before stacking")
  # not a kits result: a warning row, and the other objects are still checked
  chk <- kits_check(item_y, "hello", all = TRUE)
  expect_true(any(chk$object == "item_y"))
  expect_equal(chk$status[chk$object == "object 2"], "warn")
})

# ---- A4 DIF 的整體檢定 -------------------------------------------------------------------------------------

test_that("A4 DIF 的整體檢定：有題沒有比較時（只有一組答），檢定和原點無關，沒有 DIF 時的拒絕率接近 5%", {
  betas <- list(a1 = c(-1, 0), a2 = c(0, 1), a3 = 0.5, a4 = -0.5, a5 = 0, a6 = 1, a7 = -1)
  one <- function(r, left_out = TRUE) {
    X <- rbind(kits_sim(stats::rnorm(400), betas, prefix = "b", seed = r), kits_sim(stats::rnorm(400), betas, prefix = "g", seed = 1000 + r))
    if (left_out) X <- cbind(X, a8 = c(kits_sim(stats::rnorm(400), list(a8 = 3), prefix = "b", seed = 5000 + r)[, 1], rep(NA, 400)))
    suppressMessages(kits_dif(X, rep(c("boy", "girl"), each = 400)))
  }
  set.seed(1)
  d1 <- one(1)
  expect_equal(attr(d1, "left_out"), "a8")
  expect_gt(attr(d1, "overall")$p, 1e-3)        # 0.15.0: chi-square 67, p = 1.7e-11, although no item has DIF
  p <- vapply(1:40, \(r) attr(one(r), "overall")$p, numeric(1))
  expect_lt(mean(p < 0.05), 0.15)               # 0.15.0: 1 (every replicate rejected)
  # the statistic does not depend on the origin of either calibration
  st <- list(list(beta = c(-1, 0.2, 0.8), V = diag(c(0.02, 0.03, 0.02))), list(beta = c(-0.9, 0.1, 0.7), V = diag(0.02, 3)))
  shifted <- st
  shifted[[2]]$beta <- shifted[[2]]$beta + 5
  expect_equal(kits:::dif_overall(shifted)$stat, kits:::dif_overall(st)$stat)
  # with no item left out, the dif values and the test are those of 0.15.0 (dexter's DIF_overall): see test-against-kitsPreview.R
})

# ---- A5 CML ------------------------------------------------------------------------------------------------

test_that("A5 長的題本、固定題離原點很遠：ESF 不下溢，校準算得出來", {
  set.seed(4)
  bt <- stats::setNames(lapply(1:40, function(i) 6 + c(-1, 0, 1) + stats::rnorm(1, 0, 0.5)), sprintf("i%02d", 1:40))
  W <- kits_sim(stats::rnorm(500, 6, 1.5), bt, seed = 5)
  fixed <- data.frame(item_id = rep(names(bt)[1:5], each = 3), item_score = rep(1:3, 5), beta = unlist(bt[1:5]))
  item <- kits_item(W, fixed_params = fixed, scale = "old")   # 0.15.0: simpleError "infinite or missing values in 'x'"
  free <- coef(item)$beta_se > 0
  expect_true(all(is.finite(coef(item)$beta)))
  expect_lt(abs(mean(coef(item)$beta[free]) - mean(unlist(bt[-(1:5)]))), 0.3)
  expect_equal(item$diagnostics$status[item$diagnostics$check == "optimizer convergence"], "ok")
})

test_that("A5 沒有有限估計的步驟（只有滿分的孩子通過）：停下來（kits_error_items），不再給 SE 上萬的 OK", {
  Y <- kits_sim(stats::rnorm(400), list(a1 = -1, a2 = 0, a3 = 0.5, a4 = 1), seed = 3)
  set.seed(2)
  Y <- cbind(Y, a5 = as.integer(rowSums(Y) == 4 & stats::runif(400) < 0.9))
  expect_error(kits_item(Y), "a5:1 \\(reached by none\\)", class = "kits_error_items")
  # a calibration that exists has an "extreme estimates" check that is ok
  dg <- item_y$diagnostics
  expect_equal(dg$status[dg$check == "extreme estimates"], "ok")
})

test_that("A5 估不出來的設計：沒有共同題、每個人只答一題，停下來（kits_error_data）並說明", {
  X1 <- kits_sim(stats::rnorm(300), list(a1 = -1, a2 = 0, a3 = 1), prefix = "y", seed = 1)
  X2 <- kits_sim(stats::rnorm(300), list(q1 = 0, q2 = 1, q3 = 2), prefix = "o", seed = 2)
  X <- matrix(NA, 600, 6, dimnames = list(c(rownames(X1), rownames(X2)), c(colnames(X1), colnames(X2))))
  X[1:300, 1:3] <- X1
  X[301:600, 4:6] <- X2
  expect_error(kits_item(X, scale = "s"), "2 separate groups that share no item: \\{a1, a2, a3\\}; \\{q1, q2, q3\\}",
               class = "kits_error_data")
  one <- X1
  for (i in seq_len(nrow(one))) one[i, -((i %% 3) + 1)] <- NA
  expect_error(kits_item(one, scale = "s"), "only one item", class = "kits_error_data")
  # with a fixed item in each group, the design is identified
  fixed <- data.frame(item_id = c("a1", "q1"), item_score = 1, beta = c(-1, 0))
  expect_s3_class(kits_item(X, fixed_params = fixed, scale = "s"), "kits_item")
})

# ---- A6 常模的界 -------------------------------------------------------------------------------------------

test_that("A6 常模表和配適用同樣的界：delta 不超過 exp(2)，P3-P97 外的比例和模型相同（約 6%）", {
  set.seed(1)
  age <- stats::runif(1500, 3, 36)
  D <- 20 + 15 * log(age) + ifelse(age > 24, stats::runif(1500, -3, 3), stats::rnorm(1500, 0, 3))
  nm <- suppressWarnings(kits_norms(D, age, ages = seq(3, 36, by = 3), degrees = c(3, 2, 1, 3)))
  expect_lte(max(nm$delta), exp(2) + 1e-12)       # 0.15.0: 12.2
  old <- age > 30
  share <- mean(abs(kits_norm_z(D[old], age[old], nm)) > 1.88, na.rm = TRUE)
  expect_lt(abs(share - 0.0618), 0.01)            # 0.15.0: 0.171
  expect_true(all(nm$sigma > 0))
})

# ---- A7 換算過的題目參數 ------------------------------------------------------------------------------------

test_that("A7 kits_subj(X, coef(kits_convert(item, d_link)))：分數在校準的量尺上，標籤也對", {
  d_link <- kits_linking(item_y, anchors)
  item_D <- kits_convert(item_y, d_link)
  s <- kits_subj(young, coef(item_D))
  expect_equal(attr(s, "scale"), "young")        # 0.15.0: "D" for logits
  expect_equal(s$theta, kits_subj(young, item_y)$theta)
  expect_equal(attr(coef(item_D), "beta_scale"), "young")
  expect_equal(attr(kits_convert(s, d_link), "scale"), "D")
})

# ---- A8 領域的點都被篩掉 -----------------------------------------------------------------------------------

test_that("A8 kits_linking()：一個領域的點都被篩掉、n_boot > 0 時不再當掉", {
  set.seed(1)
  steps <- data.frame(item_id = paste0("a", 1:12), item_score = 1, thr = c(seq(-2, 2, length.out = 10), 0, 0.5),
                      domain = c(rep("lan", 10), "cog", "cog"))
  tau <- data.frame(item_id = steps$item_id, tau = 60 + 5 * steps$thr + c(stats::rnorm(10, 0, 0.2), 9, -9))
  expect_warning(link <- kits_linking(steps, tau, n_boot = 50), "cog", class = "kits_warning_scale")
  expect_true(is.finite(link$intercept_error))
  expect_true(is.na(link$domains$intercept[link$domains$domain == "cog"]))
  expect_equal(link$diagnostics$status[link$diagnostics$check == "domains without kept anchors"], "warn")
})

# ---- A10 錯誤訊息 ------------------------------------------------------------------------------------------

test_that("A10 錨題表的欄名錯了：說錨題表要 item_id、tau（不是叫我給 thr）", {
  bad <- data.frame(item = paste0("a", 1:5), difficulty = c(40, 47, 51, 55, 60))
  expect_error(kits_linking(item_y, bad), "anchor table \\(a link to D\\) needs item_id and tau", class = "kits_error_input")
})

test_that("A10 分數從 1 開始：說先減 1，不是叫我把題目拿掉", {
  Z <- young + 1L
  e <- tryCatch(kits_item(Z, scale = "z"), error = identity)
  expect_s3_class(e, "kits_error_items")
  expect_match(conditionMessage(e), "subtract 1 so that scores start at 0")
  expect_false(grepl("leave them out", conditionMessage(e)))
  Z <- young
  Z[, "a1"] <- Z[, "a1"] + 1L
  expect_match(tryCatch(kits_item(Z, scale = "z"), error = conditionMessage), "a1.*subtract 1")
})

test_that("A10 persons 的欄名不認得：列出它有的欄，或說哪幾欄沒用到", {
  p <- data.frame(visit_id = rownames(young), wave = 12, age = 12.3, gender = "girl")
  expect_error(kits_subj(young, item_y, persons = p), "it has: visit_id, wave, age, gender",
               class = "kits_error_input")
  p$wave_m <- 12
  expect_message(kits_subj(young, item_y, persons = p), "not used: wave, age, gender", class = "kits_message_data")
})

# ---- A11 KIT 專用的預設值 ----------------------------------------------------------------------------------

test_that("A11 kits_sim()：兩次呼叫的 person id 不撞名；同樣的呼叫給同樣的 id", {
  a <- kits_sim(rnorm(5), list(i1 = 0, i2 = 1), seed = 1)
  b <- kits_sim(rnorm(5), list(i1 = 0, i2 = 1), seed = 2)
  expect_length(intersect(rownames(a), rownames(b)), 0)
  expect_identical(kits_sim(c(-1, 0, 1), list(i1 = 0), seed = 3), kits_sim(c(-1, 0, 1), list(i1 = 0), seed = 3))
  expect_match(rownames(a)[1], "^p[0-9a-f]{4}_0001$")
  expect_equal(rownames(kits_sim(0, list(i1 = 0), prefix = "y", seed = 1)), "y0001")
  expect_error(kits_sim(0, list(i1 = 0), prefix = 1), class = "kits_error_input")
})

test_that("A11 kits_dif_combine() 的世代名：預設是參數的名字（或 1、2），不是 M03、M36；ref 不在 sex 裡時列出有的值", {
  dA <- data.frame(item_id = c("a1", "a2"), dif = c(0.5, 0))
  dB <- data.frame(item_id = c("a1", "a2"), dif = c(0.4, -0.5))
  expect_equal(names(kits_dif_combine(dA, dB)), c("item_id", "dif_dA", "dif_dB", "dir_dA", "dir_dB", "DIF"))
  expect_equal(names(kits_dif_combine(dA, dplyr::mutate(dB)))[2:3], c("dif_1", "dif_2"))
  X <- rbind(young[1:250, ], young[251:500, ])
  expect_error(kits_dif(X, rep(c("F", "M"), each = 250), ref = "boy"), "groups are F, M.*ref = \"F\"",
               class = "kits_error_input")
})

test_that("A11 世代和性別的名字從資料拿：plot(kits_vertical)、照組的常模的 plot()", {
  s <- data.frame(cohort = rep(c("A", "B"), each = 400),
                  wave_m = c(rep(c(24, 36), each = 200), rep(c(36, 48), each = 200)),
                  theta = rnorm(800, 65, 4))
  expect_s3_class(plot(kits_vertical(s), type = "median"), "plotly")
  age <- runif(600, 3, 24)
  sex <- rep(c("F", "M"), 300)
  fit <- suppressWarnings(kits_norms(20 + 15 * log(age) + rnorm(600, 0, 3), age, ages = 3:24,
                                     degrees = c(3, 2, 1, 1), by = data.frame(sex = sex)))
  p <- plot(fit)
  expect_s3_class(p, "plotly")
  names_of <- vapply(plotly::plotly_build(p)$x$layout$annotations, function(a) a$text, character(1))
  expect_true(all(c("F", "M") %in% names_of))
})

# ---- B 常模的警告、存圖 --------------------------------------------------------------------------------------

test_that("B 常模的警告說出參數 degrees 和建議的值；存 html 找不到 pandoc 時有 kits 的錯誤", {
  expect_match(kits:::lower_degrees_advice(c(5, 4, 3, 2)), "`degrees`.*c\\(3, 2, 1, 1\\)")
  set.seed(2)
  age <- c(stats::runif(150, 5, 7), stats::runif(150, 11, 13))
  D <- 20 + 15 * log(age) + stats::rnorm(300, 0, 3)
  w <- tryCatch(kits_norms(D, age, ages = 6:12), warning = conditionMessage)
  if (is.character(w)) expect_match(w, "degrees")
  testthat::local_mocked_bindings(pandoc_available = function(...) FALSE, .package = "rmarkdown")
  p <- plot(link_oy)
  expect_error(kits_save_plot(p, file.path(tempdir(), "x.html")), "RSTUDIO_PANDOC", class = "kits_error_input")
})

# ---- C 專家的小問題 ----------------------------------------------------------------------------------------

test_that("C kits_subj() 和 kits_subj_table(item) 的 pattern_id 相同（合併得起來）", {
  X <- kits_sim(rnorm(300), list(a1 = c(-1, 0), a2 = c(0, 1), a3 = 0.5, a4 = -0.5), seed = 1)
  X[1:100, "a4"] <- NA
  item <- kits_item(X, scale = "s")
  s <- kits_subj(X, item)
  w <- kits_subj_table(item)
  expect_equal(sort(unique(s$pattern_id)), sort(unique(w$pattern_id)))
  expect_equal(nrow(merge(s, w, by = c("pattern_id", "raw"))), 300)   # 0.15.0: 0
  expect_equal(sort(unique(item$design$pattern_id)),
               sort(unique(item$fit$inputs$design$pattern_id)))
})

test_that("C md5 只和內容有關，和字串的編碼標記無關", {
  a <- "王01"
  b <- a
  Encoding(b) <- "unknown"
  X1 <- matrix(c(0, 1, 1, 0), 2, dimnames = list(c(a, "p2"), c("i1", "i2")))
  X2 <- X1
  rownames(X2) <- c(b, "p2")
  f1 <- kits:::data_fingerprint(kits:::resp_long(X1))
  f2 <- kits:::data_fingerprint(kits:::resp_long(X2))
  expect_identical(f1$scores_md5, f2$scores_md5)   # 0.15.0: FALSE
  # ASCII data keep the md5 of 0.15.0
  expect_identical(kits:::hash_of(list("p1", 1:3)), {
    f <- tempfile()
    writeBin(serialize(list("p1", 1:3), NULL, version = 2L)[-(1:14)], f)
    unname(tools::md5sum(f))
  })
})

test_that("C 輸入檢查：n_boot 要整數、常模的 sigma > 0、分數都一樣、factor 的分數、每波一個孩子", {
  expect_error(link_of_points(1:5, c(2, 4, 6, 8, 11), id = letters[1:5], n_boot = 2.5),
               "whole number", class = "kits_error_input")
  bad_norm <- data.frame(age_m = c(6, 12), mu = 0, sigma = -1, epsilon = 0, delta = 1)
  expect_error(kits_norm_z(1, 6, bad_norm), "positive", class = "kits_error_input")
  expect_error(kits_norms(rep(5, 400), runif(400, 3, 36), ages = 3:36, degrees = c(1, 1, 1, 1)), "same", class = "kits_error_data")
  long <- data.frame(visit_id = rep(1:3, 2), item_id = rep(c("a", "b"), each = 3),
                     item_score = factor(c(0, 1, 1, 1, 0, 1)))
  expect_error(kits_item(long), "not a factor", class = "kits_error_data")
  expect_no_warning(kits_vertical(data.frame(theta = 1:3, wave_m = 1:3, se = 0.5)))
})

test_that("C pcm_probs() 只有一個定義；kits_algorithm() 只說主要的估計", {
  expect_true(exists("pcm_probs", envir = asNamespace("kits")))
  expect_equal(names(formals(kits:::pcm_probs)), c("theta", "beta_i"))
  expect_error(kits_algorithm(kits_vertical(data.frame(theta = rnorm(20), wave_m = rep(1:2, 10)))),
               "main estimates", class = "kits_error_input")
})

test_that("C 有的連結沒有重抽時，kits_convert() 說 se_link 為什麼沒有", {
  co <- kits_linking(item_o, item_y, n_boot = 20)
  so <- kits_subj(old, item_o)
  d0 <- kits_linking(item_y, anchors)
  expect_message(out <- kits_convert(so, co, d0), "link\\(s\\) 2 .*no bootstrap draws", class = "kits_message_scale")
})

test_that("C 亂數的種類固定：使用者改了 RNGkind 也給同樣的重抽，跑完還原", {
  ref <- kits_linking(item_o, item_y, n_boot = 20)$intercept_error
  old_kind <- RNGkind()
  suppressWarnings(RNGkind(sample.kind = "Rounding"))
  on.exit(RNGkind(old_kind[1], old_kind[2], old_kind[3]))
  expect_equal(kits_linking(item_o, item_y, n_boot = 20)$intercept_error, ref)
  expect_equal(RNGkind()[3], "Rounding")
})

# ---- C 學生的小問題 ----------------------------------------------------------------------------------------

test_that("C print(item, n = 20) 印出 20 個步驟；固定的步驟有標記；nlminb code 1 的說明不像錯誤", {
  out <- capture.output(print(item_y, n = 7))
  expect_true(any(grepl("a5", out)))
  expect_false(any(grepl("more steps", out)))
  expect_true(any(grepl("more steps", capture.output(print(item_y, n = 3)))))
  fixed <- dplyr::filter(coef(item_y), item_id %in% c("a3", "a4", "a5"))
  upd <- kits_item(old, fixed_params = fixed)
  out <- capture.output(print(upd, n = Inf))
  expect_true(any(grepl("Fixed \\(fixed_params, SE 0\\): 5 step\\(s\\) of 3 item\\(s\\)", out)))
  expect_true(any(grepl("fixed", out[grepl("item_id", out)])))
  conv <- item_y$diagnostics$message[item_y$diagnostics$check == "optimizer convergence"]
  expect_match(conv, "^converged")
})

test_that("C 信度是負的：說通常是題本太容易（滿分多），要看什麼", {
  easy <- kits_sim(rnorm(400, 4), list(a1 = -2, a2 = -1.5, a3 = -1), prefix = "e", seed = 9)
  it <- kits_item(rbind(easy, kits_sim(rnorm(400), list(a1 = -2, a2 = -1.5, a3 = -1), prefix = "f", seed = 8)), scale = "e")
  s <- kits_subj(easy, it, persons = data.frame(visit_id = rownames(easy), wave_m = 48))
  row <- attr(s, "diagnostics")[grepl("reliability", attr(s, "diagnostics")$check), ]
  if (row$status != "ok") expect_match(row$message, "perfect")
  expect_match(kits:::subj_diagnostics(transform(s, se = 10))$message[4], "too easy")
})

# ---- 新發現：WLE 停在加權概似的低點 -------------------------------------------------------------------------

test_that("WLE：加權概似有兩個高點時取高的（一樣高時取 theta 小的，和 dexter 相同），不停在中間的低點", {
  # two dichotomous items 4.3 logits apart, raw score 1: the weighted likelihood has maxima at -0.569 and 0.466 (equally
  # high) and a minimum at their midpoint -0.052, where 0.15.0 stopped (its start value); dexter gives -0.569
  item <- data.frame(item_id = c("d1", "d8"), item_score = 1, beta = c(-2.197564, 2.094354))
  tab <- kits_subj_table(item)
  wl <- function(t) { p <- stats::plogis(t - item$beta); t - sum(log1p(exp(t - item$beta))) + 0.5 * log(sum(p * (1 - p))) }
  low_mode <- stats::optimize(wl, c(-3, mean(item$beta)), maximum = TRUE, tol = 1e-10)$maximum
  expect_equal(tab$theta[tab$raw == 1], low_mode, tolerance = 1e-6)   # -0.569; 0.15.0: -0.052
  th <- tab$theta[2]
  expect_gt(wl(th), wl(mean(item$beta)))
  # one maximum: unchanged
  one <- kits_subj_table(data.frame(item_id = c("a", "b", "c"), item_score = 1,
                                    beta = c(-0.5, 0, 0.5)))
  expect_equal(one$theta[2], uniroot(function(t) {
    p <- stats::plogis(t - c(-0.5, 0, 0.5)); 1 - sum(p) + sum(p * (1 - p) * (1 - 2 * p)) / (2 * sum(p * (1 - p)))
  }, c(-5, 5), tol = 1e-12)$root, tolerance = 1e-8)
})
