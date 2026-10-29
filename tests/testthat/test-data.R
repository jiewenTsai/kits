# 資料直接給（像 mirt、TAM）：矩陣、寬的 data.frame、長的 data.frame；同一個轉換器檢查；人的資料用參數給
set.seed(2)

betas <- list(a1 = c(-1, 0), a2 = c(0, 1), a3 = 0.5, a4 = -0.5, a5 = c(-2, -1, 1))
X <- kits_sim(rnorm(400, 0, 1.5), betas, seed = 1)
wave <- rep(c(12, 24), each = 200)
age <- wave + runif(400, -1, 1)
sex <- sample(c("boy", "girl"), 400, TRUE)
persons <- data.frame(visit_id = rownames(X), wave_m = wave, age_m = age, sex = sex, other = 1)
item <- kits_item(X, scale = "lan M03")

wide_id <- data.frame(visit_id = rownames(X), X, row.names = NULL)
wide_rn <- as.data.frame(X)
long <- kits:::resp_long(X)
long$item_score[1] <- NA   # 長格式裡 NA 的分數是沒施測
long_full <- kits:::resp_long(X)

test_that("三種形式給同樣的校準和分數：矩陣、寬的 data.frame（visit_id 欄或列名）、長的 data.frame", {
  for (d in list(wide_id, wide_rn, long_full)) {
    it <- kits_item(d, scale = "lan M03")
    expect_equal(coef(it), coef(item), ignore_attr = "source")
    expect_equal(kits_subj(d, it)$theta, kits_subj(X, item)$theta)
    expect_equal(kits_subj(d, it)$visit_id, rownames(X))
  }
  expect_equal(kits_item(wide_id)$scale, "wide_id")
  expect_equal(kits_info(kits_item(wide_id))$input$class, "wide data.frame")
  expect_equal(kits_info(kits_item(long_full))$input$class, "long data.frame")
  expect_null(kits_info(item)$input)
  # 長格式：NA 的分數不算作答
  expect_equal(nrow(kits_subj(long, item)), 400)
  expect_equal(kits_subj(long, item)$raw[1], sum(X[1, -1]))
})

test_that("寬的 data.frame：不是數字的題、人的欄位說清楚（kits_error_input）", {
  bad <- wide_id
  bad$a3 <- as.character(bad$a3)
  expect_error(kits_item(bad), "not numbers: a3", class = "kits_error_input")
  expect_error(kits_item(data.frame(wide_id, wave_m = wave)), "person column", class = "kits_error_input")
  expect_error(kits_item(1:3), "vector", class = "kits_error_input")
  expect_error(kits_item(list(a = 1)), class = "kits_error_input")
  expect_error(kits_item(data.frame(item_id = "a1", item_score = 1)), "visit_id",
               class = "kits_error_input")
})

test_that("資料的檢查：人次和題號不重複、分數是 0 以上的整數、分數中間沒有缺（kits_error_data）", {
  Y <- X
  rownames(Y)[2] <- rownames(Y)[1]
  expect_error(kits_item(Y), "1 repeated person id\\(s\\) \\(ids not shown", class = "kits_error_data")
  w <- wide_id
  w$visit_id[2] <- w$visit_id[1]
  expect_error(kits_item(w), "repeated person id.*own `prefix`", class = "kits_error_data")
  Y <- X
  colnames(Y)[2] <- colnames(Y)[1]
  expect_error(kits_item(Y), "repeated item ids", class = "kits_error_data")
  Y <- X
  Y[1, 1] <- -1L
  expect_error(kits_item(Y), class = "kits_error_data")
  expect_error(kits_item(X * 1.5), "whole numbers", class = "kits_error_data")
  Y <- X
  Y[Y[, "a1"] == 1, "a1"] <- 2L
  Y[Y[, "a5"] == 2, "a5"] <- 3L
  expect_error(kits_item(Y), "a1 \\(0, 2\\); a5 \\(0, 1, 3\\).*recode first", class = "kits_error_data")
  expect_error(kits_item(data.frame(visit_id = rownames(Y), Y)), "recode first",
               class = "kits_error_data")
  expect_error(kits_item(long_full[c(1, 1, 2), ]), "more than once", class = "kits_error_data")
  # 用已知的參數計分時，少數人沒有出現某個分數不是錯誤
  expect_equal(nrow(kits_subj(X[X[, "a5"] %in% c(0, 2), ], item)), sum(X[, "a5"] %in% c(0, 2)))
})

test_that("沒有作答的題和人次拿掉，有訊息（kits_message_data）", {
  Y <- cbind(X, empty = NA)
  Y[5, ] <- NA
  expect_message(expect_message(it <- kits_item(Y), "empty", class = "kits_message_data"), "1 person-wave", class = "kits_message_data")
  expect_equal(it$n_persons, 399)
  # person ids are not shown by default (privacy); options(kits.show_ids = TRUE) shows them
  expect_message(kits_subj(rbind(X, zz = NA), item), "1 person-wave\\(s\\) without responses are dropped \\(ids not shown",
                 class = "kits_message_data")
  msg <- tryCatch(kits_subj(rbind(X, zz = NA), item), message = conditionMessage)
  expect_false(grepl("zz", msg))
  op <- options(kits.show_ids = TRUE)
  expect_message(kits_subj(rbind(X, zz = NA), item), "zz", class = "kits_message_data")
  lg <- rbind(long_full, data.frame(visit_id = "zz", item_id = "a9", item_score = NA))
  expect_message(expect_message(kits_subj(lg, item), "a9"), "zz")
  options(op)
})

test_that("kits_subj(persons = )：人的欄位接在 flag 後面；有波次時每個波次的信度；指紋", {
  s <- kits_subj(X, item, persons = persons[sample(400), ])
  expect_equal(names(s),
               c("pattern_id", "visit_id", "raw", "raw_max", "theta", "se", "infit", "outfit",
                 "flag",
                           "wave_m", "age_m", "sex"))
  expect_equal(s$wave_m, wave)
  expect_equal(s$theta, kits_subj(X, item)$theta)
  dg <- attr(s, "diagnostics")
  expect_true(all(c("reliability, wave 12", "reliability, wave 24") %in% dg$check))
  expect_false(any(grepl("reliability", attr(kits_subj(X, item), "diagnostics")$check)))
  fp <- kits_info(s)$persons
  expect_equal(fp$columns, c("wave_m", "age_m", "sex"))
  expect_equal(fp$n_rows, 400)
  expect_equal(kits_info(kits_subj(X, item, persons = persons))$persons, fp)   # 和列的順序無關
  expect_output(print(kits_info(s)), "Persons table: 400 person-waves, columns wave_m, age_m, sex")
  expect_null(kits_info(kits_subj(X, item))$persons)
  # 多的列沒關係；少了人、id 重複、月齡不是數字、沒有可接的欄都停下來
  more <- rbind(persons,
                data.frame(visit_id = "zz", wave_m = 36, age_m = 36, sex = "boy", other = 1))
  expect_equal(kits_subj(X, item, persons = more)$wave_m, wave)
  expect_error(kits_subj(X, item, persons = persons[-1, ]), "no row in `persons`", class = "kits_error_data")
  expect_error(kits_subj(X, item, persons = rbind(persons, persons[1, ])), "repeated", class = "kits_error_data")
  expect_error(kits_subj(X, item, persons = transform(persons, wave_m = "a")), "numbers", class = "kits_error_data")
  expect_error(kits_subj(X, item, persons = persons["visit_id"]), "none of the columns",
               class = "kits_error_input")
  expect_error(kits_subj(X, item, persons = 1:3), class = "kits_error_input")
  # kits() 也收 persons
  expect_equal(kits(X, persons = persons)$subj$wave_m, wave)
})

test_that("kits_dif() 的 group：一列一個值、或用 visit_id 命名、長格式一個作答一個值", {
  dif_of <- function(...) suppressWarnings(kits_dif(...), classes = "kits_warning_data")$dif
  ref <- dif_of(X, group = sex)
  expect_equal(dif_of(X, group = stats::setNames(rev(sex), rev(rownames(X)))), ref)
  expect_equal(dif_of(wide_id, group = sex), ref)
  lw <- long_full
  lw$sex <- sex[match(lw$visit_id, rownames(X))]
  expect_equal(dif_of(lw, group = lw$sex), ref)
  lw$sex[1] <- setdiff(c("boy", "girl"), lw$sex[1])
  expect_error(dif_of(lw, group = lw$sex), "more than one value", class = "kits_error_input")
  expect_error(dif_of(X, group = sex[-1]), "one value per row", class = "kits_error_input")
})

test_that("kits_subj_table(scores)：矩陣、寬的、長的資料相同；沒有 wave_m 時所有列是一組", {
  ref <- subj_table_of(item, X, wave_m = wave)
  persons <- data.frame(visit_id = rownames(X), wave_m = wave)
  expect_equal(kits_subj_table(kits_subj(wide_id, item, persons = persons)), ref,
               ignore_attr = "info")
  expect_equal(kits_subj_table(kits_subj(long_full, item, persons = persons)), ref,
               ignore_attr = "info")
  one <- subj_table_of(item, X)
  expect_true(all(is.na(one$wave_m)))
  expect_equal(sum(one$n), sum(ref$n))
  expect_output(print(one), "Wave all|standard form")
})

test_that("kits_itemfit、kits_drift、kits() 收三種形式", {
  expect_equal(as.data.frame(kits_itemfit(wide_id, item)), as.data.frame(kits_itemfit(X, item)), ignore_attr = "info")
  expect_equal(kits_drift(long_full, item)$items, kits_drift(X, item)$items)
  expect_equal(coef(kits(wide_rn)$item)$thr, coef(item)$thr)
})

test_that("age_pass(pass, age_m, wave_m, step) 的 p50：分數 >= step 是通過；邏輯值照用", {
  expect_equal(age50_of(X, age_m = age, wave_m = wave),
               age50_of(X >= 1, age_m = age, wave_m = wave))
  expect_equal(age50_of(X, age, wave, step = 2)[["a1"]],
               age50_of(X >= 2, age_m = age, wave_m = wave)[["a1"]])
  expect_equal(age50_of(X[, "a1"], age, wave), age50_of(X, age, wave)[["a1"]])
  # 題目表：寬的 data.frame（visit_id、age_m、wave_m 和題）、長的 data.frame 和作答表給同樣的月齡
  norm <- data.frame(domain = NA_character_, age_m = 10:26, P10 = -4, P50 = 0, P90 = 4)
  wide <- data.frame(persons[c("visit_id", "age_m", "wave_m")], X, row.names = NULL)
  tab <- kits_item_table(coef(item), norm, data = wide)
  expect_identical(tab$age50, unname(age50_of(X, age, wave)[tab$item_id]))
  expect_identical(kits_item_table(coef(item), norm, data = kits_resp(wide))$age50, tab$age50)
  long_p <- dplyr::left_join(long_full, persons[c("visit_id", "age_m", "wave_m")], by = "visit_id")
  expect_identical(kits_item_table(coef(item), norm, data = long_p)$age50, tab$age50)
  expect_error(kits_item_table(coef(item), norm, data = X), "age_m", class = "kits_error_input")
})

test_that("常模函數收 kits_subj：和向量的寫法相同；不在 D 時有訊息（kits_message_scale）", {
  set.seed(5)
  n <- 1500
  age_n <- runif(n, 3, 30)
  wave_n <- round(age_n / 3) * 3
  D <- 20 + 15 * log(age_n) + rnorm(n, 0, 3)
  s <- structure(data.frame(visit_id = paste0("p", 1:n), theta = D, se = 1, flag = "ok",
                            wave_m = wave_n, age_m = age_n,
                            cohort = "M03"), class = c("kits_subj", "kits_table", "data.frame"), scale = "D")
  expect_no_message(scr <- kits_norm_screen(s))
  expect_equal(scr$normed, kits_norm_screen(D, wave_n, cohort = rep("M03", n))$normed)
  expect_equal(scr$waves, kits_norm_screen(D, wave_n, cohort = rep("M03", n))$waves)
  ages <- seq(4, 28, by = 2)
  norm_s <- suppressWarnings(kits_norms(s, ages, degrees = c(3, 2, 1, 1)))
  norm_v <- suppressWarnings(kits_norms(D, age_n, ages, degrees = c(3, 2, 1, 1)))
  expect_equal(norm_s$P50, norm_v$P50)
  expect_equal(attr(norm_s, "scale"), "D")
  expect_true(is.na(attr(norm_v, "scale")))
  expect_equal(suppressWarnings(kits_norms(s, ages = ages, degrees = c(3, 2, 1, 1)))$P50, norm_v$P50)
  expect_equal(kits_norm_z(s, norm_s)$z, kits_norm_z(D, age_n, norm_v))
  expect_equal(kits_norm_z(s, norm = norm_s)$z, kits_norm_z(D, age_n, norm_v))
  expect_error(kits_norms(s, age_n, ages), class = "kits_error_input")
  s_theta <- s
  attr(s_theta, "scale") <- "lan M03"
  expect_no_message(kits_norm_screen(s_theta))   # 篩波次不做常模：只有 kits_norms() 說一聲
  expect_message(suppressWarnings(kits_norms(s_theta, ages, degrees = c(3, 2, 1, 1))),
                 "The norms will be on scale lan M03 \\(not D\\)", class = "kits_message_scale")
  expect_error(kits_norm_z(s_theta, norm_s), class = "kits_error_scale")
  # 真的 kits_subj，沒有 age_m 欄時說清楚
  expect_error(suppressMessages(kits_norms(kits_subj(X, item), ages)), "persons")
  # persons 接上的 kits_subj 直接給常模函數
  sp <- kits_subj(X, item, persons = persons)
  expect_equal(kits_norm_ages(sp), kits_norm_ages(age, wave))
  expect_equal(suppressMessages(kits_norm_screen(sp))$normed, kits_norm_screen(sp$theta, wave)$normed)
})
