# kits 0.22.0 第三輪：作答表 kits_resp()、分數表帶著原始分換算、kits_norm_z() 回傳分數表、ages 選填、
# 題目表 kits_item_table()、R 內建的方法（confint、vcov、logLik、nobs、update、simulate），和一些印出來的說明。

long <- kitm03$responses |>
  dplyr::filter(domain == "lan", wave_m <= 12)
resp <- suppressMessages(kits_resp(long))
item <- suppressMessages(kits_item(resp))
persons <- dplyr::distinct(long, visit_id, subj_id, cohort, sex, wave_m, age_m)

# ---- kits_resp() ---------------------------------------------------------------------------------

test_that("kits_resp()：一列一個人次，孩子的欄在前、題目在後；記得題目、領域、世代", {
  expect_s3_class(resp, "kits_resp")
  expect_identical(names(resp)[1:6], c("visit_id", "subj_id", "cohort", "sex", "wave_m", "age_m"))
  expect_identical(attr(resp, "items"), sort(unique(long$item_id)))
  expect_identical(attr(resp, "domain"), "lan")
  expect_identical(attr(resp, "cohort"), "M03")
  expect_equal(nrow(resp), dplyr::n_distinct(long$visit_id))
  one <- long[long$visit_id == resp$visit_id[5], ]
  expect_equal(unlist(resp[5, one$item_id]), stats::setNames(one$item_score, one$item_id))
  expect_output(print(resp), "person-waves, .* items")
  expect_output(print(summary(resp)), "standard form")
  expect_s3_class(plot(resp), "plotly")
  kept <- dplyr::filter(resp, wave_m == 3)
  expect_s3_class(kept, "kits_resp")
  expect_identical(attr(kept, "items"), attr(resp, "items"))
})

test_that("kits_resp()：好幾個領域、人次裡會變的欄時停下來；寬的資料加 persons；沒有編號時自動編號", {
  both <- dplyr::filter(kitm03$responses, wave_m == 3)
  expect_error(kits_resp(both), "filter\\(data, domain ==", class = "kits_error_input")
  bad <- dplyr::mutate(long, age_m = age_m + seq_along(age_m) / 1e6)
  expect_error(kits_resp(bad), "age_m", class = "kits_error_data")
  X <- kits_sim(rnorm(200), list(a1 = 0, a2 = c(-1, 1), a3 = 0.5), seed = 1)
  p <- data.frame(visit_id = rownames(X), wave_m = rep(c(6, 12), each = 100), grade = "x")
  r <- kits_resp(X, persons = p)
  expect_identical(names(r), c("visit_id", "wave_m", "grade", "a1", "a2", "a3"))
  expect_error(kits_resp(X, persons = p[-1, ]), "no row in `persons`", class = "kits_error_data")
  expect_message(r2 <- kits_resp(unname(X) |> `colnames<-`(c("a1", "a2", "a3"))),
                 "numbered 1 to 200",
                 class = "kits_message_data")
  expect_identical(r2$visit_id[1:3], c("1", "2", "3"))
  expect_message(kits_item(unname(X) |> `colnames<-`(c("a1", "a2", "a3")), scale = "x"), "numbered",
                 class = "kits_message_data")
})

test_that("作答表給 kits_item()、kits_subj()、kits_item_table()、kits_dif()：和自己對齊的結果相同", {
  expect_message(kits_item(resp), "\"lan M03\" \\(the domain and cohort",
                 class = "kits_message_scale")
  expect_identical(item$scale, "lan M03")
  direct <- kits_item(long, scale = "lan M03")
  expect_identical(coef(item), coef(direct))
  s <- kits_subj(resp, item)
  s_old <- kits_subj(long, item, persons = persons)
  expect_equal(as.data.frame(s)[order(s$visit_id), ],
               as.data.frame(s_old)[order(s_old$visit_id), ],
               tolerance = 0, ignore_attr = TRUE)
  expect_true(all(c("subj_id", "wave_m", "age_m", "sex", "cohort") %in% names(s)))
  expect_false("domain" %in% names(s))
  X <- as.matrix(as.data.frame(resp)[attr(resp, "items")])
  rownames(X) <- resp$visit_id
  norm <- data.frame(domain = NA_character_, age_m = 3:12, P10 = -4, P50 = 0, P90 = 4)
  tab <- kits_item_table(coef(item), norm, data = resp)
  ap <- kits:::age_pass(X, resp$age_m, resp$wave_m)
  expect_identical(tab$age50, ap$p50[match(tab$item_id, ap$item_id)])
  expect_identical(attr(tab, "diagnostics"), attr(ap, "diagnostics"))
  quiet_data <- function(expr) suppressWarnings(expr, classes = "kits_warning_data")   # 少數題兩組的分數不同
  d <- quiet_data(kits_dif(resp, by = "sex"))
  d2 <- quiet_data(kits_dif(X, resp$sex))
  expect_equal(d, d2, tolerance = 0, ignore_attr = c("info", "lrtest"))   # lrtest$split：怎麼分的
  expect_identical(attr(d, "lrtest")$LR, attr(d2, "lrtest")$LR)
  expect_error(kits_dif(X), "kits_resp", class = "kits_error_input")
})

# ---- 分數表：subj_id、原始分換算、z -------------------------------------------------------------------------

test_that("scores |> kits_subj_table()：換尺、篩選、疊兩個世代", {
  s <- kits_subj(resp, item)
  new <- kits_subj_table(s)
  old <- new   # 0.22.0 和以前的 kits_subj_table(item, data, wave_m) 完全相同（0.23.0 拿掉了那個寫法）
  expect_output(print(new), "the standard form of each wave of the data, with counts")
  expect_output(print(kits_subj_table(item)),
                "every raw score of every item pattern, from the item parameters")
  # 換尺：分數換過，換算表也照同樣的連結換
  later <- suppressMessages(kits_item(dplyr::filter(long, wave_m >= 6), scale = "lan later"))
  d_link <- kits_linking(item, later, n_boot = 30)
  converted <- kits_subj_table(kits_convert(s, d_link))
  expect_equal(converted, kits_convert(old, d_link), tolerance = 0,
               ignore_attr = c("info", "diagnostics"))
  expect_true("se_link" %in% names(converted))
  # 篩選：人數是現在的列
  six <- kits_subj_table(dplyr::filter(s, wave_m == 6))
  expect_identical(unique(six$wave_m), 6)
  expect_equal(six$n, old$n[old$wave_m == 6])
  # 兩個世代疊起來：各自的題目組合
  long36 <- kitm36$responses |>
    dplyr::filter(domain == "lan", wave_m <= 48)
  resp36 <- suppressMessages(kits_resp(long36))
  item36 <- suppressMessages(kits_item(resp36))
  s36 <- kits_subj(resp36, item36)
  attr(s36, "scale") <- "lan M03"   # 只為了疊起來（不換尺）
  stacked <- kits_subj_table(rbind(s, s36))
  expect_equal(dplyr::filter(stacked, wave_m <= 12), new, tolerance = 0, ignore_attr = TRUE)
  expect_equal(dplyr::filter(stacked, wave_m >= 36)$theta, kits_subj_table(s36)$theta,
               tolerance = 0)
})

test_that("kits_norm_z(scores, norm) 回傳分數表（z、pct）；s$z <- kits_norm_z(s, norm) 這一版還能用，有警告", {
  s <- kits_subj(resp, item)
  norm <- suppressWarnings(kits_norms(s, degrees = c(2, 1, 1, 1)))
  expect_identical(kits_info(norm)$settings$ages_from, "kits_norm_ages()")
  expect_equal(range(norm$age_m), range(kits_norm_ages(s)))
  expect_output(print(norm), "set by kits_norm_ages\\(\\)")
  z <- s |>
    kits_norm_z(norm)
  expect_s3_class(z, "kits_subj")
  expect_identical(z$z, kits_norm_z(s$theta, s$age_m, norm))
  expect_identical(z$pct, 100 * pnorm(z$z))
  expect_identical(attr(z, "scale"), attr(s, "scale"))
  old <- s
  expect_error(old$z <- kits_norm_z(s, norm), "s <- s \\|> kits_norm_z\\(norm\\)",
               class = "kits_error_input")
})

test_that("print(norms)：分數沒有先篩、有波次會被 kits_norm_screen() 拿掉時，說一聲", {
  set.seed(4)
  X <- kits_sim(rnorm(600), list(a1 = 0, a2 = 1, a3 = 2), seed = 1)
  p <- data.frame(visit_id = rownames(X), wave_m = rep(c(6, 12, 18), each = 200))
  p$age_m <- p$wave_m + runif(600, -1, 1)
  s <- kits_subj(X, kits_item(X, scale = "D"), persons = p)
  norm <- suppressWarnings(kits_norms(s, degrees = c(1, 1, 1, 1)))
  expect_output(print(norm), "Note: wave\\(s\\) .* kits_norm_screen\\(\\)")
})

# ---- 題目表 -----------------------------------------------------------------------------------------

test_that("kits_item_table()：和內部的 item_positions() 相同；plot() 畫成長曲線和階梯", {
  set.seed(1)
  age <- runif(800, 3, 24)
  D <- 20 + 15 * log(age) + rnorm(800, 0, 3)
  norm <- kits_norms(D, age, ages = 3:24, degrees = c(3, 2, 1, 1),
                     by = data.frame(domain = rep("lan", 800)))
  thr <- data.frame(item_id = c("a1", "a2", "a3", "a4"), item_score = 1, thr = c(40, 55, 65, 60))
  pass <- pass_resp(c(a1 = 4.5, a2 = 11, a3 = Inf, a4 = 15), age)
  dif <- data.frame(item_id = c("a2", "a3", "a4"),
                    DIF = c("favors girls", "favors boys", "inconsistent"))
  stems <- c(a1 = "Lifts head", a2 = "Says two words", a3 = "Names colours",
             a4 = "Points to pictures")
  tab <- kits_item_table(thr, norm, data = pass, dif = dif, stems = stems)
  expect_s3_class(tab, "kits_item_table")
  expect_identical(is.na(tab$age50), c(FALSE, FALSE, TRUE, FALSE))   # a3：沒有人通過
  expect_identical(tab$DIF, c("none", "favors girls", "favors boys", "inconsistent"))
  items <- data.frame(item_id = thr$item_id, score = "lan", D = thr$thr, age50 = tab$age50,
                      DIF = tab$DIF,
                      stem = unname(stems))
  # 0.21.0 的表（領域的欄叫 score）給內部的函數（以前的 kits_item_positions() 等）
  old_norm <- dplyr::rename(as.data.frame(norm), score = domain)
  pos <- old_item_positions(items, old_norm)
  expect_true("score" %in% names(pos))
  new_names <- sub("^score$", "domain", names(pos))
  expect_equal(as.data.frame(tab)[new_names], as.data.frame(pos), tolerance = 0,
               ignore_attr = TRUE)
  expect_identical(tab$domain, rep("lan", 4))
  growth <- plotly::plotly_build(plot(tab))$x$data
  old_growth <- old_plot_growth(old_norm, pos, "lan")
  expect_equal(lapply(growth, function(d) d$y),
               lapply(plotly::plotly_build(old_growth)$x$data, function(d) d$y))
  ladder <- plot(tab, type = "ladder", ages = c(6, 12, 18))
  marks <- old_ladder_marks(old_norm, "lan", c(6, 12, 18))
  old_ladder <- old_plot_ladder(pos |> transform(column = "lan"), "lan", marks)
  expect_equal(lapply(plotly::plotly_build(ladder)$x$data, function(d) d$y),
               lapply(plotly::plotly_build(old_ladder)$x$data, function(d) d$y))
  expect_length(plotly::plotly_build(plot(tab, items = c("a1", "a2")))$x$data, 5)
  one_dif <- data.frame(item_id = c("a1", "a2"), dif = c(0.5, -0.1))
  expect_identical(kits_item_table(thr, norm, dif = one_dif)$DIF,
                   c("favors girls", "none", "none", "none"))
  twice <- rbind(thr, transform(thr[1, ], thr = 99))
  expect_message(t2 <- kits_item_table(twice, norm), "a1", class = "kits_message_data")
  expect_equal(t2$D, thr$thr)
})

# ---- R 內建的方法 -------------------------------------------------------------------------------------

test_that("confint()、vcov()、logLik()、nobs()、AIC()、update()、simulate()", {
  ci <- confint(item)
  expect_equal(ci$lower, ci$thr - qnorm(0.975) * ci$thr_se)
  expect_equal(confint(item, level = 0.9)$upper, ci$thr + qnorm(0.95) * ci$thr_se)
  V <- vcov(item)
  expect_equal(sqrt(diag(V)), item$fit$coef$SE_beta, ignore_attr = TRUE)
  expect_identical(rownames(V)[1],
                   paste0(item$fit$coef$item_id[1], ":", item$fit$coef$item_score[1]))
  ll <- logLik(item)
  expect_equal(as.numeric(ll), -item$fit$est$nll)
  expect_equal(attr(ll, "df"), nrow(coef(item)) - 1)
  expect_equal(nobs(item), item$n_persons)
  set.seed(2)
  age <- runif(600, 3, 36)
  D <- 20 + 15 * log(age) + rnorm(600, 0, 3)
  norm <- suppressWarnings(kits_norms(D, age, ages = 3:36, degrees = c(3, 2, 1, 1)))
  expect_equal(as.numeric(logLik(norm)), -norm$neg_loglik[1])
  expect_equal(attr(logLik(norm), "df"), 11)
  expect_equal(nobs(norm), norm$n[1])
  lower <- suppressWarnings(update(norm, degrees = c(2, 1, 1, 1)))
  expect_identical(attr(lower, "kits_args")$degrees, c(2, 1, 1, 1))
  same <- suppressWarnings(kits_norms(D, age, ages = 3:36, degrees = c(2, 1, 1, 1)))
  expect_equal(lower$P50, same$P50, tolerance = 0)
  expect_equal(AIC(norm), 2 * norm$neg_loglik[1] + 2 * 11)
  sim <- simulate(item, theta = c(-1, 0, 1), seed = 3)
  expect_identical(colnames(sim), unique(coef(item)$item_id))
  expect_length(simulate(item, nsim = 2, theta = 0, seed = 3), 2)
  theta <- seq(-3, 3, length.out = 300)
  boys <- kits_sim(theta, list(a1 = c(-1, 0), a2 = 0, a3 = 1, a4 = 0.5), prefix = "b", seed = 1)
  girls <- kits_sim(theta, list(a1 = c(-1, 0), a2 = 0, a3 = 2, a4 = 0.5), prefix = "g", seed = 2)
  dif <- kits_dif(rbind(boys, girls), rep(c("boy", "girl"), each = 300))
  ci_dif <- confint(dif)
  expect_equal(ci_dif$upper - ci_dif$lower, 2 * qnorm(0.975) * unname(attr(dif, "se")))
  link <- kits_linking(kits_item(boys, scale = "b"), kits_item(girls, scale = "g"), n_boot = 40)
  ci_link <- confint(link)
  expect_identical(ci_link$term, c("intercept", "slope"))
  expect_equal(ci_link$lower[1], unname(quantile(link$draws[, "intercept"], 0.025)))
})

# ---- 其他：Wright map、疊表、錨題的領域 -------------------------------------------------

test_that("plot(item, type = \"wright\", scores = )；bind_rows() 也合併連結的紀錄；錨題", {
  s <- kits_subj(resp, item)
  w <- plotly::plotly_build(plot(item, type = "wright", scores = s))$x$data
  expect_identical(w[[1]]$type, "bar")
  expect_error(plot(item, type = "wright"), "scores", class = "kits_error_input")
  co <- kits_linking(item, kits_item(long[long$wave_m != 3, ], scale = "x"))   # 沒有重抽的連結
  a <- suppressMessages(kits_convert(coef(item), co))
  b <- coef(item)
  attr(b, "scale") <- "x"
  attr(b, "links_without_se") <- data.frame(link = "y -> x", reason = "no draws")
  expect_setequal(attr(rbind(a, b), "links_without_se")$link, c("lan M03 -> x", "y -> x"))
  expect_setequal(attr(purrr::list_rbind(list(a, b)), "links_without_se")$link,
                  c("lan M03 -> x", "y -> x"))
  expect_warning(dplyr::bind_rows(a, b), "rbind\\(\\) or purrr::list_rbind",
                 class = "kits_warning_scale")
  anchors <- kits_anchors()
  expect_identical(sort(unique(anchors$domain)), c("cog", "lan", "motor"))
  expect_equal(sum(anchors$domain == "motor"), sum(startsWith(anchors$item_id, "body")))
})

test_that("接到 D：沒有錨題的領域用同一條線；有錨題的數字不變；印出、settings 都說", {
  lan <- kits_item(kits_sim(rnorm(500), list(l1 = -1, l2 = 0, l3 = 1, l4 = 0.5), seed = 1),
                   scale = "lan M03")
  mot <- kits_item(kits_sim(rnorm(500), list(m1 = -1, m2 = 0, m3 = 1, m4 = 0.5), seed = 2),
                   scale = "motor M03")
  soc <- kits_item(kits_sim(rnorm(500), list(s1 = -1, s2 = 0, s3 = 1), seed = 3),
                   scale = "social M03")
  anchors <- data.frame(item_id = c("l1", "l2", "l3", "m1", "m2", "m3"),
                        tau = c(40, 47, 52, 41, 48, 53))
  two <- kits_linking(list(lan = lan, motor = mot), anchors)
  three <- kits_linking(list(lan = lan, motor = mot, social = soc), anchors)
  expect_identical(c(three$intercept, three$slope), c(two$intercept, two$slope))
  expect_identical(three$points, two$points)
  expect_identical(three$settings$no_anchors, "social")
  expect_null(two$settings$no_anchors)
  expect_output(print(three), "same line, no anchors: social")
  s_soc <- kits_subj(kits_sim(rnorm(20), list(s1 = -1, s2 = 0, s3 = 1), seed = 4), soc)
  expect_equal(kits_convert(s_soc, three)$theta, three$intercept + three$slope * s_soc$theta,
               tolerance = 0)
  expect_error(kits_linking(list(social = soc), anchors), "at least 3", class = "kits_error_items")
})

# ---- 初學者教學找到的印出問題 ------------------------------------------------------------------------

test_that("診斷表 select() 之後只印選的欄；整張表照 status、check、value、message", {
  s <- kits_subj(resp, item)
  diag <- summary(s)$diagnostics
  picked <- dplyr::select(diag, status, check, value)
  out <- capture.output(print(picked))
  expect_false(any(grepl("message", out)))
  expect_match(out[1], "^status +check +value$")
  expect_match(capture.output(print(diag))[1], "^status +check +value +message$")
})

test_that("kits_subj 選過或排過欄時，照使用者的欄和順序印，不加說明的行", {
  s <- kits_subj(resp, item)
  out <- capture.output(print(dplyr::select(s, flag, visit_id, theta)))
  expect_match(out[2], "flag +visit_id +theta")
  expect_false(any(grepl("^Next:|^Also kept:|^Other columns:", out)))
  full <- capture.output(print(s))
  expect_true(any(grepl("^Next:", full)))
  expect_true(any(grepl("^Next:", capture.output(print(dplyr::filter(s, wave_m == 3))))))
})

test_that("篩過的常模印出配的月齡範圍，和印出的列數", {
  s <- kits_subj(resp, item)
  norm <- suppressMessages(suppressWarnings(kits_norms(s, degrees = c(2, 1, 1, 1))))
  ages <- range(norm$age_m)
  head_line <- sprintf("norm fitted on %s to %s months", format(ages[1]), format(ages[2]))
  expect_output(print(norm), head_line, fixed = TRUE)
  few <- dplyr::filter(norm, age_m %in% c(6, 12))
  expect_output(print(few), paste0(head_line, " (set by kits_norm_ages()); showing 2 rows"),
                fixed = TRUE)
})

test_that("量尺不是 D 的訊息只在 kits_norms() 說一次", {
  s <- kits_subj(resp, item)
  expect_no_message(kits_norm_screen(s))
  msgs <- character()
  withCallingHandlers(suppressWarnings(kits_norms(s, degrees = c(2, 1, 1, 1))),
                      message = function(m) {
                        msgs <<- c(msgs, conditionMessage(m))
                        invokeRestart("muffleMessage")
                      })
  expect_length(grep("not D", msgs), 1)
})

test_that("用語：作答組合叫 item pattern；summary(resp) 的題號只印前幾題", {
  expect_output(print(resp), "\\d+ item pattern\\(s\\)")
  out <- capture.output(print(summary(resp)))
  expect_match(out[1], "^Item patterns of each wave")
  table_lines <- out[seq_len(grep("^All item ids", out) - 1)]
  expect_true(all(nchar(table_lines) <= 100))
  expect_true(any(grepl("without ", out)))
  full <- summary(resp)$forms$items
  expect_true(max(lengths(strsplit(full, ", "))) > 4)
  expect_identical(summary(resp)$forms$n_items, lengths(strsplit(full, ", ")))
})
