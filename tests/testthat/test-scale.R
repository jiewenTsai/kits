# 量尺：物件記得自己的量尺；kits_convert() 照順序套用連結並檢查
set.seed(1)

betas <- list(a1 = -2, a2 = -1, a3 = c(-0.5, 0.5), a4 = 0, a5 = c(0.5, 1.5), a6 = 2, a7 = 3)
X03 <- kits_sim(rnorm(500, 0), betas[1:5], seed = 1)
X36 <- kits_sim(rnorm(500, 2), betas[3:7], prefix = "q", seed = 2)
p03 <- data.frame(visit_id = rownames(X03), wave_m = rep(c(12, 24), each = 250),
                  age_m = rep(c(12, 24), each = 250) + runif(500, -1, 1), cohort = "M03")
p36 <- data.frame(visit_id = rownames(X36), wave_m = rep(c(36, 48), each = 250),
                  age_m = rep(c(36, 48), each = 250) + runif(500, -1, 1), cohort = "M36")
item03 <- kits_item(X03, scale = "lan M03")
item36 <- kits_item(X36, scale = "lan M36")
co <- kits_linking(item36, item03, n_boot = 0)
anchors <- data.frame(item_id = c("a1", "a2", "a3", "a4", "a5"), tau = c(40, 47, 51, 55, 60))
d_link <- kits_linking(item03, anchors, n_boot = 0)

test_that("兩個世代的完整流程：kits_convert(s36, co, d_link) 等於手寫的 a theta + b", {
  expect_equal(item03$scale, "lan M03")
  expect_equal(attr(coef(item03), "scale"), "lan M03")
  expect_equal(c(co$from_scale, co$to_scale), c("lan M36", "lan M03"))
  expect_equal(c(d_link$from_scale, d_link$to_scale), c("lan M03", "D"))
  s36 <- kits_subj(X36, item36, persons = p36)
  expect_equal(attr(s36, "scale"), "lan M36")
  D36 <- kits_convert(s36, co, d_link)
  expect_equal(D36$theta, d_link$intercept + d_link$slope * (co$intercept + co$slope * s36$theta))
  expect_equal(D36$se, d_link$slope * co$slope * s36$se)
  expect_equal(attr(D36, "scale"), "D")
  expect_equal(names(D36), names(s36))
  expect_equal(D36[c("visit_id", "flag", "wave_m")], s36[c("visit_id", "flag", "wave_m")],
               ignore_attr = TRUE)
  expect_equal(kits_info(D36)$conversions, c("lan M36 -> lan M03", "lan M03 -> D"))
  expect_output(print(D36[1:3, ]), "scale: D")
  # 一步一步換也一樣
  expect_equal(kits_convert(kits_convert(s36, co), d_link)$theta, D36$theta)
})

test_that("順序錯、連結錯時停下來（kits_error_scale），訊息說出兩個量尺", {
  s36 <- kits_subj(X36, item36, persons = p36)
  expect_error(kits_convert(s36, d_link, co), class = "kits_error_scale")
  expect_error(kits_convert(s36, d_link), "lan M36.*lan M03 -> D")
  s03 <- kits_subj(X03, item03, persons = p03)
  expect_error(kits_convert(s03, co), class = "kits_error_scale")
  expect_error(kits_convert(s36, co, co), class = "kits_error_scale")   # 同一個連結套兩次
  expect_error(kits_linking(item03, item03), class = "kits_error_scale")
  expect_error(kits_linking(kits_convert(item03, d_link), anchors), class = "kits_error_scale")
  expect_error(link_of_points(1:5, 1:5, letters[1:5], from_scale = "A", to_scale = "A"),
               class = "kits_error_scale")
  expect_error(kits_convert(s36), class = "kits_error_input")
  expect_error(kits_convert(s36, coef(item03)), class = "kits_error_input")
})

test_that("沒有記錄量尺時照樣換算，但警告（kits_warning_scale）；數值向量用 from", {
  th <- c(-1, 0, 1)
  expect_warning(out <- kits_convert(th, co), class = "kits_warning_scale")
  expect_equal(as.vector(out), co$intercept + co$slope * th)
  expect_equal(attr(out, "scale"), "lan M03")
  expect_no_warning(out <- kits_convert(th, co, d_link, from = "lan M36"))
  expect_equal(as.vector(out), d_link$intercept + d_link$slope * (co$intercept + co$slope * th))
  expect_error(kits_convert(th, d_link, from = "lan M36"), class = "kits_error_scale")
  plain <- data.frame(theta = th, se = 0.5)
  expect_warning(out <- kits_convert(plain, co), "not recorded")
  expect_equal(out$se, co$slope * plain$se)
  ms <- link_of_points(c(-2, -1, 0, 1, 2), c(-1, 0, 1, 2, 3), id = letters[1:5])
  expect_warning(kits_convert(kits_subj(X36, item36, persons = p36), ms), class = "kits_warning_scale")
  ms2 <- link_of_points(c(-2, -1, 0, 1, 2), c(-1, 0, 1, 2, 3), id = letters[1:5],
                        from_scale = "lan M36", to_scale = "lan M03")
  expect_equal(attr(kits_convert(th, ms2, from = "lan M36"), "scale"), "lan M03")
  expect_error(kits_convert(th, ms2, from = "lan M03"), class = "kits_error_scale")
})

test_that("kits_item 的換算：thr、thr_se 換，beta 不換；換過的 item 計分仍在校準的量尺", {
  conv <- kits_convert(item03, d_link)
  expect_equal(conv$scale, "D")
  expect_equal(conv$beta_scale, "lan M03")
  expect_equal(coef(conv)$thr, d_link$intercept + d_link$slope * coef(item03)$thr)
  expect_equal(coef(conv)$thr_se, d_link$slope * coef(item03)$thr_se)
  expect_equal(coef(conv)$beta, coef(item03)$beta)
  expect_output(kits:::print.kits_item(conv), "beta stays on lan M03")   # kitsPreview 載入時會蓋掉 print 方法
  s <- kits_subj(X03, conv)
  expect_equal(attr(s, "scale"), "lan M03")
  expect_equal(s$theta, kits_subj(X03, item03, persons = p03)$theta)
  tab <- kits_convert(coef(item03), d_link)
  expect_equal(attr(tab, "scale"), "D")
  expect_equal(attr(tab, "beta_scale"), "lan M03")
  expect_equal(tab$thr, coef(conv)$thr)
  # 固定參數用 beta：新的校準在 beta 的量尺上
  expect_equal(kits_item(X03, fixed_params = tab)$scale, "lan M03")
})

test_that("kits_subj_table()（有和沒有 data）也記得量尺，也可以換算", {
  st <- subj_table_of(item36, X36, wave_m = p36$wave_m)
  expect_equal(attr(st, "scale"), "lan M36")
  st_D <- kits_convert(st, co, d_link)
  expect_equal(st_D$theta, d_link$intercept + d_link$slope * (co$intercept + co$slope * st$theta))
  expect_equal(attr(st_D, "scale"), "D")
  expect_equal(attr(kits_subj_table(item36), "scale"), "lan M36")
})

test_that("預設的量尺標籤：資料的名字，或 scale = ；選列選欄留著量尺", {
  expect_equal(kits_item(X03)$scale, "X03")
  # data that are not a named object (code, do.call(), map()): the label is "scale 1", with a message; never data
  expect_message(sub <- kits_item(X03[1:300, ]), "scale 1", class = "kits_message_scale")
  expect_equal(sub$scale, "scale 1")
  expect_equal(suppressMessages(do.call(kits_item, list(X03)))$scale, "scale 1")
  expect_no_message(kits_item(X03[1:300, ], scale = "part"), class = "kits_message_scale")
  expect_equal(kits_item(X03, scale = "other")$scale, "other")
  expect_equal(kits(X03)$scale, "X03")
  expect_equal(kits(X03, scale = "lan M03")$subj |> attr("scale"), "lan M03")
  expect_error(kits_item(X03, scale = c("a", "b")), class = "kits_error_input")
  s <- kits_subj(X36, item36, persons = p36)
  expect_equal(attr(s[1:5, ], "scale"), "lan M36")
  expect_equal(attr(subset(s, theta > 0), "scale"), "lan M36")
  expect_equal(attr(head(coef(item03)), "scale"), "lan M03")
  expect_null(attr(s[1:5, ], "diagnostics"))
})

test_that("drift$fix 帶著參考的量尺：固定參數校準在同一條量尺上，並記下參考資料的指紋", {
  ref <- kits_item(kits_sim(seq(-3, 3, length.out = 400), list(a1 = c(-1, 0), a2 = c(0, 1), a3 = 0.5, a4 = -0.5, a5 = 1),
                            seed = 1), scale = "lan M03")
  new <- kits_sim(seq(-2, 4, length.out = 400), list(a1 = c(-1, 0), a2 = c(0, 1), a3 = 0.5, a4 = -0.5, a5 = 1, a6 = 0.5),
                  prefix = "n", seed = 2)
  drift <- kits_drift(new, ref)
  expect_equal(attr(drift$fix, "scale"), "lan M03")
  upd <- kits_item(new, fixed_params = drift$fix)
  expect_equal(upd$scale, "lan M03")
  expect_equal(kits_info(upd)$reference$scale, "lan M03")
  expect_equal(kits_info(upd)$reference$data, kits_info(ref)$data)
  expect_equal(attr(kits_subj(new, upd), "scale"), "lan M03")
  expect_error(kits_item(new, fixed_params = drift$fix, scale = "lan M48"), class = "kits_error_scale")
  expect_equal(kits_item(new, fixed_params = drift$fix, scale = "lan M03")$scale, "lan M03")
  # 子集合（subset）也留著量尺
  expect_equal(kits_item(new, fixed_params = subset(drift$fix, item_id != "a2"))$scale, "lan M03")
})

test_that("幾個領域共用一條 D 線：from_scale 是各領域的量尺", {
  X_cog <- kits_sim(rnorm(300), stats::setNames(betas[1:4], paste0("ca", 1:4)), seed = 3)
  item_cog <- kits_item(X_cog, scale = "cog M03")
  two <- kits_linking(list(lan = item03, cog = item_cog),
                        data.frame(item_id = c("a1", "a2", "a3", "ca1", "ca2", "ca3"), tau = c(50, 57, 60, 50, 57, 60)), n_boot = 0)
  expect_equal(unname(two$from_scale), c("lan M03", "cog M03"))
  expect_output(kits:::print.kits_link(two), "lan M03, cog M03 -> D")
  s_cog <- kits_subj(X_cog, item_cog)
  expect_equal(attr(kits_convert(s_cog, two), "scale"), "D")
  expect_equal(attr(kits_convert(kits_subj(X03, item03, persons = p03), two), "scale"), "D")
})

test_that("rbind()：量尺相同才能疊；不同時停下來（kits_error_scale）並說出兩個量尺；data.frame 混進來時警告", {
  s03 <- kits_subj(X03, item03, persons = p03)
  s36 <- kits_subj(X36, item36, persons = p36)
  err <- tryCatch(rbind(s03, s36), error = function(e) e)
  expect_s3_class(err, "kits_error_scale")
  expect_match(conditionMessage(err), "lan M03 and lan M36")
  both <- rbind(s03, kits_convert(s36, co))
  expect_s3_class(both, "kits_subj")
  expect_equal(attr(both, "scale"), "lan M03")
  expect_equal(nrow(both), nrow(s03) + nrow(s36))
  # the diagnostics of a stacked kits_subj are recomputed for all its rows
  expect_s3_class(attr(both, "diagnostics"), "kits_diagnostics")
  expect_equal(attr(both, "diagnostics")$check[1], "stacked tables")
  expect_null(attr(both, "info"))
  expect_equal(both$theta, c(s03$theta, kits_convert(s36, co)$theta))
  D <- rbind(kits_convert(s03, d_link), kits_convert(s36, co, d_link))
  expect_equal(attr(D, "scale"), "D")
  expect_warning(mixed <- rbind(s03, as.data.frame(s36)), class = "kits_warning_scale")
  expect_s3_class(mixed, "kits_subj")
  expect_equal(nrow(rbind(coef(item03), coef(item03))), 2 * nrow(coef(item03)))   # 門檻表也是 kits_table
})
