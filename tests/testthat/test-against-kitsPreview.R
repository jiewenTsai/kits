# 和 kitsPreview（dexter、cNORM）比：同樣的模擬資料，兩個套件的數字要相同
# 容許誤差：beta、thr、theta、se 1e-6；infit、outfit 1e-5；連結和篩波次完全相同；常模的 P 欄 1e-4
# 人的 SE：kits 0.25.0 起是 1 / sqrt(I)，kitsPreview（dexter）是 Warm 的 sqrt(I + (I'/2I)^2) / I，不會比較小

skip_if_not_installed("kitsPreview")
skip_if_not_installed("dexter")
skip_if_not_installed("cNORM")
set.seed(1)   # 資料每次都一樣（rnorm()、sample() 沒有自己的種子）

max_diff <- function(a, b) max(abs(a - b), na.rm = TRUE)

# 作答組合的名字：kits 0.16.0 起和 kits_subj() 相同（bk1, ...，C locale 的順序），kitsPreview 是 set1, ...；
# 所以比較每個作答組合的題目（不比名字）
pattern_sets <- function(design) {
  sets <- vapply(split(as.character(design$item_id), as.character(design$pattern_id)), \(x) paste(sort(x), collapse = " "), "")
  sort(unname(sets))
}

# 原始分換算表：每個作答組合（照它的題目）x 原始分一列
wle_by_set <- function(tab, design) {
  sets <- vapply(split(as.character(design$item_id), as.character(design$pattern_id)), \(x) paste(sort(x), collapse = " "), "")
  out <- data.frame(set = unname(sets[as.character(tab$pattern_id)]), raw = tab$raw, theta = tab$theta, se = tab$se)
  out[order(out$set, out$raw), ]
}

# kitsPreview 0.22.0 的原始分換算表：從分數來（以前是 kitsPreview::kits_subj_table(X, item, wave_m)）
preview_subj_table <- function(X, item, wave) {
  persons <- data.frame(visit_id = rownames(X), wave_m = wave)
  suppressMessages(kitsPreview::kits_subj(X, item, persons = persons)) |>
    kitsPreview::kits_subj_table()
}

same_item <- function(X, ..., tol = 1e-6, scale = NULL) {
  ik <- suppressMessages(kitsPreview::kits_item(X, ...))
  ir <- suppressMessages(kits_item(X, ..., scale = scale))   # 兩個世代要連結時，給不同的量尺標籤
  ck <- coef(ik)
  cr <- coef(ir)
  expect_equal(ck$item_id, cr$item_id)
  expect_equal(ck$item_score, cr$item_score)
  expect_lt(max_diff(ck$beta, cr$beta), tol)
  expect_lt(max_diff(ck$beta_se, cr$beta_se), tol)
  expect_lt(max_diff(ck$thr, cr$thr), tol)
  expect_lt(max_diff(ck$thr_se, cr$thr_se), tol)
  expect_equal(ik$left_out, ir$left_out)
  expect_equal(pattern_sets(ik$design), pattern_sets(ir$design))
  expect_equal(sort(ik$design$item_id), sort(ir$design$item_id))
  expect_equal(ik$n_persons, ir$n_persons)
  # ik$fit 是 dexter 的結果，欄名照 dexter（booklet_id）；kits 叫 pattern_id
  expect_equal(length(unique(ik$fit$inputs$design$booklet_id)), length(unique(ir$fit$inputs$design$pattern_id)))
  list(k = ik, r = ir)
}

same_subj <- function(X, item, tol = 1e-6, tol_fit = 1e-5) {
  sk <- suppressMessages(kitsPreview::kits_subj(X, item$k))
  sr <- kits_subj(X, item$r)
  expect_equal(names(sk), names(sr))
  expect_equal(sk$visit_id, sr$visit_id)
  expect_equal(sk$pattern_id, sr$pattern_id)
  expect_equal(sk$raw, sr$raw)
  expect_equal(sk$raw_max, sr$raw_max)
  expect_lt(max_diff(sk$theta, sr$theta), tol)
  expect_true(all(sr$se <= sk$se + tol))
  expect_lt(max_diff(sk$infit, sr$infit), tol_fit)
  expect_lt(max_diff(sk$outfit, sr$outfit), tol_fit)
  expect_equal(is.na(sk$infit), is.na(sr$infit))
  expect_equal(sk$flag, sr$flag)
  list(k = sk, r = sr)
}

# kitsPreview 的 kits_itemfit() 是 WLE 的題目適配：和 kits 的 type = "wle" 比
preview_itemfit <- function(X, item) suppressMessages(kitsPreview::kits_itemfit(X, item))

# 兩個套件的題目適配表：欄的數字相同；kitsPreview 有 flag 欄時，flag 完全相同
same_itemfit <- function(fk, fr, tol_fit = 1e-5) {
  expect_equal(fk$item_id, fr$item_id)
  expect_equal(fk$n, fr$n)
  expect_lt(max_diff(fk$pvalue, fr$pvalue), 1e-10)
  expect_lt(max_diff(fk$rir, fr$rir), 1e-10)
  expect_lt(max_diff(fk$infit, fr$infit), tol_fit)
  expect_lt(max_diff(fk$outfit, fr$outfit), tol_fit)
  if ("flag" %in% names(fk)) expect_identical(fk$flag, fr$flag)
}

types <- list(
  dichotomous = stats::setNames(as.list(seq(-2, 2, length.out = 8)), paste0("d", 1:8)),
  polytomous = list(a1 = c(-1, 0), a2 = c(0, 1), a3 = c(-2, -1, 1), a4 = c(-0.5, 0.5), a5 = c(0.5, 1.5, 2.5)),
  mixed = c(stats::setNames(as.list(c(-1.5, -0.5, 0.5, 1.5)), paste0("d", 1:4)),
            list(p1 = c(-1, 0), p2 = c(0, 1), p3 = c(-1, 0, 1), p4 = c(0.5, 1.5)))
)

for (type in names(types)) {
  b <- types[[type]]

  test_that(paste0(type, "：完整題本的校準、計分、題目適配、換算表和 kitsPreview 相同"), {
    X <- kits_sim(rnorm(600, 0, 1.5), b, seed = 1)
    item <- same_item(X)
    subj <- same_subj(X, item)
    same_itemfit(preview_itemfit(X, item$k), kits_itemfit(X, item$r, type = "wle"))
    wave <- rep(c(6, 12), length.out = 600)
    tk <- preview_subj_table(X, item$k, wave)
    tr <- subj_table_of(item$r, X, wave_m = wave)
    expect_equal(tk$n, tr$n)
    expect_lt(max_diff(tk$theta, tr$theta), 1e-6)
    expect_true(all(tr$se <= tk$se + 1e-6))
    expect_equal(attr(tk, "forms"), attr(tr, "forms")[c("wave_m", "form", "n", "items")])   # kits 多了 type 欄
  })

  test_that(paste0(type, "：kits_item() 的 item$itemfit 和 kitsPreview 的 item$itemfit 相同"), {
    X <- kits_sim(seq(-3, 3, length.out = 600), b, seed = 3)
    kitsPreview_item <- suppressMessages(kitsPreview::kits_item(X))
    skip_if(is.null(kitsPreview_item$itemfit), "kitsPreview's kits_item() has no itemfit yet")
    kits_item_r <- kits_item(X)
    fk <- kitsPreview_item$itemfit
    # kits 的 item$itemfit 是條件的版本；kitsPreview 的是 WLE 的：和 kits 的 type = "wle" 比
    fr <- kits_itemfit(X, kits_item_r, type = "wle")
    expect_equal(fk$item_id, fr$item_id)
    expect_equal(fk$n, fr$n)
    expect_lt(max_diff(fk$pvalue, fr$pvalue), 1e-10)   # 容許誤差和上面的 kits_itemfit() 相同
    expect_lt(max_diff(fk$rir, fr$rir), 1e-10)
    expect_lt(max_diff(fk$infit, fr$infit), 1e-5)
    expect_lt(max_diff(fk$outfit, fr$outfit), 1e-5)
    # flag（同一個規則）：kitsPreview 還沒有這一欄時只跳過這個比較
    skip_if(!"flag" %in% names(fk), "kitsPreview's item$itemfit has no flag yet")
    expect_identical(fk$flag, fr$flag)
  })

  test_that(paste0(type, "：很多題本（不完整設計）的校準、計分、換算表和 kitsPreview 相同"), {
    X <- kits_sim(rnorm(800, 0, 1.5), b, seed = 2)
    X[sample(length(X), round(length(X) * 0.35))] <- NA
    X <- X[rowSums(!is.na(X)) > 1, ]
    item <- same_item(X)
    expect_gt(length(unique(item$r$fit$inputs$design$pattern_id)), 20)
    subj <- same_subj(X, item)
    expect_true(any(subj$r$raw == 0))
    wk <- wle_by_set(kitsPreview::kits_subj_table(item$k), item$k$design)
    wr <- wle_by_set(kits_subj_table(item$r), item$r$design)
    expect_equal(wk$set, wr$set)
    expect_equal(wk$raw, wr$raw)
    expect_lt(max_diff(wk$theta, wr$theta), 1e-6)
    expect_true(all(wr$se <= wk$se + 1e-6))
    fk <- preview_itemfit(X, item$k)
    fr <- kits_itemfit(X, item$r, type = "wle")
    expect_lt(max_diff(fk$pvalue, fr$pvalue), 1e-10)
    expect_lt(max_diff(fk$outfit, fr$outfit), 1e-5)
    expect_lt(max_diff(fk$rir, fr$rir), 0.05)   # 很小的題本（其他題總分沒有變異）dexter 的 rir 是浮點雜訊
    # 把人數少於 5 的題本拿掉後，rir 完全相同
    book <- apply(!is.na(X), 1, function(a) paste(which(a), collapse = " "))
    big <- X[book %in% names(table(book))[table(book) >= 5], ]
    expect_lt(max_diff(preview_itemfit(big, item$k)$rir, kits_itemfit(big, item$r, type = "wle")$rir), 1e-10)
    wave <- rep(c(6, 12, 24), length.out = nrow(X))
    tk <- preview_subj_table(X, item$k, wave)   # kits 只有標準題本
    tr <- subj_table_of(item$r, X, wave_m = wave)
    expect_equal(tk$n, tr$n)
    expect_lt(max_diff(tk$theta, tr$theta), 1e-6)
  })
}

test_that("固定參數（含沒放進校準的題）和 kitsPreview 相同", {
  betas <- list(a1 = c(-1, 0), a2 = c(0, 1), a3 = 0.5, a4 = -0.5, a5 = c(-2, -1, 1), a6 = 1.5)
  X <- kits_sim(rnorm(400, 0, 1.5), betas, seed = 1)
  X[sample(length(X), 300)] <- NA
  X <- X[rowSums(!is.na(X)) > 0, ]
  item <- same_item(X)
  Y <- X
  Y[, "a4"] <- pmax(Y[, "a4"], 1L)
  fixed <- subset(coef(item$k), item_id %in% c("a1", "a4"))
  upd <- same_item(Y, fixed_params = fixed)
  expect_equal(upd$r$left_out, "a4")
  same_subj(Y, upd)
  wk <- wle_by_set(kitsPreview::kits_subj_table(upd$k), upd$k$design)
  wr <- wle_by_set(kits_subj_table(upd$r), upd$r$design)
  expect_equal(wk$set, wr$set)
  expect_lt(max_diff(wk$theta, wr$theta), 1e-6)
  all_fixed <- same_item(Y, fixed_params = coef(item$k))
  expect_equal(coef(all_fixed$r)$beta_se, rep(0, nrow(coef(all_fixed$r))))
})

test_that("kits_drift 和 drift$fix 的固定參數校準和 kitsPreview 相同", {
  betas <- list(a1 = c(-1, 0), a2 = c(0, 1), a3 = 0.5, a4 = -0.5, a5 = c(-0.5, 0.5), a6 = 1.5, a7 = -3.5)
  old <- kits_sim(seq(-3, 3, length.out = 400), betas, seed = 1)
  new_betas <- c(modifyList(betas, list(a3 = 2)), list(a8 = c(1, 2)))
  new <- kits_sim(seq(-2, 4, length.out = 400), new_betas, prefix = "n", seed = 2)
  ref <- same_item(old)
  dk <- kitsPreview::kits_drift(new, ref$k)
  dr <- kits_drift(new, ref$r)
  expect_equal(dk$items$item_id, dr$items$item_id)
  expect_equal(dk$items$drift, dr$items$drift)
  expect_equal(dk$items$reason, dr$items$reason)
  expect_equal(dk$items$suggest_fix, dr$items$suggest_fix)
  expect_lt(max_diff(dk$items$max_abs_resid, dr$items$max_abs_resid), 1e-6)
  expect_lt(abs(dk$intercept - dr$intercept), 1e-8)
  expect_lt(abs(dk$slope - dr$slope), 1e-8)
  expect_equal(dk$points$dropped, dr$points$dropped)
  expect_lt(max_diff(dk$fix$beta, dr$fix$beta), 1e-6)
  expect_lt(max_diff(dk$steps_free$thr, dr$steps_free$thr), 1e-6)
  upd <- same_item(new, fixed_params = dr$fix)
  expect_true(all(c("a3", "a8") %in% coef(upd$r)$item_id))
})

test_that("世代連結、D 連結（含 n_boot 的連結誤差、共用斜率）和 kitsPreview 完全相同", {
  betas <- list(a1 = -2, a2 = -1, a3 = c(-0.5, 0.5), a4 = 0, a5 = c(0.5, 1.5), a6 = 2, a7 = 3)
  young <- kits_sim(rnorm(500, 0), betas[1:5], seed = 1)
  old <- kits_sim(rnorm(500, 2), betas[3:7], prefix = "q", seed = 2)
  iy <- same_item(young, scale = "young")
  io <- same_item(old, scale = "old")
  lk <- kitsPreview::kits_linking(io$k, iy$k, n_boot = 200)
  lr <- kits_linking(io$r, iy$r, n_boot = 200)
  expect_lt(abs(lk$intercept - lr$intercept), 1e-8)
  expect_lt(abs(lk$slope - lr$slope), 1e-8)
  expect_lt(abs(lk$intercept_error - lr$intercept_error), 1e-8)
  expect_lt(abs(lk$slope_error - lr$slope_error), 1e-8)
  expect_equal(lk$n_failed, lr$n_failed)
  expect_equal(lk$dropped, lr$dropped)
  expect_equal(lk$n_items, lr$n_items)
  expect_equal(lk$points$id, lr$points$id)
  # 同樣的點：完全相同
  from <- c(-2, -1, 0, 1, 2, 3, 0.5, 1.5)
  to <- 1 + from
  to[8] <- to[8] + 3
  points <- data.frame(id = paste0("p", 1:8), item = rep(c("i1", "i2", "i3", "i4"), 2), from = from, to = to)
  mk <- points |>
    kitsPreview::kits_link_screen() |>
    kitsPreview::kits_link_line() |>
    kitsPreview::kits_link_error(n_boot = 100)
  mr <- link_of_points(from, to, id = paste0("p", 1:8), item = rep(c("i1", "i2", "i3", "i4"), 2),
                       n_boot = 100)
  same <- c("intercept", "slope", "intercept_error", "slope_error", "n_failed", "intercept0", "slope0", "n",
            "dropped")
  expect_identical(unclass(mk)[same], unclass(mr)[same])
  set.seed(1)
  steps <- data.frame(item_id = paste0("a", 1:12), item_score = 1, thr = rnorm(12), domain = rep(c("lan", "cog"), each = 6))
  anchors <- data.frame(item_id = steps$item_id, tau = 64 + 7 * steps$thr + rnorm(12, 0, 0.5))
  dk <- kitsPreview::kits_linking_d(steps, anchors, n_boot = 200)
  dr <- kits_linking(steps, anchors, n_boot = 200)
  expect_identical(c(dk$intercept, dk$slope, dk$intercept_error, dk$slope_error),
                   c(dr$intercept, dr$slope, dr$intercept_error, dr$slope_error))
  expect_identical(as.data.frame(dk$domains), as.data.frame(dr$domains))
})

test_that("篩波次、月齡範圍和 kitsPreview 完全相同；常模（預設次數、小次數）的 P 欄在 1e-4 內", {
  set.seed(5)
  age <- runif(2000, 3, 36)
  wave <- round(age / 6) * 6
  D <- 20 + 15 * log(age) + rnorm(2000, 0, 3)
  D[wave == 6][1:60] <- min(D[wave == 6])
  cohort <- ifelse(age > 15, "M36", "M03")
  sk <- kitsPreview::kits_norm_screen(D, wave, cohort = cohort)
  sr <- kits_norm_screen(D, wave, cohort = cohort)
  expect_identical(sk$normed, sr$normed)
  expect_identical(as.data.frame(sk$waves), as.data.frame(sr$waves))
  expect_identical(kitsPreview::kits_norm_ages(age[sk$normed], wave[sk$normed]), kits_norm_ages(age[sr$normed], wave[sr$normed]))

  ages <- seq(4, 35, by = 0.5)
  for (deg in list(c(5, 4, 3, 2), c(3, 2, 1, 1))) {
    nk <- suppressWarnings(kitsPreview::kits_norms(D, age, ages = ages, degrees = deg, pad = 0))
    nr <- suppressWarnings(kits_norms(D, age, ages = ages, degrees = deg, pad = 0))
    expect_equal(nk$age_m, nr$age_m)
    expect_equal(nk$n, nr$n)
    expect_equal(names(nk), names(nr))
    for (p in c("P3", "P10", "P50", "P90", "P97")) expect_lt(max_diff(nk[[p]], nr[[p]]), 1e-4)
    expect_lt(abs(nk$neg_loglik[1] - nr$neg_loglik[1]), 1e-4)
    expect_equal(attr(nk, "kits_args"), attr(nr, "kits_args"))
    expect_lt(max_diff(kitsPreview::kits_norm_z(D, age, nk), kits_norm_z(D, age, nr)), 1e-4)
  }
})

test_that("cNORM 的 optim 停在 code 52 的常模：kits（nlminb）不警告，和 kitsPreview（cNORM）的數字相同", {
  set.seed(1)
  age <- runif(2000, 3, 36)
  D <- 20 + 15 * log(age) + rnorm(2000, 0, 3)
  ages <- seq(4, 35, by = 0.5)
  for (deg in list(c(3, 2, 1, 1), c(3, 2, 2, 1))) {
    nk <- suppressWarnings(kitsPreview::kits_norms(D, age, ages = ages, degrees = deg, pad = 0))
    expect_no_warning(nr <- kits_norms(D, age, ages = ages, degrees = deg, pad = 0))
    for (p in c("P3", "P50", "P97")) expect_lt(max_diff(nk[[p]], nr[[p]]), 1e-4)
    expect_lt(abs(nk$neg_loglik[1] - nr$neg_loglik[1]), 1e-4)
  }
})

test_that("kits_sim 和 kits 產生同樣的資料；kits() 的結果相同", {
  betas <- list(a1 = c(-1, 0), a2 = c(0, 1), a3 = 0.5)
  expect_identical(kits_sim(c(-1, 0, 1), betas, prefix = "p", seed = 3), kitsPreview::kits_sim(c(-1, 0, 1), betas, seed = 3))
  expect_identical(unname(kits_sim(c(-1, 0, 1), betas, seed = 3)), unname(kitsPreview::kits_sim(c(-1, 0, 1), betas, seed = 3)))
  X <- kits_sim(rnorm(200), betas, seed = 4)
  kk <- kitsPreview::kits(X)
  kr <- kits(X)
  expect_lt(max_diff(coef(kk$item)$thr, coef(kr$item)$thr), 1e-6)
  expect_lt(max_diff(kk$subj$theta, kr$subj$theta), 1e-6)
})

test_that("性別 DIF：kits_dif 和 kits_dif_combine 和 kitsPreview（dexter 的 DIF()）相同", {
  # 內建 DIF：a2 對女孩容易、a5 對女孩難；兩個世代，不完整設計，二點和多點計分
  b <- list(a1 = c(-1, 0), a2 = c(0, 1), a3 = 0.5, a4 = -0.5, a5 = 0, a6 = c(-0.5, 0.5, 1.5), a7 = 1)
  girl_b <- modifyList(b, list(a2 = c(-0.6, 0.4), a5 = 0.8))
  cohort <- function(n, shift, seed) {
    set.seed(seed)
    X <- rbind(kits_sim(rnorm(n, shift), b, prefix = paste0("b", seed, "_")),
               kits_sim(rnorm(n, shift), girl_b, prefix = paste0("g", seed, "_")))
    X[sample(length(X), round(0.15 * length(X)))] <- NA
    X[rowSums(!is.na(X)) > 1, ]
  }
  X1 <- cohort(500, 0, 1)
  X2 <- cohort(500, 0.5, 2)[, -1]                        # the second cohort has no a1
  sex_of <- function(X) ifelse(startsWith(rownames(X), "b"), "boy", "girl")
  for (ref in c("boy", "girl")) {
    k1 <- kitsPreview::kits_dif(X1, sex_of(X1), ref = ref)
    r1 <- kits_dif(X1, sex_of(X1), ref = ref)
    expect_equal(r1$item_id, k1$item_id)
    expect_lt(max_diff(r1$dif, k1$dif), 1e-6)
  }
  k1 <- kitsPreview::kits_dif(X1, sex_of(X1))
  k2 <- kitsPreview::kits_dif(X2, sex_of(X2))
  r1 <- kits_dif(X1, sex_of(X1))
  r2 <- kits_dif(X2, sex_of(X2))
  expect_lt(max_diff(r2$dif, k2$dif), 1e-6)
  expect_lt(r1$dif[r1$item_id == "a5"], -0.3)            # favors boys
  expect_gt(r1$dif[r1$item_id == "a2"], 0.3)             # favors girls
  # the overall test of dexter's DIF(): the same when no item is left out (here every item is in both groups with the
  # same scores); with a left-out item kits centres the compared steps and dexter does not (see test-dif.R)
  expect_length(attr(r1, "left_out"), 0)
  resp <- kits:::resp_long(X1)
  resp$sex <- sex_of(X1)[match(resp$visit_id, rownames(X1))]
  names(resp)[names(resp) == "visit_id"] <- "person_id"   # dexter 的名字
  dx <- dexter::DIF(resp, person_property = "sex")
  expect_lt(abs(attr(r1, "overall")$stat - dx$DIF_overall$stat), 1e-5)
  expect_equal(attr(r1, "overall")$df, dx$DIF_overall$df)
  expect_lt(abs(attr(r1, "overall")$p - dx$DIF_overall$p), 1e-8)
  # combined over two cohorts: the same table
  for (cut in c(0.1, 0.3)) {
    ck <- kitsPreview::kits_dif_combine(k1, k2, cut = cut)
    cr <- kits_dif_combine(r1, r2, cut = cut, cohorts = c("M03", "M36"))   # kitsPreview's default names
    expect_equal(names(cr), names(ck))
    expect_equal(cr$item_id, ck$item_id)
    expect_equal(cr[c("dir_M03", "dir_M36", "DIF")], ck[c("dir_M03", "dir_M36", "DIF")], ignore_attr = TRUE)
    expect_lt(max_diff(cr$dif_M03, ck$dif_M03), 1e-6)
    expect_lt(max_diff(cr$dif_M36, ck$dif_M36), 1e-6)
  }
  expect_true(all(c("none", "favors girls", "favors boys") %in% kits_dif_combine(r1, r2)$DIF))
  # an item with other scores in one group: left out of the comparison in both packages
  Y <- X1
  Y[sex_of(Y) == "boy" & Y[, "a6"] %in% 3, "a6"] <- 2L
  expect_warning(ry <- kits_dif(Y, sex_of(Y)), "a6", class = "kits_warning_data")
  utils::capture.output(ky <- suppressWarnings(suppressMessages(kitsPreview::kits_dif(Y, sex_of(Y)))))
  expect_equal(ry$item_id, ky$item_id)
  expect_false("a6" %in% ry$item_id)
  expect_lt(max_diff(ry$dif, ky$dif), 1e-6)   # the size measure is the same also with a left-out item
})
