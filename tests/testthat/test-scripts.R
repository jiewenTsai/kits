# 教學程式（inst/scripts/，教學檔 source 它們）：在示範資料上、乾淨的環境裡跑，數字和 kits 相同。
# compare 之後的段是其他情況（固定題、換算、其他篩點規則、照組的常模），整份程式跑的時候也跑，每段和 kits 比。
# 分段的名字和 kits_algorithm() 列的段在 test-algorithm.R。容許誤差：門檻、theta、連結 1e-6；常模 1e-4

SCRIPTS <- c("cml.R", "wle.R", "link.R", "norm_screen.R", "norms.R", "dif.R", "drift.R", "information.R")

script_file <- function(f) system.file("scripts", f, package = "kits", mustWork = TRUE)

# 分段：第一個標記之前是開頭（說明和 library()），之後每個標記一段（名字是標記裡的名字）
script_sections <- function(lines) {
  mark <- "^# ==== (part: )?([a-z_]+) ====$"
  is_mark <- grepl(mark, lines)
  sections <- split(lines, factor(cumsum(is_mark), levels = 0:sum(is_mark)))
  names(sections) <- c("preamble", sub(mark, "\\2", lines[is_mark]))
  sections
}

# 學生做的事：在乾淨的環境裡 source() 整份程式
run_script <- function(f) {
  env <- new.env(parent = globalenv())
  utils::capture.output(suppressMessages(source(script_file(f), local = env)))
  env
}

# 一段裡 "# With kits" 開頭的幾行（到空行為止）：其他情況的段用 kits 準備輸入和比較
kits_chunks <- function(x) {
  chunk <- cumsum(x == "")
  heads <- vapply(split(x, chunk), \(l) c(l[l != ""], "")[1], "")
  chunk %in% as.integer(names(heads)[grepl("^# With kits", heads)])
}

# 教學的段落（inputs 和 compare 以外，"# With kits" 的幾行也不算）
teaching_lines <- function(f) {
  s <- script_sections(readLines(script_file(f)))
  s <- s[!names(s) %in% c("preamble", "inputs", "compare")]
  unlist(lapply(s, \(x) x[!kits_chunks(x)]))
}

# compare 之前的段（預設的情況）和之後的段（其他情況）
sections_of <- function(f) names(script_sections(readLines(script_file(f))))
default_parts <- function(f) {
  s <- sections_of(f)
  setdiff(s[seq_along(s) < match("compare", s)], c("preamble", "inputs"))
}
case_parts <- function(f) {
  s <- sections_of(f)
  s[seq_along(s) > match("compare", s)]
}

test_that("教學程式：只有 ASCII、每行最多 100 字、沒有分號；教學的段落不用 kits，導數不自己寫", {
  for (f in SCRIPTS) {
    lines <- readLines(script_file(f))
    expect_true(all(utf8ToInt(paste(lines, collapse = "\n")) < 128), info = f)
    expect_lte(max(nchar(lines)), 100)
    expect_false(any(grepl(";", sub("#.*$", "", lines))), info = f)
    expect_true(any(grepl("^# ==== inputs ====$", lines)), info = f)
    expect_true(any(grepl("^# ==== compare ====$", lines)), info = f)
    teach <- teaching_lines(f)
    expect_false(any(grepl("kits::|kits_[a-z]+\\(|library\\(kits\\)", sub("#.*$", "", teach))), info = f)
    expect_false(any(grepl("dexter|cNORM::|TAM::", sub("#.*$", "", teach))), info = f)
    s <- script_sections(lines)
    for (p in default_parts(f)) expect_false(any(kits_chunks(s[[p]])), info = p)
    for (p in case_parts(f)) expect_true(any(kits_chunks(s[[p]])), info = p)   # 每段和 kits 比
  }
})


test_that("cml.R：RTMB 的 CML 和門檻重現 kits_item()（示範資料 kitm03 的語言題）", {
  env <- run_script("cml.R")
  want <- coef(env$item)
  got <- env$thresholds
  expect_identical(got$item_id, want$item_id)
  expect_equal(got$item_score, want$item_score)
  for (v in c("beta", "beta_se", "thr", "thr_se")) expect_lt(max(abs(got[[v]] - want[[v]])), 1e-6)
  expect_equal(env$n_ge, env$item$stats$steps$n)                   # part "stats" = item$stats
  expect_equal(sum(unlist(env$patterns$n_raw)), sum(env$item$stats$raw$n))
  expect_true(all(env$diffs < 1e-6))
  # part "fixed": fixed betas set the origin; the fixed steps have SE 0
  got <- env$thresholds_fixed
  want <- env$want_fixed
  expect_identical(got$item_id, want$item_id)
  for (v in c("beta", "beta_se", "thr")) expect_lt(max(abs(got[[v]] - want[[v]])), 1e-6)
  expect_true(all(want$beta_se[want$item_id %in% env$fixed$item_id] == 0))
  expect_gt(abs(sum(want$beta)), 1)                                 # no longer sum 0
  # part "fixed_out": the item with one score only is left out, its thresholds copied, SE 0
  got <- env$thresholds_easy
  want <- env$want_easy
  expect_identical(env$item_easy$left_out, "lana00")
  expect_identical(got$item_id, want$item_id)
  for (v in c("beta", "beta_se", "thr")) expect_lt(max(abs(got[[v]] - want[[v]])), 1e-6)
  expect_equal(want$thr_se[want$item_id == "lana00"], c(0, 0))
})

test_that("wle.R：RTMB 的 WLE、SE = 1/√I、個人適配重現 kits_subj()", {
  env <- run_script("wle.R")
  want <- kits_subj(env$X, env$item)
  got <- env$subj[match(want$visit_id, env$subj$visit_id), ]
  expect_equal(got$raw, want$raw)
  expect_equal(got$raw_max, want$raw_max)
  for (v in c("theta", "se", "infit", "outfit")) {
    expect_lt(max(abs(got[[v]] - want[[v]]), na.rm = TRUE), 1e-6)
  }
  expect_equal(is.na(got$infit), is.na(want$infit))
  expect_equal(got$flag, want$flag)
  # part "convert": the line of the link (lan M03 -> lan M36) on theta and se; se_link, se_total
  want <- kits_convert(want, env$link)
  got <- env$subj_36[match(want$visit_id, env$subj_36$visit_id), ]
  for (v in c("theta", "se", "se_link", "se_total")) expect_lt(max(abs(got[[v]] - want[[v]])), 1e-6)
  expect_identical(kits:::scale_of(want), "lan M36")
})

test_that("link.R：點、篩點、mean/sigma、連結誤差、領域重現 kits_linking()", {
  env <- run_script("link.R")
  want <- env$link
  expect_identical(env$points$id, want$points$id)
  expect_lt(max(abs(env$line - c(want$intercept, want$slope))), 1e-6)
  expect_identical(env$points$id[env$dropped], want$dropped)
  expect_lt(abs(sd(env$draws[env$ok, "intercept"]) - want$intercept_error), 1e-6)
  expect_lt(abs(sd(env$draws[env$ok, "slope"]) - want$slope_error), 1e-6)
  expect_equal(sum(!env$ok), want$n_failed)
  theta <- c(-1, 0, 1)
  V <- want$vcov
  se <- sqrt(V[1, 1] + theta^2 * V[2, 2] + 2 * theta * V[1, 2])
  expect_lt(max(abs(env$se_link(theta) - se)), 1e-6)
  expect_lt(abs(env$slope_s - want$domains$slope[1]), 1e-6)
  expect_lt(max(abs(env$domains$intercept - want$domains$intercept)), 1e-6)
  # part "slope": between two calibrations, the common steps, the line, slope_error and z
  want <- env$link36
  expect_identical(want$settings$kind, "calibration")
  expect_identical(env$points36$id, want$points$id)
  expect_identical(env$points36$id[env$dropped36], want$dropped)
  expect_lt(max(abs(env$line36 - c(want$intercept, want$slope))), 1e-6)
  expect_lt(abs(env$slope_error36 - want$slope_error), 1e-6)
  expect_lt(abs(env$z_slope - (want$slope - 1) / want$slope_error), 1e-6)
  # part "rules": spread = "mad", "se", unit = "item"
  expect_setequal(paste(env$rules$spread, env$rules$unit), c("mad step", "se step", "sd item", "mad item"))
  for (i in seq_len(nrow(env$rules))) {
    want <- kits_linking(env$item36, env$item_lan, spread = env$rules$spread[i], unit = env$rules$unit[i])
    expect_identical(env$rules$dropped[[i]], want$points$dropped, info = i)
    expect_lt(max(abs(c(env$rules$intercept[i], env$rules$slope[i]) - c(want$intercept, want$slope))), 1e-6)
  }
  expect_true(all(env$rules$n_dropped != sum(env$dropped36)))       # each rule drops other points
})

test_that("norm_screen.R：篩波次和月齡範圍重現 kits_norm_screen()、kits_norm_ages()", {
  env <- run_script("norm_screen.R")
  want <- kits_norm_screen(env$D, env$wave_m, cohort = env$cohort)
  expect_identical(env$normed, want$normed)
  expect_identical(env$waves$floor_pct, want$waves$floor_pct)
  expect_false(all(env$waves$normed))                               # the demo leaves a wave out
  expect_identical(env$ages, kits_norm_ages(env$age_m[want$normed], env$wave_m[want$normed]))
})

test_that("norms.R：RTMB 的 shash（nlminb）重現 kits_norms()，1e-4 以內", {
  env <- run_script("norms.R")
  want <- env$want
  expect_equal(env$norm$age_m, want$age_m)
  for (v in c("mu", "sigma", "epsilon", "delta", "P3", "P50", "P97")) {
    expect_lt(max(abs(env$norm[[v]] - want[[v]])), 1e-4)
  }
  z <- kits_norm_z(env$D, env$age_m, want)
  expect_lt(max(abs(env$z - z), na.rm = TRUE), 1e-4)
  expect_equal(is.na(env$z), is.na(z))
  expect_lt(abs(env$best$objective - want$neg_loglik[1]), 1e-4)
  # part "by": one norm per sex, and each person-wave on the curves of its sex
  want <- env$want_by
  expect_identical(attr(want, "by"), "sex")
  expect_identical(env$norm_by$sex, want$sex)
  expect_equal(env$norm_by$age_m, want$age_m)
  for (v in c("mu", "sigma", "epsilon", "delta")) expect_lt(max(abs(env$norm_by[[v]] - want[[v]])), 1e-4)
  z <- kits_norm_z(env$D, env$age_m, want, by = data.frame(sex = env$sex))
  expect_lt(max(abs(env$z_by - z), na.rm = TRUE), 1e-4)
  expect_equal(is.na(env$z_by), is.na(z))
})

test_that("dif.R：每組的 CML、dif 和整體檢定重現 kits_dif()", {
  env <- run_script("dif.R")
  want <- kits_dif(env$X, env$sex, ref = env$ref)
  expect_identical(env$dif$item_id, want$item_id)
  expect_lt(max(abs(env$dif$dif - want$dif)), 1e-6)
  expect_lt(abs(env$overall[["stat"]] - attr(want, "overall")$stat), 1e-6)
  boys <- kits_cml(kits_stats(env$X[which(env$sex == env$ref), ]))
  expect_lt(max(abs(env$fit_ref$beta - boys$coef$beta)), 1e-6)
})

test_that("drift.R：線、偏移的題和建議固定的題重現 kits_drift()", {
  env <- run_script("drift.R")
  want <- env$want
  expect_lt(max(abs(env$line - c(want$intercept, want$slope))), 1e-6)
  expect_equal(env$points$judged, want$points$judged)
  expect_lt(max(abs(env$points$resid - want$points$resid)), 1e-6)
  expect_identical(env$items$drift, want$items$drift)
  expect_identical(env$items$suggest_fix, want$items$suggest_fix)
  expect_identical(env$items$reason, want$items$reason)
  expect_true(any(env$items$drift))
})

test_that("information.R：測驗訊息、期望原始分重現 kits；期望原始分的斜率（RTMB）是訊息", {
  env <- run_script("information.R")
  theta <- seq(-5, 5, by = 0.5)
  expect_lt(max(abs(env$information(theta) - kits_information(env$item)(theta))), 1e-10)
  expect_lt(max(abs(env$expected_score(theta) - kits_expected_score(env$item)(theta))), 1e-10)
  expect_lt(max(abs(env$slope - env$information(env$theta))), 1e-10)
})
