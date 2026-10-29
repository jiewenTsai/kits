# kits_algorithm()：每個主題、每種結果印 40 行以內；有步驟、式子、文獻；教學程式在、段的名字在程式裡；
# 結果多一行它自己的事實。

set.seed(1)
betas <- list(a1 = -2, a2 = -1, a3 = c(-0.5, 0.5), a4 = 0, a5 = c(0.5, 1.5), a6 = 2, a7 = 3)
X03 <- kits_sim(rnorm(400, 0), betas[1:5], seed = 1)
X36 <- kits_sim(rnorm(400, 2), betas[3:7], prefix = "q", seed = 2)
item03 <- kits_item(X03, scale = "lan M03")
item36 <- kits_item(X36, scale = "lan M36")
anchors <- data.frame(item_id = c("a1", "a2", "a3", "a4", "a5"), tau = c(40, 47, 51, 55, 60))
co <- kits_linking(item36, item03, n_boot = 20)
d_link <- kits_linking(item03, anchors)
s36 <- suppressMessages(kits_convert(kits_subj(X36, item36), co, d_link))
sex <- rep(c("boy", "girl"), 200)
age <- runif(400, 2, 6)

results <- list(
  item = item03, kits = kits(X03, scale = "lan M03"), stats = item03$stats, subj = s36, link = co,
  anchors = d_link, norms = suppressWarnings(kits_norms(s36$theta, age_m = age, ages = 2:6, degrees = c(1, 1, 1, 1))),
  norm_screen = kits_norm_screen(s36$theta, rep(c(3, 6), 200)), dif = kits_dif(X03, sex),
  drift = suppressMessages(kits_drift(X03, item03)),
  precision = kits_precision(item03, list(all = paste0("a", 1:5))), information = kits_information(item03),
  expected_score = kits_expected_score(item03)
)

script_parts <- function(file) {
  lines <- readLines(system.file("scripts", file, package = "kits", mustWork = TRUE))
  sub("^# ==== part: ([a-z_]+) ====$", "\\1", grep("^# ==== part: [a-z_]+ ====$", lines, value = TRUE))
}

check_algorithm <- function(a, info) {
  out <- utils::capture.output(print(a))
  expect_s3_class(a, "kits_algorithm")
  expect_lte(length(out), 40)
  expect_gte(nrow(a$steps), 3)
  expect_true(all(nzchar(a$steps$formula)), info = info)
  refs <- out[seq(match("References:", out) + 1, match("Read more:", out) - 1)]
  expect_true(length(refs) >= 2 && length(refs) <= 4, info = info)
  expect_true(all(grepl("\\(\\d{4}\\)", refs) | grepl("^  kits:", refs)), info = info)
  expect_true(any(grepl("|>", a$components, fixed = TRUE)), info = info)
  if (a$topic != "norm_screen") {   # 篩波次是封閉式，沒有 algo_ 的函數
    expect_true(any(grepl("algo_", a$components, fixed = TRUE)), info = info)
  }
  expect_true(nzchar(system.file("scripts", a$script, package = "kits")), info = info)
  expect_true(all(setdiff(a$parts, "warm") %in% script_parts(a$script)), info = info)   # wle.R 的段 warm 另一組正在加
  expect_true(any(grepl("^  tutorial_algorithms/0[1-6]_[a-z]+\\.html$", out)), info = info)
  out
}

test_that("每個主題：40 行以內，步驟和式子、兩層的函數、2 到 4 篇文獻、教學程式和它所有的段", {
  topics <- c("cml", "wle", "link", "norms", "norm_screen", "dif", "drift", "information")
  for (t in topics) {
    a <- NULL
    out <- utils::capture.output(a <- kits_algorithm(t))
    expect_identical(a$topic, t)
    expect_null(a$facts)
    expect_identical(check_algorithm(a, t), out)
  }
  # 每份教學程式的段都列出來
  all_parts <- map(topics, \(t) {
    utils::capture.output(a <- kits_algorithm(t))
    a
  })
  for (f in unique(map_chr(all_parts, "script"))) {
    listed <- unlist(map(Filter(\(a) a$script == f, all_parts), "parts"))
    expect_setequal(setdiff(listed, "warm"), setdiff(script_parts(f), "warm"))
  }
})

test_that("每種結果：主題對、多一行事實、40 行以內；回傳 invisible", {
  want <- c(item = "cml", kits = "cml", stats = "cml", subj = "wle", link = "link", anchors = "link",
            norms = "norms", norm_screen = "norm_screen", dif = "dif", drift = "drift",
            precision = "information", information = "information", expected_score = "information")
  for (n in names(results)) {
    a <- NULL
    out <- utils::capture.output(a <- kits_algorithm(results[[n]]))
    expect_identical(a$topic, want[[n]], info = n)
    expect_length(a$facts, 1)
    expect_true(startsWith(out[2], "This result: "), info = n)
    check_algorithm(a, n)
  }
  utils::capture.output(expect_invisible(kits_algorithm("cml")))
  utils::capture.output(expect_invisible(kits_algorithm(item03)))
})

test_that("事實：題數、步數、人數、量尺；連結的量尺和規則；換算的鏈", {
  algo_line_of <- function(x) utils::capture.output(kits_algorithm(x))[2]
  expect_match(algo_line_of(item03), "5 items, 7 steps, 400 person-waves, 1 item pattern(s); scale lan M03", fixed = TRUE)
  expect_match(algo_line_of(s36), "converted lan M36 -> lan M03 -> D", fixed = TRUE)
  expect_match(algo_line_of(co), "lan M36 -> lan M03: 5 points", fixed = TRUE)
  expect_match(algo_line_of(co), "|e| > 2 SD(e)", fixed = TRUE)
  expect_match(algo_line_of(results$dif), "groups boy (ref), girl; overall test", fixed = TRUE)
})

test_that("不認得的輸入說清楚", {
  expect_error(kits_algorithm("item"), "one topic", class = "kits_error_input")
  expect_error(kits_algorithm(c("cml", "wle")), "one topic", class = "kits_error_input")
  expect_error(kits_algorithm(data.frame(a = 1)), "main estimates", class = "kits_error_input")
})

test_that("wle：Jeffreys 的眾數估計（PCM 時就是 Warm 的 WLE），algo_jme()，Magis 2015", {
  utils::capture.output(a <- kits_algorithm("wle"))
  expect_true(any(grepl("Jeffreys modal estimate", a$steps$step, fixed = TRUE)))
  expect_true(any(grepl("algo_jme(", a$components, fixed = TRUE)))
  expect_false(any(grepl("algo_wle", a$components, fixed = TRUE)))
  expect_true(any(startsWith(a$references, "Magis, D. (2015)")))
  expect_true("warm" %in% a$parts)
})

test_that("Components 裡 algo_ 的每一行照印的樣子跑得動（一個主題一個環境，照順序）", {
  sim <- list(
    X = X03, X_new = X36, sex = sex, item = item03, theta = c(-1, 0, 1), beta = betas$a3,
    betas = betas[1:5], from = item36, to = item03, link = co, D = s36$theta, age_m = 3 + 33 * runif(400),
    thr_new = co$points$from, thr_ref = co$points$to, item_id = co$points$item
  )
  for (t in names(kits:::ALGO_TOPICS)) {
    env <- list2env(sim, parent = environment())
    for (line in kits:::ALGO_TOPICS[[t]]$components) {
      if (!grepl("algo_", line, fixed = TRUE) && !startsWith(line, "pts <- ")) next   # pts：link 的 algo_ 行用到
      code <- sub("\\s+#.*$", "", line)
      expect_no_error(suppressWarnings(eval(parse(text = code), env)), message = paste(t, code))
    }
  }
})
