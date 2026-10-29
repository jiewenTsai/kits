# 圖（plotly）：都能建出來、資料不是空的；輸入檢查
set.seed(1)

test_that("圖的英文字：代碼有名稱；沒有對照的代碼照印", {
  expect_equal(names(LAB$domain), c("lan", "motor", "cog", "social", "emotion"))
  expect_equal(lab_name("domain", c("lan", "reading")), c("Language", "reading"))
})

test_that("short_stem 去掉括號並截短", {
  expect_equal(short_stem("Stacks two blocks (for example, wooden)"), "Stacks two blocks")
  long <- short_stem("Child can say the names of at least three family members")
  expect_lte(nchar(long), 37)
  expect_true(endsWith(long, "…"))
  expect_equal(nchar(short_stem("會說兩個字（例如媽媽）")), 5)   # 括號去掉，留 5 個字
})

test_that("比常模最大月齡還難、又沒有 age50 的題：圖會列出題號；輸入有問題時說清楚", {
  norm <- data.frame(score = "lan", age_m = c(6, 12, 18, 24), P10 = 0, P50 = c(40, 52, 60, 66), P90 = 80)
  items <- data.frame(item_id = c("a1", "hard"), score = "lan", D = c(45, 80), age50 = c(8.5, NA),
                      DIF = "none", stem = c("x", "y"))
  out <- old_item_positions(items, norm)
  expect_true(is.na(out$x[2]))
  expect_message(p <- old_plot_growth(norm, out, "lan"), "hard")
  notes <- vapply(plotly::plotly_build(p)$x$layout$annotations, function(a) a$text, character(1))
  expect_true(any(grepl("hard", notes)))
  expect_error(old_plot_growth(norm, transform(out, DIF = "girls"), "lan"), "DIF")
  expect_error(old_plot_growth(norm, out, "language"), "language")
})

test_that("kits_plot_ladder：沒有名稱的欄代碼照印", {
  items <- data.frame(item_id = "a1", score = "motor", column = "motor", D = 40, DIF = "none", stem = "x",
                      age50 = 5, age_reach = 5)
  p <- old_plot_ladder(items, columns = "motor", marks = data.frame(age_m = 6, D = 45))
  expect_equal(as.character(plotly::plotly_build(p)$x$layout$xaxis$ticktext), "motor")
})

test_that("kits_save_plot 的檔名要以 .html 結尾；p 要是 plotly", {
  expect_error(kits_save_plot(plotly::plot_ly(), "fig.htm"), ".html")
  expect_error(kits_save_plot(1, "fig.html"), "plotly")
})

test_that("四種報告的圖都能建出 plotly 物件", {
  norm <- data.frame(score = "lan", age_m = 3:24)
  norm$P50 <- 20 + 15 * log(norm$age_m)
  for (q in c(3, 10, 90, 97)) norm[[paste0("P", q)]] <- norm$P50 + qnorm(q / 100) * 4
  items <- data.frame(item_id = c("a1", "a2"), score = "lan", column = "lan", D = c(40, 60),
                      age50 = c(4.5, NA), DIF = c(NA, "favors boys"), stem = c("Lifts head", "Names colours"))
  items <- old_item_positions(items, norm)
  s <- data.frame(cohort = rep(c("M03", "M36"), each = 300), wave_m = 24, theta = rnorm(600, 63, 4))
  age <- runif(600, 3, 24)
  D <- 20 + 15 * log(age) + rnorm(600, 0, 3)
  fit <- suppressWarnings(kits_norms(D, age, ages = 3:24, degrees = c(3, 2, 1, 1)))
  plots <- list(
    old_plot_growth(norm, items, "lan"),
    old_plot_ladder(items, "lan", old_ladder_marks(norm, "lan", ages = c(6, 12))),
    plot(kits_vertical(s), type = "median"),
    plot(fit)
  )
  for (p in plots) {
    expect_s3_class(p, "plotly")
    expect_gt(length(plotly::plotly_build(p)$x$data), 0)
  }
})

test_that("圖 C：實際月齡和波次名稱差超過 1 個月時，標題下註明", {
  s <- data.frame(cohort = "M36", wave_m = rep(c(72, 84), each = 300),
                  age_m = rep(c(72.0, 86.5), each = 300),
                  theta = rnorm(600, 95, 4))
  title <- plotly::plotly_build(plot(kits_vertical(s), type = "median"))$x$layout$title
  if (is.list(title)) title <- title$text
  expect_match(title, "84-month wave 86.5")
  expect_no_match(title, "72-month wave")
})

test_that("plot() 的方法：Wright map、連結、常模都是 plotly", {
  betas <- list(a1 = c(-1, 0), a2 = c(0, 1), a3 = 0.5, a4 = -0.5)
  X <- kits_sim(seq(-3, 3, length.out = 300), betas, seed = 1)
  k <- kits(X)
  p <- plot(k)
  expect_s3_class(p, "plotly")
  built <- plotly::plotly_build(p)
  expect_gt(length(built$x$data), 1)
  from <- c(-2, -1, 0, 1, 2, 3, 0.5, 1.5)
  to <- 1 + from
  to[8] <- to[8] + 3   # 第 8 點會被篩掉
  link <- link_of_points(from, to, id = paste0("p", 1:8))
  expect_equal(link$dropped, "p8")
  pl <- kits:::plot.kits_link(link)   # 測試時 kitsPreview 也載入，直接叫這裡的方法
  expect_s3_class(pl, "plotly")
  names_of <- vapply(plotly::plotly_build(pl)$x$data, function(d) if (is.null(d$name)) "" else as.character(d$name), character(1))
  expect_true(any(grepl("Before screening", names_of)))
  expect_true(any(grepl("Screening band", names_of)))
  expect_true(all(c("Kept", "Dropped") %in% names_of))
  norm <- data.frame(age_m = c(6, 12, 18), mu = c(40, 55, 62), sigma = 5, epsilon = 0, delta = 1,
                     P3 = c(31, 46, 53), P50 = c(40, 55, 62), P97 = c(49, 64, 71))
  class(norm) <- c("kits_norms", "kits_table", "data.frame")
  expect_s3_class(plot(norm, pct = c(3, 50, 97)), "plotly")
  expect_s3_class(plot(norm, pct = c(3, 50, 97),
                       scores = data.frame(theta = c(40, 50), age_m = c(6, 12))), "plotly")
  expect_error(plot(norm, pct = 1), "P1")
  expect_error(plot(norm, scores = data.frame(theta = 1)), "age_m")
})

test_that("Wright map：人的長條是灰色（不透明），題目步驟是藍色", {
  betas <- list(a1 = c(-1, 0), a2 = c(0, 1), a3 = 0.5, a4 = -0.5)
  X <- kits_sim(seq(-3, 3, length.out = 300), betas, seed = 1)
  built <- plotly::plotly_build(kits:::plot.kits(kits(X)))   # 測試時 kitsPreview 也載入，直接叫這裡的方法
  bars <- Filter(function(d) identical(d$type, "bar"), built$x$data)
  steps <- Filter(function(d) identical(d$mode, "markers"), built$x$data)
  expect_equal(bars[[1]]$marker$color, PAL$gray)
  expect_null(bars[[1]]$marker$opacity)
  expect_true(length(steps) > 0)
  expect_true(all(vapply(steps, function(d) d$marker$color == PAL$blue, logical(1))))
})

test_that("DIF 的顏色：kits_dif -> kits_dif_combine -> 成長曲線和階梯；none 灰、girls 橘、boys 藍、inconsistent 黑", {
  pal <- c(none = "#7F7F7F", `favors girls` = "#D55E00", `favors boys` = "#0072B2", inconsistent = "#000000")
  b <- list(a1 = c(-1, 0), a2 = 0.5, a3 = -0.5, a4 = 0, a5 = 1)
  th <- seq(-3, 3, length.out = 300)
  sim <- function(girl, seed) {
    X <- rbind(kits_sim(th, b, prefix = "b", seed = seed), kits_sim(th, modifyList(b, girl), prefix = "g", seed = seed + 1))
    list(X = X, sex = rep(c("boy", "girl"), each = 300))
  }
  c1 <- sim(list(a2 = -0.5, a4 = 1, a5 = 2), 1)    # a2 easier for girls, a4 and a5 harder
  c2 <- sim(list(a2 = -0.5, a4 = -1, a5 = 2), 3)   # a4 the other way: inconsistent
  dif <- kits_dif_combine(kits_dif(c1$X, c1$sex), kits_dif(c2$X, c2$sex))
  expect_equal(dif$DIF[match(c("a2", "a4", "a5"), dif$item_id)], c("favors girls", "inconsistent", "favors boys"))
  expect_true("none" %in% dif$DIF)
  norm <- data.frame(score = "lan", age_m = 3:24)
  norm$P50 <- 20 + 15 * log(norm$age_m)
  norm$P10 <- norm$P50 - 5
  norm$P90 <- norm$P50 + 5
  items <- data.frame(item_id = dif$item_id, score = "lan", column = "lan", D = c(40, 50, 55, 60, 65),
                      age50 = c(4, 8, 11, 15, 20), DIF = dif$DIF, stem = paste("Item", dif$item_id))
  items <- old_item_positions(items, norm)
  built <- plotly::plotly_build(old_plot_growth(norm, items, "lan"))$x$data
  markers <- Filter(function(tr) identical(tr$mode, "markers"), built)
  got <- vapply(markers, function(tr) tr$marker$color[1], character(1))
  names(got) <- vapply(markers, function(tr) tr$legendgroup, character(1))
  expect_setequal(names(got), names(pal))
  expect_equal(toupper(got[names(pal)]), toupper(unname(pal)), ignore_attr = TRUE)
  marks <- old_ladder_marks(norm, "lan", ages = c(6, 12))
  ladder <- plotly::plotly_build(old_plot_ladder(items, "lan", marks))$x$data
  lines <- Filter(function(tr) identical(tr$mode, "lines") && !is.null(tr$legendgroup), ladder)
  got <- vapply(lines, function(tr) tr$line$color, character(1))
  names(got) <- vapply(lines, function(tr) tr$legendgroup, character(1))
  expect_equal(toupper(got[names(pal)]), toupper(unname(pal)), ignore_attr = TRUE)
})

test_that("Wright map：兩格之間只有一點空隙，刻度還在", {
  k <- kits(kits_sim(rnorm(300), list(a1 = c(-1, 0), a2 = 0.5, a3 = -0.5), seed = 1))
  lay <- plotly::plotly_build(kits:::plot.kits(k))$x$layout
  gap <- lay$xaxis2$domain[1] - lay$xaxis$domain[2]
  expect_gt(gap, 0)
  expect_lt(gap, 0.02)
  expect_equal(length(lay$xaxis2$ticktext), 3)
})
