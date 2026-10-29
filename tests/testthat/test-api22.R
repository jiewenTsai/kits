# kits 0.22.0 的整理（0.23.0 拿掉了被取代的函數）：scores |> kits_subj_table()、summary(link, theta = )、
# 圖收物件、kits_vertical() 的 P25、P75、泛型的 kits_information()、kits_expected_score()。

betas <- list(a1 = -2, a2 = -1, a3 = c(-0.5, 0.5), a4 = 0, a5 = c(0.5, 1.5), a6 = 2, a7 = 3)
X03 <- kits_sim(rnorm(500, 0), betas[1:5], seed = 1)
X36 <- kits_sim(rnorm(500, 2), betas[3:7], prefix = "q", seed = 2)
item03 <- kits_item(X03, scale = "lan M03")
item36 <- kits_item(X36, scale = "lan M36")
anchors <- data.frame(item_id = c("a1", "a2", "a3", "a4", "a5"), tau = c(40, 47, 51, 55, 60))
co <- kits_linking(item36, item03, n_boot = 50)
d_link <- kits_linking(item03, anchors, n_boot = 50)

test_that("kits_link_screen()：自己的點才給 from_scale、to_scale；kits_link_points() 的點已經有量尺", {
  pts <- kits_link_points(item36, item03)
  expect_error(kits_link_screen(pts, from_scale = "x"), "only with your own points",
               class = "kits_error_input")
  own <- data.frame(id = letters[1:5], item = letters[1:5], from = c(-2, -1, 0, 1, 2),
                    to = c(-1, 0, 1, 2, 3.2))
  expect_error(kits_link_screen(own, from_scale = "A", to_scale = "A"), class = "kits_error_scale")
  link <- kits_link_line(kits_link_screen(own, to_scale = "D"))
  expect_true(is.na(link$from_scale))
  expect_identical(link$to_scale, "D")
})

# ---- kits_subj_table(item, data, wave_m) ---------------------------------------------------------

test_that("kits_subj_table(item) 不用 data：每個題本、每個原始分的 WLE，和 kits_subj() 相同", {
  tab <- kits_subj_table(item03)
  s <- kits_subj(X03, item03)
  both <- dplyr::inner_join(as.data.frame(s), as.data.frame(tab),
                            by = c("pattern_id", "raw"))
  expect_equal(both$theta.x, both$theta.y)
  expect_equal(both$se.x, both$se.y)
  expect_identical(attr(tab, "scale"), "lan M03")
  expect_match(attr(tab, "info")$call, "^kits_subj_table\\(")
})

# ---- 圖收物件 ----------------------------------------------------------------------------------------

set.seed(5)
s_lan <- dplyr::tibble(cohort = rep(c("M03", "M36"), c(600, 400)),
                       wave_m = c(rep(c(12, 24, 36), each = 200), rep(c(36, 48), each = 200))) |>
  dplyr::mutate(age_m = wave_m + runif(dplyr::n(), -1, 1),
                theta = 40 + 8 * log(wave_m) + rnorm(dplyr::n(), 0, 4))
v_lan <- kits_vertical(s_lan)
v_motor <- kits_vertical(dplyr::mutate(s_lan, theta = theta + 3))

test_that("kits_vertical() 多了 P25、P75（quantile() 的 type 7）；其他的數字不變", {
  g <- s_lan$theta[s_lan$cohort == "M03" & s_lan$wave_m == 24]
  row <- v_lan[v_lan$cohort == "M03" & v_lan$wave_m == 24, ]
  expect_equal(row$P25, unname(quantile(g, 0.25)))
  expect_equal(row$P75, unname(quantile(g, 0.75)))
  expect_equal(row$median, median(g))
  expect_identical(names(v_lan)[5:9], c("mean", "sd", "median", "P25", "P75"))
})

test_that("幾個領域的 kits_vertical() 用 list_rbind(names_to = \"domain\") 疊起來；量尺不同時停下來", {
  both <- purrr::list_rbind(list(lan = v_lan, motor = v_motor), names_to = "domain")
  expect_s3_class(both, "kits_vertical")
  expect_identical(unique(both$domain), c("lan", "motor"))
  expect_equal(nrow(both), 2 * nrow(v_lan))
  on_d <- kits_vertical(structure(as.data.frame(s_lan), scale = "D"))
  expect_identical(attr(on_d, "scale"), "D")
  on_other <- kits_vertical(structure(as.data.frame(s_lan), scale = "lan M03"))
  expect_error(purrr::list_rbind(list(lan = on_d, motor = on_other), names_to = "domain"),
               class = "kits_error_scale")
})

test_that("plot(kits_vertical)：type = \"mean\" 或 \"median\"；一個領域一格", {
  both <- purrr::list_rbind(list(lan = v_lan, motor = v_motor), names_to = "domain")
  for (p in list(plot(v_lan), plot(v_lan, type = "median"), plot(both),
                 plot(both, type = "median"))) {
    expect_s3_class(p, "plotly")
  }
  built <- plotly::plotly_build(plot(both, type = "median"))
  panels <- unique(vapply(built$x$data, function(d) d$xaxis, character(1)))
  expect_length(panels, 2)
  expect_equal(as.vector(built$x$data[[1]]$y), v_lan$median[v_lan$cohort == "M03"])
  notes <- vapply(built$x$layout$annotations, function(a) a$text, character(1))
  expect_true(all(c("Language", "Motor") %in% notes))
  expect_error(plot(v_lan, type = "mode"), "type", class = "kits_error_input")
  one <- kits_vertical(data.frame(theta = rnorm(40), wave_m = rep(1:2, 20)))
  expect_s3_class(plot(one, type = "median"), "plotly")
})

test_that("plot(kits_vertical)：月齡的註記、hover 的月齡、cohorts 的順序、age_ticks、y_range", {
  late <- dplyr::mutate(s_lan, age_m = ifelse(wave_m == 48, 50.5, wave_m))
  built <- plotly::plotly_build(plot(kits_vertical(late), type = "median",
                                     cohorts = c("M36", "M03"),
                                     age_ticks = c(12, 36, 50), y_range = c(30, 90)))
  title <- built$x$layout$title
  if (is.list(title)) title <- title$text
  expect_match(title, "48-month wave 50.5")
  expect_no_match(title, "36-month wave")
  expect_identical(built$x$data[[1]]$name, "M36")
  expect_match(built$x$data[[1]]$hovertext[1], "actual median age")
  expect_equal(built$x$layout$xaxis$tickvals, c(12, 36, 50))
  expect_equal(built$x$layout$yaxis$range, c(30, 90))
  expect_error(plot(v_lan, cohorts = "M99"), "M99", class = "kits_error_input")
})

# ---- 泛型、套件裡面不叫被取代的函數 -----------------------------------------------------------------------------

test_that("kits_information()、kits_expected_score() 是泛型：kits_item、kits、kits_fit、參數表的結果相同", {
  k <- kits(X03, scale = "lan M03")
  th <- c(-2, 0, 2)
  for (f in list(kits_information, kits_expected_score)) {
    ref <- f(k$item)(th)
    expect_equal(f(k)(th), ref)
    expect_equal(f(k$item$fit)(th), ref)
    expect_equal(f(coef(k$item))(th), ref)
  }
  expect_true(all(c("kits_information.kits_item", "kits_information.kits") %in%
                    as.character(utils::methods("kits_information"))))
  expect_error(kits_information(k$item, step = 1), "does not have", class = "kits_error_input")
})

test_that("0.22.0 停用的函數拿掉了；平常的用法沒有 kits_warning_deprecated", {
  gone <- c("kits_age50", "kits_mean_sigma", "kits_link_se", "kits_link_screening",
            "kits_wle_table",
            "kits_plot_cohort", "kits_plot_norms", "kits_item_positions", "kits_ladder_marks",
            "kits_plot_growth", "kits_plot_ladder",
            "kits_inter",   # 0.23.1 拿掉
            "kits_lrtest", "kits_age_pass",   # 0.25.0 併進 kits_dif()、kits_item_table(data = )
            "kits_code", "kits_equations", "kits_rtmb_code",   # 0.25.0 拿掉
            "kits_raw_equate")   # 0.25.0：kits_subj_table() 的 WLE 加 kits_expected_score()
  expect_false(any(gone %in% getNamespaceExports("kits")))
  # 0.23.0：加 kits_esf、kits_item_restscore、kits_lrtest、kits_pcm_prob；0.23.1：拿掉 kits_inter；
  # 0.25.0：kits_algorithm() 取代 kits_code()、kits_equations()、kits_rtmb_code()；拿掉 kits_raw_equate()
  expect_length(grep("^algo_", getNamespaceExports("kits"), invert = TRUE, value = TRUE), 42)
  expect_length(grep("^algo_", getNamespaceExports("kits"), value = TRUE), 12)   # 0.25.0：演算法層
  no_old <- function(expr) {
    stop_old <- function(w) stop("deprecated warning: ", conditionMessage(w))
    withCallingHandlers(expr, kits_warning_deprecated = stop_old)
  }
  wave <- rep(c(12, 24), each = 250)
  scores <- kits_subj(X03, item03, persons = data.frame(visit_id = rownames(X03), wave_m = wave))
  expect_no_error(no_old(kits_subj_table(scores)))
  expect_no_error(no_old(kits_precision(item03, kits_subj_table(scores))))
  expect_no_error(no_old(suppressMessages(kits_drift(X03, item03))))
  expect_no_error(no_old(summary(co, theta = 0)))
  expect_no_error(no_old(plot(v_lan, type = "median")))
})
