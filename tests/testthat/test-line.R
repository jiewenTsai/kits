# 直線代數（R/line.R），和換算過的分數的 summary()：合成的線、連結誤差，和 kits_convert() 相同

lan03 <- suppressMessages(kits_resp(dplyr::filter(kitm03$responses, domain == "lan")))
lan36 <- suppressMessages(kits_resp(dplyr::filter(kitm36$responses, domain == "lan")))
mot03 <- suppressMessages(kits_resp(dplyr::filter(kitm03$responses, domain == "motor")))
i03 <- suppressMessages(kits_item(lan03, scale = "lan M03"))
i36 <- suppressMessages(kits_item(lan36, scale = "lan M36"))
im03 <- suppressMessages(kits_item(mot03, scale = "motor M03"))
co <- suppressMessages(kits_linking(i36, i03, n_boot = 50))
d_link <- suppressMessages(kits_linking(list(lan = i03, motor = im03), kitm03$anchors,
                                        n_boot = 50))
s36 <- kits_subj(lan36, i36)
s36_d <- suppressMessages(s36 |> kits_convert(co) |> kits_convert(d_link))
theta <- seq(-2, 6, by = 0.25)

test_that("line algebra: apply, compose, invert, mean/sigma", {
  l1 <- list(intercept = 1.5, slope = 2)
  l2 <- list(intercept = -3, slope = 0.25)
  x <- c(-1, 0, 2.5)
  expect_identical(kits:::line_apply(l1, x), 1.5 + 2 * x)
  expect_identical(kits:::line_resid(l1, x, x), x - 1.5 - 2 * x)
  both <- kits:::line_compose(l1, l2)
  expect_identical(both, list(intercept = -3 + 0.25 * 1.5, slope = 0.5))
  expect_equal(kits:::line_apply(both, x), kits:::line_apply(l2, kits:::line_apply(l1, x)))
  back <- algo_line(l1, invert = TRUE)
  expect_equal(kits:::line_apply(back, kits:::line_apply(l1, x)), x)
  from <- c(-1, 0.3, 2, 4)
  to <- c(40, 47, 51, 60)
  ms <- kits:::algo_mean_sigma(from, to)
  expect_identical(ms, kits:::mean_sigma_moments(mean(from), sd(from), mean(to), sd(to)))
  expect_identical(ms$slope, sd(to) / sd(from))
  expect_identical(ms$intercept, mean(to) - ms$slope * mean(from))
})

test_that("summary(): the chain of kitm36 lan -> lan M03 -> D, in order", {
  conv <- summary(s36_d)$conversion
  expect_identical(conv$chain$step, 1:2)
  expect_identical(conv$chain$from, c("lan M36", "lan M03"))
  expect_identical(conv$chain$to, c("lan M03", "D"))
  expect_identical(conv$chain$intercept, c(co$intercept, d_link$intercept))
  expect_identical(conv$chain$slope, c(co$slope, d_link$slope))
  expect_identical(c(conv$intercept, conv$slope),
                   c(d_link$intercept + d_link$slope * co$intercept, d_link$slope * co$slope))
  expect_identical(c(conv$from_scale, conv$to_scale), c("lan M36", "D"))
  # one kits_convert(co, d_link) gives the same chain
  one <- suppressMessages(kits_convert(s36, co, d_link))
  expect_identical(summary(one)$conversion[c("chain", "intercept", "slope")],
                   conv[c("chain", "intercept", "slope")])
})

test_that("summary(): the composed line and its se_link agree with kits_convert()", {
  conv <- summary(s36_d)$conversion
  by_convert <- kits_convert(data.frame(theta = theta), co, d_link, from = "lan M36")
  composed <- kits:::line_apply(conv, theta)
  expect_equal(composed, by_convert$theta, tolerance = 1e-12)
  expect_equal(kits:::line_apply(conv, s36$theta), s36_d$theta, tolerance = 1e-12)
  # se_link at the mean theta of the starting scale
  expect_equal(conv$theta_mean, mean(s36$theta), tolerance = 1e-12)
  at_mean <- kits_convert(data.frame(theta = conv$theta_mean), co, d_link, from = "lan M36")
  expect_identical(conv$se_link, at_mean$se_link)
  expect_identical(conv$value_mean, at_mean$theta)
  out <- capture.output(print(summary(s36_d)))
  expect_true(any(grepl("Conversion lan M36 -> D", out)))
  expect_true(any(grepl("composed: new value = .* x old value", out)))
  expect_true(any(grepl("Chained links compose as lines", out)))
  expect_true(any(grepl("linking error at the mean theta", out)))
})

test_that("summary(): links without draws give no se_link; unconverted scores are unchanged", {
  co0 <- suppressMessages(kits_linking(i36, i03))
  conv0 <- summary(suppressMessages(kits_convert(s36, co0)))$conversion
  expect_true(is.na(conv0$se_link))
  expect_identical(conv0$chain$to, "lan M03")
  plain <- summary(s36)
  expect_null(plain$conversion)
  expect_identical(names(plain), c("object", "diagnostics", "digits"))
  expect_identical(capture.output(print(plain)),
                   capture.output(print(kits:::kits_summary(s36))))
})

test_that("summary(): stacked scores of two calibrations, one conversion each", {
  s03_d <- suppressMessages(kits_subj(lan03, i03) |> kits_convert(d_link))
  conv <- summary(rbind(s03_d, s36_d))$conversion
  expect_length(conv, 2)
  expect_identical(map_chr(conv, "from_scale"), c("lan M03", "lan M36"))
  expect_identical(conv[[2]]$chain, summary(s36_d)$conversion$chain)
})

test_that("kits_vertical(): the linking error of a mean comes from the links of the scores", {
  s03 <- suppressMessages(kits_subj(lan03, i03) |> kits_convert(d_link))
  v <- suppressWarnings(kits_vertical(rbind(s03, s36_d)))
  g <- dplyr::filter(as.data.frame(s36_d), wave_m == 48)
  line <- kits:::line_compose(co, d_link)
  x0 <- algo_line(line, mean(g$theta), invert = TRUE)
  expect_identical(v$se_link[v$cohort == "M36" & v$wave_m == 48],
                   sqrt(kits:::link_moments(list(co, d_link), x0)$var))
  # the back-solve (scores without their links) agrees up to rounding
  plain <- as.data.frame(rbind(s03, s36_d))
  plain <- plain[c("cohort", "wave_m", "theta", "se", "se_link")]
  v_plain <- suppressWarnings(kits_vertical(plain))
  expect_equal(v_plain$se_link, v$se_link, tolerance = 1e-12)
})
