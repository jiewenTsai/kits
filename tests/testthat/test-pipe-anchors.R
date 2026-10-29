# 0.17.0：kits_convert() 用管線一行一條連結（和一次給所有連結相同）；kits_anchors() 的錨題表

betas_p <- list(a1 = -2, a2 = -1, a3 = c(-0.5, 0.5), a4 = 0, a5 = c(0.5, 1.5), a6 = 2, a7 = 3)
X03_p <- kits_sim(rnorm(500, 0), betas_p[1:5], seed = 1)
X36_p <- kits_sim(rnorm(500, 2), betas_p[3:7], prefix = "q", seed = 2)
item03_p <- kits_item(X03_p, scale = "lan M03")
item36_p <- kits_item(X36_p, scale = "lan M36")
co_p <- kits_linking(item36_p, item03_p, n_boot = 50)
anchors_p <- data.frame(item_id = c("a1", "a2", "a3", "a4", "a5"), tau = c(40, 47, 51, 55, 60))
d_link_p <- kits_linking(item03_p, anchors_p, n_boot = 50)

test_that("管線：一行一條連結，和 kits_convert(x, co, d_link) 相同（分數、門檻、數字向量）", {
  s36 <- kits_subj(X36_p, item36_p)
  piped <- s36 |>
    kits_convert(co_p) |>
    kits_convert(d_link_p)
  once <- kits_convert(s36, co_p, d_link_p)
  expect_identical(piped$theta, once$theta)
  expect_identical(piped$se, once$se)
  expect_equal(piped$se_link, once$se_link, tolerance = 1e-14)
  expect_equal(piped$se_total, once$se_total, tolerance = 1e-14)
  expect_identical(attr(piped, "scale"), "D")
  expect_identical(attr(piped, "scale"), attr(once, "scale"))

  thr_piped <- coef(item36_p) |>
    kits_convert(co_p) |>
    kits_convert(d_link_p)
  thr_once <- kits_convert(coef(item36_p), co_p, d_link_p)
  expect_identical(thr_piped$thr, thr_once$thr)
  expect_identical(thr_piped$thr_se, thr_once$thr_se)
  expect_equal(thr_piped$se_link, thr_once$se_link, tolerance = 1e-14)
  expect_equal(thr_piped$se_total, thr_once$se_total, tolerance = 1e-14)
  expect_identical(attr(thr_piped, "scale"), attr(thr_once, "scale"))

  v_piped <- c(-1, 0, 1) |>
    kits_convert(co_p, from = "lan M36") |>
    kits_convert(d_link_p)
  v_once <- kits_convert(c(-1, 0, 1), co_p, d_link_p, from = "lan M36")
  expect_identical(as.vector(v_piped), as.vector(v_once))
  expect_identical(attr(v_piped, "scale"), "D")
  expect_identical(attr(v_once, "scale"), "D")
})

test_that("管線：順序錯了還是停下來（kits_error_scale），數字向量也是", {
  s36 <- kits_subj(X36_p, item36_p)
  expect_error(s36 |> kits_convert(d_link_p), class = "kits_error_scale")
  expect_error(s36 |> kits_convert(co_p) |> kits_convert(co_p), class = "kits_error_scale")
  expect_error(c(-1, 0, 1) |> kits_convert(co_p, from = "lan M36") |> kits_convert(co_p), class = "kits_error_scale")
  expect_error(coef(item36_p) |> kits_convert(d_link_p), class = "kits_error_scale")
})

test_that("kits_anchors()：兩個 key 的列數和欄；gsed2510 有警告；兩個 key 疊在一起時 kits_linking() 停下來", {
  a <- kits_anchors()
  expect_s3_class(a, "data.frame")
  expect_identical(class(a), "data.frame")
  expect_equal(nrow(a), 29)
  expect_equal(names(a), c("item_id", "gsed_item", "key", "tau", "note", "domain"))
  expect_true(all(a$key == "gsed2406"))
  expect_false(anyDuplicated(a$item_id) > 0)
  expect_identical(kits_anchors("gsed2406"), a)
  expect_warning(b <- kits_anchors("gsed2510"), class = "kits_warning_anchors")
  expect_equal(nrow(b), 8)
  expect_equal(names(b), names(a))
  expect_true(all(b$item_id %in% a$item_id))
  expect_error(kits_anchors("gsed2212"), class = "kits_error_input")
  both <- rbind(a, b)
  ids <- a$item_id[1:6]
  betas <- stats::setNames(as.list(seq(-2, 2, length.out = 6)), ids)
  X <- kits_sim(rnorm(500, 0, 1.5), betas, seed = 3)
  item <- kits_item(X, scale = "sim M03")
  expect_error(kits_linking(item, both), "choose one key", class = "kits_error_input")
})

test_that("kits_anchors()：用 KIT 題號的模擬資料連到 D；連結記下 key", {
  a <- kits_anchors()
  ids <- a$item_id[1:8]
  betas <- stats::setNames(as.list((a$tau[1:8] - 55) / 10), ids)
  X <- kits_sim(rnorm(800, 0, 1.5), betas, seed = 4)
  item <- kits_item(X, scale = "sim M03")
  d_link <- kits_linking(item, a)
  expect_equal(d_link$n, 8)
  expect_equal(d_link$to_scale, "D")
  expect_equal(d_link$settings$key, "gsed2406")
  expect_equal(kits_info(d_link)$settings$key, "gsed2406")
  expect_output(print(d_link), "dscore key gsed2406")
  manual <- kits_linking(item, data.frame(item_id = a$item_id, tau = a$tau))
  expect_identical(c(d_link$intercept, d_link$slope), c(manual$intercept, manual$slope))
  expect_null(manual$settings$key)
  expect_gt(d_link$slope, 5)                                  # 1 logit is about 10 D here
})
