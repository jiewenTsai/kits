# 0.20.0：kits_convert() 用管線逐條換算時，和一次給所有連結一樣決定要不要加 se_link。
# 前面有一條連結沒有重抽時，結果記下來（attr(, "links_without_se")），之後的換算不加 se_link，並給同樣的訊息。

betas_c <- list(a1 = -2, a2 = -1, a3 = c(-0.5, 0.5), a4 = 0, a5 = c(0.5, 1.5), a6 = 2, a7 = 3)
X03_c <- kits_sim(rnorm(500, 0), betas_c[1:5], seed = 1)
X36_c <- kits_sim(rnorm(500, 2), betas_c[3:7], prefix = "q", seed = 2)
item03_c <- kits_item(X03_c, scale = "lan M03")
item36_c <- kits_item(X36_c, scale = "lan M36")
anchors_c <- data.frame(item_id = c("a1", "a2", "a3", "a4", "a5"), tau = c(40, 47, 51, 55, 60))
co0_c <- kits_linking(item36_c, item03_c)
cob_c <- kits_linking(item36_c, item03_c, n_boot = 30)
d0_c <- kits_linking(item03_c, anchors_c)
db_c <- kits_linking(item03_c, anchors_c, n_boot = 30)
s36_c <- kits_subj(X36_c, item36_c)

link_cols <- c("se_link", "se_total")

# 一次給兩條連結，和用管線逐條換算（訊息不管）
convert_both_ways <- function(x, l1, l2, from = NULL) {
  once <- suppressMessages(kits_convert(x, l1, l2, from = from))
  piped <- suppressMessages(kits_convert(kits_convert(x, l1, from = from), l2))
  list(once = once, piped = piped)
}

expect_same_link_columns <- function(r) {
  has_once <- link_cols %in% names(r$once)
  has_piped <- link_cols %in% names(r$piped)
  expect_identical(has_piped, has_once)
  for (col in link_cols[has_once]) {
    expect_equal(r$piped[[col]], r$once[[col]], tolerance = 1e-12, info = col)
  }
}

test_that("分數：四種順序，管線和一次呼叫在有沒有 se_link、se_total 和它們的值上一致", {
  orders <- list(list(co0_c, db_c), list(cob_c, d0_c), list(co0_c, d0_c), list(cob_c, db_c))
  for (o in orders) {
    r <- convert_both_ways(s36_c, o[[1]], o[[2]])
    expect_same_link_columns(r)
    expect_identical(r$piped$theta, r$once$theta)
  }
  expect_false("se_link" %in% names(convert_both_ways(s36_c, co0_c, db_c)$piped))
  expect_true("se_link" %in% names(convert_both_ways(s36_c, cob_c, db_c)$piped))
})

test_that("沒有重抽、再有重抽：管線給和一次呼叫同樣的訊息，並說出前面那條連結", {
  expect_message(kits_convert(s36_c, co0_c, db_c), "link\\(s\\) 1 .*no bootstrap draws", class = "kits_message_scale")
  first <- suppressMessages(kits_convert(s36_c, co0_c))
  rec <- attr(first, "links_without_se")
  expect_identical(rec$link, "lan M36 -> lan M03")
  expect_identical(rec$reason, "no draws")
  expect_message(out <- kits_convert(first, db_c), "earlier conversion.*lan M36 -> lan M03.*no bootstrap draws",
                 class = "kits_message_scale")
  expect_false(any(link_cols %in% names(out)))
  expect_identical(attr(out, "links_without_se"), rec)
  # 整串都沒有重抽：兩種寫法都不說話
  expect_no_message(kits_convert(s36_c, co0_c, d0_c), class = "kits_message_scale")
  expect_no_message(kits_convert(first, d0_c), class = "kits_message_scale")
})

test_that("有重抽、再沒有重抽：原本的 se_link 拿掉，有訊息，之後再有重抽也不加回來", {
  first <- kits_convert(s36_c, cob_c)
  expect_true(all(link_cols %in% names(first)))
  expect_null(attr(first, "links_without_se"))
  expect_message(second <- kits_convert(first, d0_c), "link\\(s\\) 1 .*no bootstrap draws", class = "kits_message_scale")
  expect_false(any(link_cols %in% names(second)))
  expect_identical(attr(second, "links_without_se")$link, "lan M03 -> D")
})

test_that("門檻表（coef(item)）：四種順序都一致；kits_item 換算過的門檻不帶連結誤差，之後也不加", {
  orders <- list(list(co0_c, db_c), list(cob_c, d0_c), list(co0_c, d0_c), list(cob_c, db_c))
  for (o in orders) {
    r <- convert_both_ways(coef(item36_c), o[[1]], o[[2]])
    expect_same_link_columns(r)
    expect_identical(r$piped$thr, r$once$thr)
  }
  item_conv <- suppressMessages(kits_convert(item36_c, cob_c))
  expect_identical(attr(coef(item_conv), "links_without_se")$reason, "kits_item")
  expect_message(thr <- kits_convert(coef(item_conv), db_c), "applied to a kits_item", class = "kits_message_scale")
  expect_false(any(link_cols %in% names(thr)))
})

test_that("數字向量：兩種寫法的值相同，都沒有 se_link", {
  orders <- list(list(co0_c, db_c), list(cob_c, d0_c), list(cob_c, db_c))
  for (o in orders) {
    r <- convert_both_ways(c(-1, 0, 1), o[[1]], o[[2]], from = "lan M36")
    expect_identical(as.vector(r$piped), as.vector(r$once))
    expect_null(attr(r$piped, "se_link"))
    expect_null(attr(r$once, "se_link"))
  }
})

test_that("紀錄跟著表走：選列、dplyr、rbind() 留著；bind_rows() 只留第一張表的", {
  first <- suppressMessages(kits_convert(s36_c, co0_c))
  expect_false("se_link" %in% names(suppressMessages(kits_convert(first[1:5, ], db_c))))
  expect_false("se_link" %in% names(suppressMessages(kits_convert(dplyr::filter(first, theta > 0), db_c))))
  with_se <- kits_convert(s36_c, cob_c)
  stacked <- rbind(with_se, first)
  expect_identical(attr(stacked, "links_without_se")$link, "lan M36 -> lan M03")
  expect_false("se_link" %in% names(suppressMessages(kits_convert(stacked, db_c))))
  expect_warning(bound <- dplyr::bind_rows(first, with_se), "first table only",
                 class = "kits_warning_scale")
  expect_identical(attr(bound, "links_without_se")$link, "lan M36 -> lan M03")
  expect_identical(attr(purrr::list_rbind(list(with_se, first)), "links_without_se")$link,
                   "lan M36 -> lan M03")
  expect_identical(attr(vctrs::vec_rbind(with_se, first), "links_without_se")$link, "lan M36 -> lan M03")
})
