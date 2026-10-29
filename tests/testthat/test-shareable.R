# 模型物件（題目參數、連結、常模……）只存彙總的數字，可以分享：不能帶任何人次或孩子的編號

ids_inside <- function(x, ids) {
  s <- rawToChar(serialize(x, NULL, ascii = TRUE))
  ids[vapply(ids, \(i) grepl(i, s, fixed = TRUE), logical(1))]
}

test_that("題目參數物件不帶人次或孩子的編號（item$design 的列名也不行）", {
  resp <- kitm03$responses[kitm03$responses$domain == "lan", ]
  item <- suppressMessages(kits_item(kits_resp(resp), scale = "lan M03"))
  expect_length(ids_inside(item, unique(resp$visit_id)), 0)
  expect_length(ids_inside(item, unique(resp$subj_id)), 0)
  expect_identical(rownames(item$design), as.character(seq_len(nrow(item$design))))
})
