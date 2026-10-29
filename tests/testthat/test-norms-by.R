# 照組的常模（kits_norms(..., by = )）：每組和只用這組的 kits_norms() 完全相同；同齡 z 用自己那組的曲線；
# 圖一組一格。

set.seed(11)
age <- runif(1200, 3, 30)
sex <- rep(c("boy", "girl"), 600)
domain <- rep(c("lan", "motor"), each = 600)
D <- 20 + 15 * log(age) + ifelse(sex == "girl", 1, 0) + ifelse(domain == "motor", 3, 0) +
  2.5 * rt(1200, 6)
by_sex <- kits_norms(D, age, ages = 3:30, degrees = c(3, 2, 1, 1), by = data.frame(sex = sex))

# 一組的列，拿掉組的欄
group_rows <- function(norm, keep) {
  as.data.frame(norm)[keep, setdiff(names(norm), attr(norm, "by"))]
}

test_that("每組的常模和只用這組的 kits_norms() 完全相同（容許誤差 0）；組的欄在前", {
  for (g in c("boy", "girl")) {
    alone <- kits_norms(D[sex == g], age[sex == g], ages = 3:30, degrees = c(3, 2, 1, 1))
    mine <- group_rows(by_sex, by_sex$sex == g)
    expect_equal(mine, as.data.frame(alone), tolerance = 0, ignore_attr = TRUE)
  }
  expect_identical(names(by_sex)[1:2], c("sex", "age_m"))
  expect_s3_class(by_sex, "kits_norms")
  expect_identical(attr(by_sex, "by"), "sex")
  expect_identical(kits_info(by_sex)$settings$by, "sex")
  expect_length(attr(by_sex, "groups"), 2)
  expect_true(all(grepl("^sex = (boy|girl): ", attr(by_sex, "diagnostics")$check)))
})

test_that("kits_subj 的 by 是欄名：領域和性別兩欄；和數字的 by 相同", {
  s <- structure(data.frame(visit_id = paste0("p", 1:1200), theta = D, age_m = age, sex = sex,
                            domain = domain),
                 class = c("kits_subj", "kits_table", "data.frame"), scale = "D")
  two <- kits_norms(s, ages = 3:30, degrees = c(3, 2, 1, 1), by = c("domain", "sex"))
  expect_identical(names(two)[1:3], c("domain", "sex", "age_m"))
  num <- kits_norms(D, age, ages = 3:30, degrees = c(3, 2, 1, 1),
                    by = data.frame(domain = domain, sex = sex))
  expect_equal(as.data.frame(two), as.data.frame(num), tolerance = 0, ignore_attr = TRUE)
  keep <- domain == "motor" & sex == "girl"
  alone <- kits_norms(D[keep], age[keep], ages = 3:30, degrees = c(3, 2, 1, 1))
  mine <- group_rows(two, two$domain == "motor" & two$sex == "girl")
  expect_equal(mine, as.data.frame(alone), tolerance = 0, ignore_attr = TRUE)
  expect_error(kits_norms(s, ages = 3:30, by = "cohort"), "cohort", class = "kits_error_input")
  expect_error(kits_norms(D, age, ages = 3:30, by = "sex"), "data.frame",
               class = "kits_error_input")
})

test_that("一組太小時，錯誤說出是哪一組；by 有 NA 時停下來", {
  small <- c(rep("big", 1190), rep("tiny", 10))
  err <- expect_error(kits_norms(D, age, ages = 3:30, degrees = c(3, 2, 1, 1),
                                 by = data.frame(g = small)),
                      class = "kits_error_data")
  expect_match(conditionMessage(err), "^g = tiny: ")
  expect_error(kits_norms(D, age, ages = 3:30, by = data.frame(sex = replace(sex, 3, NA))), "NA",
               class = "kits_error_data")
  expect_error(kits_norms(D, age, ages = 3:30, by = data.frame(sex = sex[-1])), "one row per score",
               class = "kits_error_input")
})

test_that("kits_norm_z()：每個孩子用自己那組的曲線；沒有組的欄時說清楚；沒有常模的組是 NA", {
  z <- kits_norm_z(D, age, by_sex, by = data.frame(sex = sex))
  for (g in c("boy", "girl")) {
    curve <- group_rows(by_sex, by_sex$sex == g)
    expect_identical(z[sex == g], kits_norm_z(D[sex == g], age[sex == g], curve))
  }
  expect_error(kits_norm_z(D, age, by_sex), "by = data.frame", class = "kits_error_input")
  s <- structure(data.frame(visit_id = paste0("p", 1:1200), theta = D, age_m = age, sex = sex),
                 class = c("kits_subj", "kits_table", "data.frame"))
  expect_identical(kits_norm_z(s, by_sex)$z, z)
  expect_error(kits_norm_z(s[c("visit_id", "theta", "age_m")], by_sex), "sex",
               class = "kits_error_input")
  other <- replace(sex, 1:2, "other")
  expect_message(z2 <- kits_norm_z(D, age, by_sex, by = data.frame(sex = other)), "sex = other",
                 class = "kits_message_data")
  expect_true(all(is.na(z2[1:2])))
})

test_that("plot()：一個分數 x 一組一格；分數也是組時一列一個分數；孩子只畫在一個常模後面", {
  built <- plotly::plotly_build(plot(by_sex))
  expect_length(built$x$data, 10)
  colours <- unique(vapply(built$x$data, function(d) d$line$color, character(1)))
  expect_identical(colours, c("#0072B2", "#D55E00"))
  notes <- vapply(built$x$layout$annotations, function(a) a$text, character(1))
  expect_true(all(c("Boys", "Girls") %in% notes))
  two <- kits_norms(D, age, ages = 3:30, degrees = c(3, 2, 1, 1),
                    by = data.frame(domain = domain, sex = sex))
  expect_length(plotly::plotly_build(plot(two))$x$data, 20)
  motor <- dplyr::filter(two, domain == "motor")
  expect_length(plotly::plotly_build(plot(motor, pct = c(10, 50, 90)))$x$data, 6)
  children <- data.frame(theta = D, age_m = age, sex = sex)
  with_children <- plotly::plotly_build(plot(by_sex, scores = children))$x$data
  expect_length(with_children, 12)   # 每一格：孩子的格子 + 5 條線
  expect_error(plot(by_sex, scores = children[c("theta", "age_m")]), "sex",
               class = "kits_error_input")
  single <- kits_norms(D, age, ages = 3:30, degrees = c(3, 2, 1, 1))
  expect_s3_class(plot(single, scores = children), "plotly")
  expect_error(plot(single, d = D, age_m = age), "does not have", class = "kits_error_input")
})
