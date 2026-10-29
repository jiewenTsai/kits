# 每題 p 通過的月齡（內部用，給 kits_item_table(data = )）：每題配通過對月齡的 logistic 迴歸。
# pass：矩陣（一列一個人次、一欄一題；NA 是沒施測）或一個向量；分數 >= step 是通過（step 一題一個或共用一個），
# 邏輯值照用。
# 回傳一題一列：item_id（向量時 NA）、每個 p 一欄（p50、p75……）、separated（完全分離；沒配迴歸時 NA）、
# age_last_fail、age_first_pass（沒通過的最大月齡、通過的最小月齡）；attr "diagnostics" 記完全分離的題。
# 沒跨過 p、不隨月齡上升或不到 2 個月齡時是 NA。
age_pass <- function(pass, age_m, wave_m, step = 1, p = 0.5) {
  pass <- as.matrix(pass)
  if (is.numeric(pass)) pass <- sweep(pass, 2, rep_len(step, ncol(pass)), ">=")
  items <- if (is.null(colnames(pass))) rep(NA_character_, ncol(pass)) else colnames(pass)
  rows <- map(seq_len(ncol(pass)), \(j) age_pass_one(pass[, j], age_m, wave_m, p))
  ages <- matrix(unlist(map(rows, "age")), ncol = length(p), byrow = TRUE, dimnames = list(NULL, sprintf("p%g", 100 * p)))
  out <- data.frame(item_id = items, ages, check.names = FALSE) |>
    mutate(separated = map_lgl(rows, "separated"),
           age_last_fail = map_dbl(rows, "last_fail"),
           age_first_pass = map_dbl(rows, "first_pass"))
  out <- kits_table(out)
  attr(out, "diagnostics") <- diag_table(list(separation_row(out)))
  out
}

# 診斷：完全分離的題（年齡只知道在最後一個沒通過和第一個通過之間）
separation_row <- function(out) {
  separated <- out$item_id[out$separated %in% TRUE]
  n <- length(separated)
  status <- if (n > 0) "note" else "ok"
  message <- if (n > 0) {
    paste0("perfect separation (the age is only known between age_last_fail and age_first_pass): ",
           some_ids(separated, 10))
  } else {
    "none"
  }
  diag_row("items separated by age", n, status, message)
}

# 一題在 p 通過的月齡（p 可以是向量，迴歸只配一次）：age（每個 p 一個）、separated、last_fail、first_pass
age_pass_one <- function(pass, age_m, wave_m, p = 0.5) {
  none <- list(age = rep(NA_real_, length(p)), separated = NA, last_fail = NA_real_, first_pass = NA_real_)
  d <- data.frame(pass = as.numeric(pass), age_m = age_m, wave_m = wave_m) |>
    filter(!is.na(pass), is.finite(age_m), !is.na(wave_m))
  if (nrow(d) == 0) {
    return(none)
  }
  none$last_fail <- max_or_na(d$age_m[d$pass == 0])
  none$first_pass <- min_or_na(d$age_m[d$pass == 1])
  pass <- d$pass
  age_m <- d$age_m
  if (length(unique(age_m)) < 2) {
    return(none)
  }
  by_wave <- summarise(d, rate = mean(pass), .by = wave_m)$rate   # 每個波次的通過率
  crosses <- max(by_wave) >= p & min(by_wave) <= p
  if (!any(crosses)) {
    return(none)
  }
  fit <- logistic_fit(pass, age_m)
  none$separated <- fit$separated
  slope <- stats::coef(fit$glm)[2]
  if (!is.finite(slope) || slope <= 0) {
    return(none)
  }
  age <- unname((stats::qlogis(p) - stats::coef(fit$glm)[1]) / slope)
  age[!crosses | age < min(age_m) | age > max(age_m)] <- NA_real_
  c(list(age = age), none[-1])
}

# 通過對月齡的 logistic 迴歸。glm() 對完全分離的兩個警告（沒有收斂、配適的機率到 0 或 1）不印出來，
# 改記成 separated（其他警告照常發出）。判斷和 glm.fit 相同：機率在 10 x 機器精度以內就算到了 0 或 1
logistic_fit <- function(pass, age_m) {
  separation_warnings <- gettext(c("glm.fit: algorithm did not converge",
                                   "glm.fit: fitted probabilities numerically 0 or 1 occurred"),
                                 domain = "R-stats")   # 翻譯過的訊息也認得
  fit <- withCallingHandlers(stats::glm(pass ~ age_m, family = stats::binomial), warning = function(w) {
    if (conditionMessage(w) %in% separation_warnings) invokeRestart("muffleWarning")
  })
  eps <- 10 * .Machine$double.eps
  at_bound <- any(fit$fitted.values < eps | fit$fitted.values > 1 - eps)
  list(glm = fit, separated = !fit$converged || at_bound)
}

max_or_na <- function(x) if (length(x) == 0) NA_real_ else max(x)
min_or_na <- function(x) if (length(x) == 0) NA_real_ else min(x)
