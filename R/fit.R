# 每個作答在 PCM 下的期望值 E 與變異數 W（內部用，給 kits_subj()、kits_itemfit()）。
# resp：visit_id、item_id、item_score；abilities：visit_id、theta（每人一列）。
pcm_moments <- function(resp, parms, abilities) {
  beta_of <- betas_by_item(parms)
  resp |>
    left_join(abilities, by = "visit_id") |>
    mutate(item_moments(theta, beta_of[[item_id[1]]]), .by = item_id)
}

# 一題的 E = sum_k k P_k 和 W = sum_k (k - E)^2 P_k（中央動差），一列一個作答；來自 algo_pcm()
item_moments <- function(theta, beta) {
  m <- algo_pcm(theta, beta)
  data.frame(E = m$E, W = m$Var)
}

# infit、outfit 均方（by 是 visit_id 或 item_id）：infit = sum (x - E)^2 / sum W；outfit = (x - E)^2 / W 的平均
mean_squares <- function(r, by) {
  r |>
    summarise(n = n(), infit = sum((item_score - E)^2) / sum(W), outfit = mean((item_score - E)^2 / W), .by = {{ by }})
}


#' Item fit: infit and outfit mean squares, and the item-rest correlation
#'
#' For each item, the infit and outfit mean squares of the residuals under the PCM (expected value
#' 1; above 1 more random, below 1 more predictable than the model expects) and the item-rest
#' correlation `rir` (weighted mean over item patterns). By default the residuals are conditional
#' on the raw score within the item pattern (Christensen, Kreiner & Mesbah, 2013), so no ability
#' estimate is needed; persons with a zero or a perfect score are left out.
#'
#' Details: tutorial_algorithms/05_checks.html.
#'
#' @inheritParams kits_subj
#' @param type `"conditional"` (the default; given the raw score within the item pattern) or
#'   `"wle"` (at each person's WLE; biased downward on short forms). Recorded in `attr(x, "type")`.
#' @param ... Not used: an argument the function does not have stops with an error (class
#'   `kits_error_input`).
#' @return A data.frame (class `kits_itemfit`) with one row per item: `item_id`, `n` (responses
#'   used), `pvalue` (the mean score as a share of the maximum), `rir`, `infit`, `outfit`, and
#'   `flag`, against the range 0.6 to 1.4 (Wright & Linacre, 1994):
#'   * `"underfit"`: infit or outfit above 1.4 (it wins when an item is also below 0.6).
#'   * `"overfit"`: otherwise, infit or outfit below 0.6.
#'   * `"ok"`: both within 0.6 to 1.4, the limits included; `NA` when infit and outfit are NA.
#' @references
#' Christensen, K. B., Kreiner, S., & Mesbah, M. (Eds.). (2013). *Rasch models in health*.
#' ISTE and Wiley.
#'
#' Wright, B. D., & Linacre, J. M. (1994). Reasonable mean-square fit values.
#' *Rasch Measurement Transactions, 8*(3), 370.
#' @examples
#' betas <- list(a1 = c(-1, 0), a2 = c(0, 1), a3 = 0.5, a4 = -0.5, a5 = 0)
#' X <- kits_sim(seq(-3, 3, length.out = 400), betas, seed = 1)
#' X[, "a5"] <- sample(0:1, 400, replace = TRUE)   # a5 is pure noise
#' item <- kits_item(X)
#' fit <- kits_itemfit(X, item)
#' fit
#' identical(fit$infit, item$itemfit$infit)   # kits_item() keeps the same table
#'
#' # At each person's WLE: lower on short forms
#' kits_itemfit(X, item, type = "wle")
#' @export
kits_itemfit <- function(data, item, type = c("conditional", "wle"), ...) {
  no_extra_args(..., fun = "kits_itemfit")
  type <- match_option(type, c("conditional", "wle"), "type")
  call <- call_text(sys.call())
  resp <- resp_long(data, gaps = FALSE)
  itemfit_table(resp, item, call, type = type)
}

# 題目適配的表（內部用，給 kits_itemfit() 和 kits_item() 的 item$itemfit）：一列一題的彙總，沒有個人的值。
# resp：長格式的作答；call：記在 info 的呼叫；data：resp 的指紋（kits_item() 已經算過，不用再算一次）；
# type：conditional（給定題本內的原始分）或 wle（每個人在他的 WLE）
# sets：person_sets(resp)，kits_item() 算過時直接傳進來（同樣的資料只算一次）
itemfit_table <- function(resp, item, call, data = data_fingerprint(resp), type = "conditional",
                          sets = person_sets(resp)) {
  parms <- param_frame(item)
  need_items(resp, parms)
  abort_if(length(unique(resp$item_id)) < 2, "Item fit needs at least 2 items.")
  moments <- if (type == "conditional") {
    cond_moments(resp, parms, sets)
  } else {
    wle_moments(resp, parms)
  }
  inner <- unique(moments$visit_id)
  abort_if(length(inner) == 0,
           "Every person has a zero or a perfect score, so the fit cannot be computed.", class = "kits_error_fit")
  fit <- mean_squares(moments, item_id)
  tia <- item_totals(resp, sets) |>
    summarise(pvalue = booklet_mean(pvalue, n_persons), rir = booklet_mean(rir, n_persons), .by = item_id)
  out <- data.frame(item_id = sort(unique(resp$item_id))) |>
    left_join(fit, by = "item_id") |>
    left_join(tia, by = "item_id") |>
    mutate(n = coalesce(n, 0L), flag = itemfit_flag(infit, outfit)) |>
    select(item_id, n, pvalue, rir, infit, outfit, flag) |>
    kits_table("kits_itemfit")
  source <- source_of(item)
  same <- !is.null(source) && identical(source$scores_md5, data$scores_md5)
  msg <- if (is.null(source)) {
    "no record of the data the item parameters came from"
  } else if (same) {
    "the item parameters were calibrated on these data"
  } else {
    sprintf("the item parameters were calibrated on other data (%d person-waves, %d items); fine when intended",
            source$n_rows, length(source$items))
  }
  same_status <- if (same) "ok" else "note"
  attr(out, "diagnostics") <- diag_table(list(
    diag_row("item parameters from these data", as.numeric(same), same_status, msg),
    diag_row("persons used", length(inner), "ok", sprintf("%d person-waves without a zero or perfect score", length(inner)))))
  attr(out, "info") <- make_info(call, list(scale = beta_scale_of(item), type = type), data, item_data = source)
  attr(out, "type") <- type
  out
}

# 每個作答在 WLE 的 E、W（type = "wle"）；0 分和滿分的人次不算
wle_moments <- function(resp, parms) {
  subj <- person_sets(resp, parms) |>
    add_wle(parms)   # 每個人次的原始分、滿分和 WLE：kits_subj() 的程式（不用它的個人適配、診斷）
  inner <- filter(subj, raw > 0, raw < raw_max)
  resp |>
    filter(visit_id %in% inner$visit_id) |>
    pcm_moments(parms, select(inner, visit_id, theta))
}

# 每個作答的條件 E、W（type = "conditional"）：給定這個人次在他的題本（題目組合 b）裡的原始分 r，
# 這題分數的條件分布（Christensen, Kreiner & Mesbah, 2013；iarm::out_infit()，這裡在每個題本內）：
#   P(X_i = x | r, b) = e_ix gamma^(-i)_{b, r - x} / gamma_{b, r}，e_ix = exp(-sum_{g <= x} beta_ig)
# gamma 是題本的基本對稱函數（CML 的多項式，algo_esf()），gamma^(-i) 是拿掉第 i 題的。不需要 theta。
# 給定 r 時這題只可能有一個分數（0 分、滿分、只有一題的題本）：條件變異數是 0，這些作答不算（和 iarm 相同）
# sets：每個人次的題目組合和原始分（person_sets(resp)；kits_item() 只算一次，傳進來）
cond_moments <- function(resp, parms, sets = person_sets(resp)) {
  beta_of <- betas_by_item(parms)
  tables <- map(unique(sets$set), \(key) cond_table(key, beta_of)) |>
    list_rbind()
  resp |>
    left_join(select(sets, visit_id, set, raw), by = "visit_id") |>
    left_join(tables, by = c("set", "raw", "item_id")) |>
    filter(free) |>
    select(-set, -raw, -free)
}

# 一個題本（題目組合）：每個原始分 r、每一題的條件 E、W；free = 給定 r 時這題可能有兩個以上的分數。
# gamma、gamma^(-i) 和條件分布：algo_esf()（algo.R；和 CML 的 cml_gamma_matrix() 同一段程式）
cond_table <- function(key, beta_of) {
  items <- set_items(key)
  esf <- algo_esf(beta_of[items])
  r <- 0:esf$max_score
  data.frame(set = key, raw = rep(r, length(items)), item_id = rep(items, each = length(r)),
             E = as.vector(esf$expected), W = as.vector(esf$variance),
             free = as.vector(esf$free))
}

# plot() 的題目適配畫的可接受範圍，也是 kits_item() 印出和診斷用的範圍（只在 plot.kits_itemfit() 的預設寫一次）
itemfit_range <- function() eval(formals(plot.kits_itemfit)$range)

# 題目的 flag（Rasch 的用語，Linacre；Wright & Linacre 1994）：infit 或 outfit 高於範圍是 underfit（雜訊，
# 比較嚴重，兩個都有時算它）；不然低於範圍是 overfit（太好預測）；其他 ok；剛好在界線上算 ok；兩個都 NA 時 NA
itemfit_flag <- function(infit, outfit, range = itemfit_range()) {
  high <- (infit > range[2]) %in% TRUE | (outfit > range[2]) %in% TRUE
  low <- (infit < range[1]) %in% TRUE | (outfit < range[1]) %in% TRUE
  case_when(is.na(infit) & is.na(outfit) ~ NA_character_,
            high ~ "underfit",
            low ~ "overfit",
            .default = "ok")
}

# 每個題本每一題：人數、通過率（平均分 / 這題觀察到的最高分）、和題本其他題總分的相關（內部用，給 kits_itemfit()）。
# 和 dexter::tia_tables()$items 的 n_persons、pvalue、rir 相同。題本和題照名字排。
item_totals <- function(resp, sets = person_sets(resp)) {
  resp |>
    mutate(max_score = max(item_score), .by = item_id) |>
    left_join(rename(sets, book = set, total = raw), by = "visit_id") |>
    summarise(n_persons = n(), pvalue = mean(item_score) / max_score[1], rir = rest_cor(item_score, total - item_score),
              .by = c(book, item_id)) |>
    arrange(sort_rank(book), sort_rank(item_id))
}

# 題本的值的加權平均（權重是人數）。題本很小時算不出來（例如只有一個人），只用有值的題本；都沒有時 NA。
booklet_mean <- function(v, n) {
  ok <- is.finite(v)
  if (!any(ok)) NA_real_ else sum((v * n)[ok]) / sum(n[ok])
}

# 題目和其他題總分的相關。題目沒有變異、只有一個人、或其他題的總分沒有變異時是 NA（0/0），不計入加權平均。
# （dexter 在其他題總分沒有變異的小題本裡算出來的是 0/0 的浮點雜訊：NA、正負 1 或約 1e-8；這裡一律 NA。）
rest_cor <- function(x, rest) {
  if (length(x) < 2 || stats::sd(x) == 0 || stats::sd(rest) == 0) {
    return(NA_real_)
  }
  stats::cor(x, rest)
}
