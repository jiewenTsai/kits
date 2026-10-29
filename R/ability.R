#' Person abilities (WLE) that can always be computed
#'
#' Estimates the weighted likelihood estimate (WLE; Warm, 1989) of every person-wave and its SE,
#' finite also for a zero or a perfect score, plus person fit. theta is on the scale of `item`
#' (recorded in `attr(, "scale")`); convert it with [kits_convert()].
#' The result keeps the raw-score tables of its item patterns, so `scores |> kits_subj_table()`
#' gives the table of each wave.
#' For the PCM the WLE equals the mode of the likelihood penalized by the Jeffreys prior (Firth,
#' 1993; Magis, 2015); kits computes it by maximizing log L(theta) + 1/2 log I(theta) with RTMB
#' ([algo_jme()]). With two maxima (few items with very different difficulties) the higher one is
#' taken (the objective value decides); with two equally high maxima, the smaller theta.
#' The SE is 1 / sqrt(I(theta)), the asymptotic SE of the WLE (Warm, 1989).
#'
#' Details: `kits_algorithm("wle")` and the tutorial, section wle;
#' the models tutorial (tutorial_models/kits_models.html), section 4.
#'
#' @param data Responses, as in [kits_item()]: a matrix, a wide or long data.frame, or a
#'   [kits_resp()] table (its columns of the child are then the `persons`). One domain only;
#'   gaps in the observed scores are allowed.
#' @param item The item parameters: the result of [kits_item()] (or its `fit`), or a [kits()],
#'   or a data.frame with `item_id`, `item_score`, `beta`.
#' @param persons Optional person information to add to the result: a data.frame with `visit_id` and
#'   any of `subj_id` (the child id, the same in every wave), `wave_m` (wave, months), `age_m` (age,
#'   months), `sex` and `cohort`; other columns are dropped. Every person-wave needs a row.
#' @return A data.frame with one row per person, in the order of `data`:
#'   * `pattern_id`, `visit_id`: the item pattern (the combination of items answered) and the id.
#'   * `raw`: the raw score.
#'   * `raw_max`: the highest possible raw score on the items answered (`raw == raw_max`: perfect).
#'   * `theta`, `se`: the WLE and its SE, 1 / sqrt(I(theta)).
#'   * `infit`, `outfit`: person fit, the mean squares of this person's residuals (NA for a zero or
#'   a perfect score; much below 1 with few items, so compare persons of the same form).
#'   * `flag`: `"zero"` (zero score), `"perfect"` (perfect score), `"misfit"` (infit or outfit above
#'   2), or `"ok"`.
#'   * With `persons`, then its columns `subj_id`, `wave_m`, `age_m`, `sex`, `cohort` (those
#'   present).
#' @examples
#' betas <- list(a1 = c(-1, 0), a2 = c(0, 1), a3 = 0.5, a4 = -0.5)
#' X <- kits_sim(seq(-3, 3, length.out = 300), betas, seed = 1)
#' item <- kits_item(X)
#' head(kits_subj(X, item))
#'
#' # Add the person information: the norm functions then take the result directly
#' library(dplyr)
#' persons <- tibble(visit_id = rownames(X), wave_m = rep(c(12, 24), each = 150)) |>
#'   mutate(age_m = wave_m + runif(n(), -0.5, 0.5))
#' s <- kits_subj(X, item, persons = persons)
#' head(s)
#' summary(s)
#'
#' # Or the responses and the child's columns in one table (kits_resp)
#' resp <- kitm03$responses |>
#'   filter(domain == "lan", wave_m <= 6) |>
#'   kits_resp()
#' item_resp <- kits_item(resp, scale = "lan M03")
#' head(kits_subj(resp, item_resp))
#'
#' # dplyr verbs keep the scale of the scores
#' s |>
#'   filter(flag == "ok") |>
#'   summarise(mean = mean(theta), .by = wave_m)
#' @export
kits_subj <- function(data, item, persons = NULL, ...) {
  no_extra_args(..., fun = "kits_subj")
  call <- call_text(sys.call())
  data_expr <- substitute(data)
  if (is.null(persons)) persons <- resp_persons(data)   # 作答表：孩子的欄跟著走
  resp <- resp_long(data, gaps = FALSE)
  parms <- param_frame(item)
  need_items(resp, parms)

  sets <- person_sets(resp, parms) |>
    add_wle(parms) |>
    add_booklet_ids()
  est <- sets |>
    add_person_fit(resp, parms) |>
    join_persons(persons)   # 人的資料接在 flag 後面
  est <- kits_table(est, "kits_subj")
  attr(est, "scale") <- beta_scale_of(item)
  attr(est, "raw_theta") <- list(raw_theta_source(parms, beta_scale_of(item), sets))
  attr(est, "diagnostics") <- subj_diagnostics(est)
  attr(est, "info") <- make_info(call, list(scale = beta_scale_of(item)), data_fingerprint(resp), item_data = source_of(item),
                                 persons = if (!is.null(persons)) {
                                   persons_fingerprint(est[c("visit_id", person_cols(est))])
                                 },
                                 input = data_input(data, data_expr))
  est
}

# WLE 和 SE 只和「題目組合 × 原始分」有關，所以每個組合只算一次，再接回每個人次
add_wle <- function(est, parms) {
  combos <- distinct(est, set, raw)
  combos <- bind_cols(combos, algo_jme(betas_by_item(parms), combos$raw, map(combos$set, set_items)))
  left_join(est, combos, by = c("set", "raw"))
}

# 題本的名字和 dexter 相同（number_booklets()）
add_booklet_ids <- function(est) {
  books <- select(number_booklets(est), set, pattern_id)
  left_join(est, books, by = "set")
}

# 個人適配（infit、outfit）和 flag；0 分和滿分的人沒有適配，給 NA
add_person_fit <- function(est, resp, parms) {
  fit <- pcm_moments(resp, parms, select(est, visit_id, theta)) |>
    mean_squares(visit_id)
  est |>
    left_join(select(fit, visit_id, infit, outfit), by = "visit_id") |>
    mutate(inner = raw > 0 & raw < raw_max,
           infit = if_else(inner, infit, NA_real_),
           outfit = if_else(inner, outfit, NA_real_),
           flag = case_when(raw == 0 ~ "zero",
                            raw == raw_max ~ "perfect",
                            .default = if_else(pmax(infit, outfit) > 2, "misfit", "ok"))) |>
    select(pattern_id, visit_id, raw, raw_max, theta, se, infit, outfit, flag)
}

# 分數表記下的原始分 -> theta 的來源（給 kits_subj_table(scores)）：題目參數（彙總的數字）、量尺、
# 每個題本是哪些題、這個校準的人次，和之後的換算（kits_convert() 一次一步）。沒有孩子的作答
raw_theta_source <- function(parms, scale, sets) {
  list(parms = parms, scale = scale, booklets = distinct(as.data.frame(sets), pattern_id, set),
       persons = unique(as.character(sets$visit_id)), steps = list())
}

# kits_subj() 可以接上的人的欄
person_cols <- function(d) intersect(c("wave_m", "age_m", "sex", "cohort"), names(d))

# 分數表留下的孩子的欄：subj_id（孩子的編號）和上面的欄
KEEP_PERSON_COLS <- c("subj_id", "wave_m", "age_m", "sex", "cohort")

# 作答表（kits_resp）的孩子的欄，當成 persons；不是作答表，或沒有這些欄時 NULL
resp_persons <- function(data) {
  if (!inherits(data, "kits_resp")) {
    return(NULL)
  }
  cols <- intersect(KEEP_PERSON_COLS, names(data))
  if (length(cols) == 0) {
    return(NULL)
  }
  as.data.frame(data)[c("visit_id", cols)]
}

# persons 的欄（wave_m、age_m、sex、cohort）照 visit_id 接上
join_persons <- function(est, persons) {
  if (is.null(persons)) {
    return(est)
  }
  persons <- old_ids(persons, "persons")
  if (!is.data.frame(persons)) {
    kits_abort("`persons` must be a data.frame with visit_id and any of wave_m, age_m, sex, ",
               "cohort.",
               class = "kits_error_input")
  }
  need_cols(persons, "visit_id", "persons")
  cols <- intersect(KEEP_PERSON_COLS, names(persons))
  unused <- setdiff(names(persons), c("visit_id", cols))
  if (length(cols) == 0) {
    kits_abort("`persons` has none of the columns subj_id, wave_m, age_m, sex, cohort (it has: ",
               paste(names(persons), collapse = ", "),
               "). Rename them, e.g. dplyr::rename(persons, wave_m = wave, age_m = age).", class = "kits_error_input")
  }
  if (length(unused) > 0) {
    kits_inform("Column(s) of `persons` not used: ", paste(unused, collapse = ", "),
                ". kits_subj() adds only subj_id (child id), wave_m (wave, months), ",
                "age_m (age, months), sex and cohort; rename a column to use it.", class = "kits_message_data")
  }
  persons <- persons |>
    as.data.frame() |>
    mutate(visit_id = as.character(visit_id)) |>
    select(visit_id, all_of(cols))
  ids <- persons$visit_id
  if (anyNA(ids)) {
    kits_abort("`persons$visit_id` has ", sum(is.na(ids)),
               " missing id(s); give one row per person-wave.",
               class = "kits_error_data")
  }
  if (anyDuplicated(ids)) {
    dup <- unique(ids[duplicated(ids)])
    kits_abort("`persons$visit_id` has ", length(dup), " repeated id(s)", ids_text(dup),
               "; give one row per person-wave.",
               class = "kits_error_data")
  }
  for (v in intersect(c("wave_m", "age_m"), cols)) {
    if (!is.numeric(persons[[v]])) {
      kits_abort("`persons$", v, "` must be numbers (months).", class = "kits_error_data")
    }
  }
  miss <- setdiff(est$visit_id, ids)
  if (length(miss) > 0) {
    kits_abort(length(miss), " person-wave(s) of the data have no row in `persons`", ids_text(miss, 3), ".",
               class = "kits_error_data")
  }
  left_join(est, persons, by = "visit_id")
}

# kits_subj() 的診斷：0 分、滿分、不適配的比例；知道波次時，每個波次的信度
subj_diagnostics <- function(est) {
  pct <- function(f) 100 * mean(est$flag == f)
  share <- function(f, cut, what) {
    v <- pct(f)
    above <- v > cut
    status <- if (above) "note" else "ok"
    detail <- if (above) sprintf(" (above %g%%)", cut) else ""
    diag_row(paste("share of", what), v, status, sprintf("%.1f%% of person-waves%s", v, detail))
  }
  rows <- list(share("zero", 10, "zero scores"), share("perfect", 10, "perfect scores"), share("misfit", 5, "misfit"))
  if ("wave_m" %in% names(est)) {
    waves <- est |>
      filter(!is.na(wave_m)) |>
      summarise(n = n(), var_theta = stats::var(theta), mean_se2 = mean(se^2), .by = wave_m) |>
      filter(n >= 3, var_theta != 0) |>
      arrange(wave_m) |>
      mutate(rel = 1 - mean_se2 / var_theta)
    rows <- c(rows, map2(waves$wave_m, waves$rel, reliability_row))
  }
  diag_table(rows)
}

# 一個波次的信度：1 - mean(se^2) / var(theta)；低於 0 警告、低於 0.5 註記
reliability_row <- function(w, rel) {
  status <- case_when(rel < 0 ~ "warn", rel < 0.5 ~ "note", .default = "ok")
  detail <- case_when(
    rel < 0 ~ paste(": below 0, the wave cannot separate children. Usually the items are too easy (or too hard)",
                    "for this wave, so many children have a perfect (or zero) score, whose se is large;",
                    "see the shares of perfect and zero scores (flag) and kits_precision()"),
    rel < 0.5 ~ ": below 0.5; check the shares of perfect and zero scores (flag)",
    .default = "")
  diag_row(paste("reliability, wave", format(w)), rel, status, paste0("1 - mean(se^2) / var(theta)", detail))
}

#' @rdname kits_subj
#' @param x A `kits_subj`.
#' @param digits Decimals to print (only the printout is rounded; the object keeps full precision).
#' @param n Rows to print: by default all rows of a table of up to 30 rows, else the first 10.
#' @param ... In `kits_subj()`: not used (an unknown argument stops with an error).
#'   In `print()`: passed to `print()` for a data.frame.
#' @export
print.kits_subj <- function(x, digits = 2, n = NULL, ...) {
  need_digits(digits)
  need_rows(n)
  cat(sprintf("<kits_subj> %d person-waves; scale: %s\n", nrow(x), scale_label(scale_of(x))))
  if (!whole_subj(x)) {   # 使用者選過或排過欄：照使用者的欄和順序印，不加說明的行
    print_rows(as.data.frame(x), digits, n, ...)
    return(invisible(x))
  }
  main <- intersect(c("visit_id", "wave_m", "age_m", "theta", "se", "z", "pct", "flag"), names(x))
  print_rows(as.data.frame(x)[main], digits, n, ...)
  others <- setdiff(names(x), main)
  if (length(others) > 0) cat(sprintf("Other columns: %s\n", paste(others, collapse = ", ")))
  if (!is.null(attr(x, "raw_theta"))) {
    cat("Also kept: the raw score -> theta of each item pattern ",
        "(kits_subj_table(x) gives the table of each wave).\n", sep = "")
  }
  cat("Next: kits_convert(x, link), kits_norms(x), kits_norm_z(x, norm), kits_subj_table(x), ",
      "kits_vertical(x).\n", sep = "")
  checks_line(attr(x, "diagnostics"))
  invisible(x)
}


# kits_subj() 的整張表：計分的欄都在，順序也沒變（select() 或 relocate() 之後不是）
whole_subj <- function(x) {
  core <- c("pattern_id", "visit_id", "raw", "raw_max", "theta", "se", "infit",
            "outfit", "flag")
  all(core %in% names(x)) && identical(intersect(names(x), core), core)
}

# 每個題本、每個原始分（0 到滿分）的 WLE 表（內部用；kits_subj_table(item)、kits_precision() 用）。
# 和 kitsPreview::kits_wle_table()、dexter::ability_tables(method = "WLE") 相同。
# design：題本（pattern_id、item_id）；NULL 時 kits_item 用它的題本，data.frame 的所有題是一個題本。
# call：記在 kits_info() 的呼叫
wle_table <- function(item, design = NULL, call = NULL) {
  parms <- item
  if (inherits(parms, "kits")) parms <- parms$item
  if (inherits(parms, "kits_item")) {
    if (is.null(design)) design <- parms$design   # 包括沒放進校準的固定題
    parms <- parms$thresholds
  }
  if (is.null(design)) {
    all_items <- \() data.frame(pattern_id = "all", item_id = unique(param_frame(parms)$item_id))
    fit_design <- inherits(parms, "kits_fit")
    design <- if (fit_design) parms$inputs$design[, c("pattern_id", "item_id")] else all_items()
  }
  design <- as.data.frame(design)
  need_cols(design, c("pattern_id", "item_id"), "design")
  parms <- param_frame(parms)
  need_params(design$item_id, parms$item_id, "Item(s) in `design` have no parameters in `parms`: ")
  need_items(data.frame(visit_id = "p", item_id = parms$item_id, item_score = 0L), parms)
  item_max <- summarise(parms, item_max = max(item_score), .by = item_id)

  # 題本照名字排序（和 dexter 相同），每個題本 0 到滿分各一列
  rows <- design |>
    transmute(pattern_id = as.character(pattern_id), item_id = as.character(item_id)) |>
    distinct() |>
    left_join(item_max, by = "item_id") |>
    summarise(items = list(sort(item_id)), max_raw = sum(item_max), .by = pattern_id) |>
    arrange(sort_rank(pattern_id)) |>
    reframe(items = rep(items, max_raw + 1), raw = 0:max_raw, .by = pattern_id)
  fill <- algo_jme(betas_by_item(parms), rows$raw, rows$items)
  tab <- kits_table(data.frame(pattern_id = rows$pattern_id, raw = rows$raw, theta = fill$theta,
                               se = fill$se), "kits_pattern_table")
  attr(tab, "scale") <- beta_scale_of(item)
  attr(tab, "info") <- make_info(call, list(scale = beta_scale_of(item)), item_data = source_of(item))
  tab
}


# kits_item、它的 fit、或 data.frame → item_id、item_score、beta
param_frame <- function(parms) {
  if (inherits(parms, "kits")) parms <- parms$item
  if (inherits(parms, "kits_item")) parms <- parms$thresholds   # 也包括沒放進校準的固定題
  if (inherits(parms, "kits_fit")) {
    need_cml(parms)
    parms <- stats::coef(parms)
  }
  if (!is.data.frame(parms)) {
    kits_abort("The item parameters must be the result of kits_item(), or a data.frame with item_id, item_score, beta.",
               class = "kits_error_input")
  }
  need_cols(parms, c("item_id", "item_score", "beta"), "parms")
  parms <- data.frame(item_id = as.character(parms$item_id),
                      item_score = suppressWarnings(as.integer(as.character(parms$item_score))),
                      beta = parms$beta)
  if (!is.numeric(parms$beta) || any(!is.finite(parms$beta))) {
    kits_abort("`parms$beta` must be finite numbers (no NA or Inf).", class = "kits_error_items")
  }
  if (any(is.na(parms$item_score) | parms$item_score < 1)) {
    kits_abort("`parms$item_score` must be whole numbers 1, 2, ... (one row per step).", class = "kits_error_items")
  }
  if (anyDuplicated(parms[c("item_id", "item_score")])) {
    kits_abort("`parms` has repeated item steps.", class = "kits_error_items")
  }
  parms
}

# 只支援 CML 的校準結果
need_cml <- function(fit) {
  if (!identical(fit$inputs$method, "CML")) {
    kits_abort("Only CML calibrations are supported (the default method = \"CML\" of kits_item()).", class = "kits_error_input")
  }
}

# resp 的題目都要有參數
need_items <- function(resp, parms) {
  need_params(unique(as.character(resp$item_id)), as.character(parms$item_id),
              "Item(s) in the responses have no parameters in `parms`: ", ". Include only the items of this domain.")
  steps <- parms |>
    summarise(item_max = max(item_score), no_gap = identical(sort(as.integer(item_score)), seq_len(max(item_score))),
              .by = item_id)
  high <- resp |>
    distinct(item_id, item_score) |>
    left_join(steps, by = "item_id") |>
    filter(item_score > item_max)
  if (nrow(high) > 0) {
    bad <- paste0(high$item_id, ":", high$item_score)
    kits_abort("Score(s) above the highest step in `parms`: ", paste(utils::head(bad, 10), collapse = ", "),
               ". Recode the data, or calibrate again.", class = "kits_error_items")
  }
  gap <- steps |>
    filter(!no_gap) |>
    arrange(sort_rank(item_id)) |>
    pull(item_id)
  if (length(gap) > 0) {
    kits_abort("Item(s) in `parms` with missing steps (steps must be 1, 2, ...): ", paste(utils::head(gap, 10), collapse = ", "),
               ".", class = "kits_error_items")
  }
}


# 每個「題目組合 × 原始分」的 WLE（Jeffreys 眾數估計，PCM 時兩者相同）和 SE：algo_jme()（algo.R；給 kits_subj()、wle_table()、kits_link_check()）。
# 下面是它用的步驟：
#   只寫模型，導數交給 RTMB：f = log L(theta) + 1/2 log I(theta)（Warm, 1989），
#   log S_i 和 Var_i 來自 PCM 的核心 algo_pcm()。SE = 1 / sqrt(I(theta))：WLE 的漸近 SE（Warm, 1989）。
#   起始值和 TAM 的 tam.wle() 相同：這組題的 beta 的平均，加上 log((r + 0.5) / (滿分 - r + 1))。
#   加權概似有兩個高點時（題目少、難度差很遠），取比較高的那個（目標函數 f 的值決定）：
#   在格點上算 f（wle_grid()），格點上有比解高出 1e-8 以上的點時，從格點最高的點再解一次。
#   牛頓法從剛好在低點（兩個高點中間）的起始值出發時，梯度是 0，會停在低點；這時也是格點比較高，重解。
#   一樣高的兩個高點（對稱的題目），格點不會比較高，這個起始值偏向低的一邊，所以取 theta 小的那個
#   （和 dexter、0.23.0 相同）。

# 組合 idx 的 WLE：-f 加起來做成 RTMB 的 tape，tape$newton() 從 start 一次解出來
# （RTMB 解不出來時回傳 NaN；一個組合失敗時整批都是 NaN）
wle_solve <- function(setup, raw, idx, start = setup$theta0[idx]) {
  sub <- wle_subset(setup, idx)
  neg_f <- \(d) -sum(wle_f(start + d, sub, raw[idx]))   # d = theta - start
  solve_f <- RTMB::MakeTape(neg_f, numeric(length(idx)))$newton(seq_along(idx), grad_tol = 1e-12)
  start + suppressWarnings(solve_f(numeric(0)))
}

# 每個組合在格點上 f 最高的點（theta）和它的 f（只算模型的值）。格點是所有組合共用的：從最小的步驟參數 - 2
# 到最大的 + 2，每 0.05 一點；每題在格點上的 log S_i、Var_i 算一次，再用「組合 × 題」的 0/1 矩陣加起來
wle_grid <- function(setup, raw) {
  n <- length(raw)
  theta <- seq(min(unlist(setup$beta_of)) - 2, max(unlist(setup$beta_of)) + 2, by = 0.05)
  m <- map(setup$beta_of, \(b) algo_pcm(theta, b))
  A <- matrix(vapply(setup$rows_of, \(rows) seq_len(n) %in% rows, logical(n)), n)   # 組合 j 答第 i 題嗎
  sum_items <- \(v) A %*% do.call(rbind, map(m, v))   # 組合 × 格點：答的題加起來
  f <- outer(raw, theta) - sum_items("log_S") + 0.5 * log(sum_items("Var"))
  f[is.na(f)] <- -Inf
  best <- max.col(f, ties.method = "first")
  list(theta = theta[best], f = f[cbind(seq_len(n), best)])
}

# 每個組合在 theta 的 f = r theta - sum_i log S_i + 1/2 log I（theta 可以是數值或 RTMB 的 advector）
wle_f <- function(theta, setup, raw) {
  parts <- wle_parts(theta, setup)
  raw * theta - parts$log_S + 0.5 * log(parts$info)
}

# 只留組合 idx 的 setup
wle_subset <- function(setup, idx) {
  rows_of <- lapply(setup$rows_of, \(rows) match(rows[rows %in% idx], idx))
  list(beta_of = setup$beta_of, rows_of = rows_of[lengths(rows_of) > 0])
}

# 每個組合在 theta 的 sum_i log S_i 和 I = sum_i Var_i：一題一題加到答這題的組合上
# （theta 可以是數值或 RTMB 的 advector）
wle_parts <- function(theta, setup) {
  "[<-" <- RTMB::ADoverload("[<-")   # 套件裡要明說：advector 放進數值向量
  log_S <- numeric(length(theta))
  info <- numeric(length(theta))
  for (i in names(setup$rows_of)) {
    rows <- setup$rows_of[[i]]
    m <- algo_pcm(theta[rows], setup$beta_of[[i]])
    log_S[rows] <- log_S[rows] + m$log_S
    info[rows] <- info[rows] + m$Var
  }
  list(log_S = log_S, info = info)
}

# 每題的步驟參數（beta_of：一題一個向量，有名字）、答這題的組合，和起始值 theta0
wle_setup <- function(beta_of, sets, raw) {
  max_raw <- map_dbl(sets, \(s) sum(lengths(beta_of[s])))
  mean_beta <- map_dbl(sets, \(s) mean(unlist(beta_of[s])))
  items <- sort(unique(unlist(sets)))
  rows_of <- map(items, \(i) which(map_lgl(sets, \(s) i %in% s)))
  names(rows_of) <- items
  list(beta_of = beta_of[items], rows_of = rows_of,
       theta0 = mean_beta + log((raw + 0.5) / (max_raw - raw + 1)))
}
