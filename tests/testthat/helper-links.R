# 測試用的小函數。以前的 kits_mean_sigma()、kits_age50() 等（0.22.0 停用，0.23.0 拿掉）的測試改走公開的寫法。

# 自己準備的點的連結：參數和以前的 kits_mean_sigma() 相同，
# 走 kits_link_screen() |> kits_link_line() |> kits_link_error()
link_of_points <- function(from, to, id, item = id, cut = 2, n_boot = 0, seed = 1,
                           from_scale = NULL, to_scale = NULL, spread = "sd", unit = "step",
                           se_from = NULL, se_to = NULL) {
  points <- data.frame(id = id, item = item, from = from, to = to)
  if (!is.null(se_from)) points$se_from <- se_from
  if (!is.null(se_to)) points$se_to <- se_to
  points |>
    kits_link_screen(cut = cut, spread = spread, unit = unit, from_scale = from_scale,
                     to_scale = to_scale) |>
    kits_link_line() |>
    kits_link_error(n_boot = n_boot, seed = seed)
}

# 50% 通過的月齡：內部的 age_pass() 的 p50，用題號當名字（以前的 kits_age50() 的樣子）
age50_of <- function(pass, age_m, wave_m, step = 1) {
  out <- kits:::age_pass(pass, age_m, wave_m, step = step)
  if (anyNA(out$item_id)) {
    return(out$p50)
  }
  stats::setNames(out$p50, out$item_id)
}

# 有月齡的作答（給 kits_item_table(data = )）：題 j 的通過機率 plogis((age - at[j]) / 2)；
# at 是 Inf 的題沒有人通過
pass_resp <- function(at, age) {
  pass <- sapply(at, \(a) stats::rbinom(length(age), 1, stats::plogis((age - a) / 2)))
  data.frame(visit_id = paste0("p", seq_along(age)), age_m = age, wave_m = round(age / 3) * 3, pass)
}

# 選定的 theta 經過一串連結：換算值和連結誤差（以前的 kits_link_se()；現在換算一張 theta 的表）
link_se_of <- function(..., theta, from) {
  out <- kits_convert(data.frame(theta = theta), ..., from = from)
  data.frame(theta = theta, value = out$theta, se_link = out$se_link)
}

# 每一波的原始分換算表（kits 0.22.0 起的寫法）：先計分（每一列的波次放在 persons），再 kits_subj_table()
subj_table_of <- function(item, X, wave_m = NULL) {
  persons <- if (is.null(wave_m)) NULL else data.frame(visit_id = rownames(X), wave_m = wave_m)
  kits_subj(X, item, persons = persons) |>
    kits_subj_table()
}

# 題目表的內部的函數（kits_item_table() 和它的 plot() 用）：0.21.0 的表（領域的欄叫 score）先改名
score_dom <- function(d) {
  rp <- attr(d, "rp")
  d <- as.data.frame(d)
  attr(d, "rp") <- rp
  if ("score" %in% names(d) && !"domain" %in% names(d)) names(d)[names(d) == "score"] <- "domain"
  d
}
old_item_positions <- function(items, norm, rp = 0.5) {
  out <- kits:::item_positions(score_dom(items), score_dom(norm), rp)
  names(out)[names(out) == "domain"] <- "score"
  out
}
old_ladder_marks <- function(norm, scores, ages) kits:::ladder_marks(score_dom(norm), scores, ages)
old_plot_growth <- function(norm, items, scores, age_ticks = AGE_TICKS, y_range = NULL) {
  kits:::growth_plot(score_dom(norm), score_dom(items), scores, age_ticks, y_range)
}
old_plot_ladder <- function(items, columns,
                            marks) kits:::ladder_plot(score_dom(items), columns, marks)

# bind_rows() 疊 kits 的表：只留第一張表的紀錄，會警告（這些測試測的是別的事）
quiet_records <- function(expr) {
  withCallingHandlers(expr, kits_warning_scale = function(w) {
    if (grepl("first table only", conditionMessage(w))) invokeRestart("muffleWarning")
  })
}
