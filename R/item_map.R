# 常模中位數第一次到達某個 D 的月齡（內部用，給 kits_item_table()）。
# 沿常模表的月齡往上找，回傳第一個 P50 >= D 的月齡；不內插。
# 有好幾段（segment 欄）時，從月齡最小的一段開始找。找不到時回傳 NA。
reach_age <- function(D, norm) {
  for (seg in norm_segments(norm)) {
    hit <- which(seg$P50 >= D)
    if (length(hit) > 0) {
      return(seg$age_m[hit[1]])
    }
  }
  NA_real_
}


# 題目在成長曲線上的位置（內部用，給 kits_item_table()）：
# 有經驗的通過月齡就用它（empirical），沒有時用常模中位數第一次到達這個 D 的月齡（age_reach）
item_positions <- function(items, norm, rp = 0.5) {
  need_rp(rp)
  age <- age_col(rp)
  need_cols(items, c("domain", "D", age), "items")
  need_cols(norm, c("domain", "age_m", "P50"), "norm")
  miss <- setdiff(unique(items$domain), unique(norm$domain))
  abort_if(length(miss) > 0, "Domain(s) in `items` with no norm in `norm`: ", paste(miss, collapse = ", "), ".")
  norm <- as.data.frame(norm)
  p50_range <- norm |>
    summarise(lowest_p50 = P50[which.min(age_m)], highest_p50 = max(P50), .by = domain)
  reach <- \(d, sc) reach_age(d, filter(norm, domain %in% sc))
  out <- items |>
    as.data.frame() |>
    left_join(p50_range, by = "domain") |>
    mutate(age_reach = map2_dbl(D, domain, reach),
           below_range = D < lowest_p50,
           above_range = D > highest_p50,
           empirical = !is.na(.data[[age]]),
           x = if_else(empirical, .data[[age]], age_reach)) |>
    select(-lowest_p50, -highest_p50) |>
    kits_table()
  rownames(out) <- rownames(items)   # left_join() 不留列名
  if (rp != 0.5) attr(out, "rp") <- rp
  out
}

# 經驗月齡的欄名：rp = .5 是 age50，.8 是 age80
age_col <- function(rp) sprintf("age%g", 100 * rp)


# 階梯右邊的刻度（內部用）：這些月齡上，幾個領域的常模中位數（P50）的平均；每個領域都有常模的月齡才留
ladder_marks <- function(norm, domains, ages) {
  need_cols(norm, c("domain", "age_m", "P50"), "norm")
  domains <- unique(domains)
  miss <- setdiff(domains, norm$domain)
  abort_if(length(miss) > 0, "Domain(s) not in `norm`: ", paste(miss, collapse = ", "), ".")
  out <- norm |>
    as.data.frame() |>
    filter(domain %in% domains) |>
    mutate(age_m = round(age_m, 6)) |>                         # 3 和 3.0000000000000004 當成同一個月齡
    filter(age_m %in% round(ages, 6)) |>
    summarise(P50 = mean(P50), .by = c(domain, age_m)) |>      # 兩段共用的月齡：每個領域先平均
    filter(n() == length(domains), .by = age_m) |>             # 每個領域都有常模的月齡
    arrange(age_m, sort_rank(domain)) |>
    summarise(D = mean(P50), .by = age_m)
  abort_if(nrow(out) == 0, "None of `ages` has a norm in every domain.")
  kits_table(out)
}
