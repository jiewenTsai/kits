# PCM 下，能力 theta 的人在一題得到 k 分以上的機率（內部用，給門檻和它的 SE）。
# theta 可以是向量，theta 和 betas 也可以是 RTMB 的 advector。類別機率來自 pcm_kernel()（pcm_core.R）：P(X >= k) = 第 k..m 欄的和 / 列和
p_at_least <- function(theta, betas, k) {
  abort_if(k < 1 || k > length(betas), "`k` must be between 1 and the number of steps (", length(betas), ").")
  q <- pcm_kernel(theta, betas)
  row_sum(q$ex[, (k + 1):(length(betas) + 1), drop = FALSE]) / q$S
}

# 一題一個步驟的門檻（Thurstonian threshold；內部用）：一個能力值，在這個能力得到 k 分以上的機率
# 剛好是 rp（預設 .5，RP50）。uniroot() 解 P(X >= k | theta) = rp（不用導數）。回傳 c(thr = )
step_threshold <- function(betas, k = 1, rp = 0.5) {
  c(thr = stats::uniroot(function(t) p_at_least(t, betas, k) - rp, c(-40, 40), tol = 1e-10)$root)
}

# 步驟參數加門檻的表（內部用，給 kits_item()、kits_thresholds()）。一列一題的一個步驟。
# fit 是 cml_estimate() 的結果（題照名字排，每題的步驟照順序）；門檻和 SE 是 algo_thresholds()（algo.R）。
# 沒放進校準的固定題（fit$inputs$left_out）：步驟參數照抄，門檻從固定的 beta 算，SE 是 0；再照題號排。
item_thresholds <- function(fit) {
  abort_if(!inherits(fit, "kits_fit"),
           "`fit` must be the result of the CML calibration (kits_cml(), or kits_item()$fit).")
  thr <- algo_thresholds(list(coef = rename(fit$coef, se = SE_beta), vcov = fit$est$acov.beta))
  left_out <- fit$inputs$left_out
  if (is.null(left_out) || nrow(left_out) == 0) {
    return(thr)
  }
  bind_rows(thr, copied_thresholds(left_out, unique(left_out$item_id))) |>
    arrange(sort_rank(item_id), item_score)
}

# 沒放進校準的固定題：步驟參數照抄，門檻從固定的 beta 算（step_thresholds()），SE 是 0
copied_thresholds <- function(fixed, bad) {
  fixed |>
    filter(item_id %in% bad) |>
    arrange(item_id, item_score) |>
    mutate(beta_se = 0, thr = step_thresholds(beta)$thr, thr_se = 0, .by = item_id) |>
    select(item_id, item_score, beta, beta_se, thr, thr_se)
}

# 一題每個步驟的門檻和 SE（V 是這題步驟參數的共變異矩陣；NULL 時 SE 是 NA）。
# SE：delta method sqrt(J V t(J))，J = d thr / d beta 交給 RTMB（threshold_jacobian()）
step_thresholds <- function(betas, V = NULL, rp = 0.5) {
  thr <- map_dbl(seq_along(betas), \(k) step_threshold(betas, k, rp)[["thr"]])
  if (is.null(V)) {
    return(data.frame(thr = thr, thr_se = NA_real_))
  }
  J <- threshold_jacobian(betas, thr, rp)
  data.frame(thr = thr, thr_se = sqrt(pmax(rowSums((J %*% V) * J), 0)))
}

# 門檻對步驟參數的導數（一列一個門檻，一欄一個 beta），不用手推：
# 門檻是 P(X >= k | theta_k) = rp 的解，寫成對 theta_1..theta_m 的最小化 sum_k (P(X >= k | theta_k) - rp)^2；
# RTMB 的 tape$newton() 把它變成「beta -> 門檻」的 tape（從 uniroot() 的解出發），
# 這個 tape 的 jacobian 就是 d thr / d beta（隱函數的微分由 RTMB 做）
threshold_jacobian <- function(betas, thr, rp) {
  m <- length(betas)
  gap <- function(x) {
    theta <- x[seq_len(m)]
    beta <- x[m + seq_len(m)]
    total <- 0
    for (k in seq_len(m)) total <- total + (p_at_least(theta[k], beta, k) - rp)^2
    total
  }
  solve_thr <- RTMB::MakeTape(gap, c(thr, betas))$newton(seq_len(m), grad_tol = 1e-14)
  matrix(solve_thr$jacobian(betas), m, m)
}

# 門檻表的 RP（沒有記錄時是 RP50）
rp_of <- function(x) {
  rp <- attr(x, "rp")
  if (is.null(rp)) 0.5 else rp
}

need_rp <- function(rp, arg = "rp") {
  need_number(rp, arg, "one probability between 0 and 1, for example 0.5 or 0.8", above = 0, below = 1)
}

# 每個步驟在 RP rp 的門檻和 SE（內部用，給 coef(item, rp = )）。beta 和它的共變異來自校準；
# 沒放進校準的固定題（不在 fit 裡）SE 是 0。
item_thresholds_rp <- function(item, rp) {
  V <- item$fit$est$acov.beta
  fit_rows <- item$fit$coef |>
    transmute(item_id, item_score, fit_row = row_number())
  thr <- item$thresholds |>
    left_join(fit_rows, by = c("item_id", "item_score")) |>
    mutate(rp_thresholds(beta, fit_row, V, rp), .by = item_id) |>
    select(-fit_row)
  attr(thr, "rp") <- rp
  thr
}

# 一題在 RP rp 的門檻（步驟照順序）；這題的步驟不全在校準裡時（沒放進校準的固定題）SE 是 0
rp_thresholds <- function(beta, fit_row, V, rp) {
  Vi <- if (anyNA(fit_row) || is.null(V)) NULL else V[fit_row, fit_row, drop = FALSE]
  out <- step_thresholds(beta, Vi, rp)
  if (is.null(Vi)) out$thr_se <- 0
  out
}
