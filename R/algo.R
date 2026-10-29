# 演算法層（algo_）：數字進、數字出，不記量尺、少檢查；kits_ 的函數管輸入、量尺的標籤、診斷和印出，
# 計算交給這裡（像 torch 的使用者函數和 kernel 分開）。每個演算法只寫一次：kits_ 的函數呼叫這些函數。
#   algo_pcm()        PCM 的類別機率和動差（kits_pcm_prob()、WLE、適配、tcc()）
#   algo_esf()        CML 的 gamma 和 P(X_i = x | r)（kits_esf()、條件適配、題目–其餘分數）
#   algo_cml()        CML（kits_cml()、kits_item()）       -> algo_thresholds()：門檻和 SE
#   algo_jme()        Jeffreys 眾數估計（PCM 時就是 Warm 的 WLE）和 SE（kits_subj()、kits_precision()）
#   algo_screen()     篩一次（kits_link_screen()）       -> algo_mean_sigma()：直線（kits_link_line()）
#                                                        -> algo_link_error()：重抽（kits_link_error()）
#   algo_line()       直線的套用、反函數、合成（kits_convert() 的換算摘要、kits_precision()）
#   algo_shash()      shash 常模的最大概似（kits_norm_fit()） -> algo_shash_z()：同齡 z（常模的診斷）
#   algo_lw()         Lord-Wingersky 的原始分分布（kits_precision()、觀察分數等化）

#' PCM category probabilities and moments of one item
#'
#' The partial credit model for one item at abilities `theta`: the category probabilities, the
#' expected score and the variance (the item information). It does not track scales or check its
#' input; [kits_pcm_prob()] does.
#'
#' @param theta Abilities (a numeric vector; also an RTMB advector).
#' @param beta The step parameters of the item, in step order.
#' @param third Also return the third central moment `M3`.
#' @return A list: `P` (a matrix, one row per theta, one column per score 0..m), `E`, `Var`,
#'   `log_S` (the log of the normalizing sum), and `M3` when `third = TRUE`.
#' @references
#' Masters, G. N. (1982). A Rasch model for partial credit scoring. *Psychometrika, 47*(2),
#' 149-174. \doi{10.1007/BF02296272}
#' @family algo
#' @examples
#' p <- c(-1, 0, 1) |>
#'   algo_pcm(beta = c(-0.5, 0.5))
#' p$P
#' p$E
#' @export
algo_pcm <- function(theta, beta, third = FALSE) {
  q <- pcm_kernel(theta, beta)
  P <- q$ex / q$S
  k <- seq_len(ncol(P)) - 1
  E <- drop(P %*% k)
  dev <- outer(-E, k, "+")                  # k - E
  Var <- row_sum(P * dev^2)                 # 中央動差：極端的 theta 不會相消
  if (!is_ad(Var)) Var <- pmax(Var, 0)      # 數值的捨入不讓它小於 0（RTMB 不能比大小，也不需要）
  out <- list(P = P, E = E, Var = Var, log_S = q$log_S)
  if (third) out$M3 <- row_sum(P * dev^3)
  out
}

#' Elementary symmetric functions and the conditional item-score distributions
#'
#' For one set of items, the elementary symmetric functions gamma_r of CML and, for each item, the
#' distribution of its score given the raw score r, with its mean and variance. It does not track
#' scales or check its input; [kits_esf()] does.
#'
#' @param betas The step parameters: a named list, one numeric vector per item (in step order).
#' @return A list: `items`, `m` (steps per item), `max_score`, `centre`, `log_gamma` (r = 0..R),
#'   `log_gamma_rest` (one row per item: gamma without that item), `prob` (one matrix per item:
#'   P(X_i = x | r), rows r, columns x), `expected`, `variance` and `free` (rows r, columns items).
#' @family algo
#' @examples
#' esf <- list(a1 = c(-1, 0.5), a2 = 0.3) |>
#'   algo_esf()
#' exp(esf$log_gamma)
#' esf$prob$a1
#' @export
algo_esf <- function(betas) {
  # gamma 用 CML 的 cml_book_design() 和 cml_gamma_matrix()（cml.R）算：第 1 個「題本」是整個組合，
  # 第 1 + i 個是拿掉第 i 題的組合；所有列用同一個置中 c（整個組合的步驟平均），所以可以互相除
  items <- names(betas)
  n <- length(items)
  m <- unname(lengths(betas))
  beta <- unlist(betas, use.names = FALSE)
  inp <- list(items = items, m = m, step_item = rep(items, m),
              book_items = c(list(seq_len(n)), map(seq_len(n), \(i) seq_len(n)[-i])))
  design <- cml_book_design(inp)
  c_all <- cml_centre(beta, design)[1]   # CML 的 c_b（整個組合）：和 cml_nll() 同一個函數
  centre <- rep(c_all, n + 1)
  # 第 1 列 gamma(beta - c)；第 1 + i 列 gamma^(-i)(beta - c)
  G <- cml_gamma_matrix(beta, centre, inp, design)
  top <- sum(m)
  r <- 0:top
  # CML 的 log gamma_br：gamma_r(beta) = exp(-r c) gamma_r(beta - c)
  log_gamma <- log(G[1, ]) - r * c_all
  log_rest <- log(G[-1, , drop = FALSE]) - rep(r, each = n) * c_all
  rownames(log_rest) <- items
  prob <- map(seq_len(n), \(j) esf_item_prob(betas[[j]], c_all, G[1 + j, ], G[1, ], m[j], top))
  names(prob) <- items
  expected <- vapply(prob, \(P) drop(P %*% (seq_len(ncol(P)) - 1)), numeric(top + 1))
  variance <- vapply(seq_len(n), \(j) {
    P <- prob[[j]]
    rowSums(P * outer(expected[, j], seq_len(ncol(P)) - 1, \(e, x) (x - e)^2))
  }, numeric(top + 1))
  free <- vapply(m, \(mj) pmax(0, r - (top - mj)) < pmin(mj, r), logical(top + 1))
  dim(expected) <- dim(variance) <- dim(free) <- c(top + 1, n)
  colnames(expected) <- colnames(variance) <- colnames(free) <- items
  list(items = items, m = m, max_score = top, centre = c_all, log_gamma = log_gamma,
       log_gamma_rest = log_rest, prob = prob, expected = expected, variance = variance,
       free = free)
}

#' Conditional maximum likelihood of the PCM
#'
#' Maximizes the conditional likelihood of the partial credit model from the sufficient
#' statistics, with the gradient and Hessian from RTMB. It does not track scales or check its
#' input; [kits_cml()] and [kits_item()] do.
#'
#' @param stats The sufficient statistics: a [kits_stats()] result (or `item$stats`).
#' @return A list: `coef` (one row per step: `item_id`, `item_score`, `beta`, `se`), `betas` (the
#'   betas as a named list, one vector per item), `vcov` (the covariance of the betas), `nll`
#'   (-log L) and `optim` (the convergence record).
#' @family algo
#' @examples
#' betas <- list(a1 = c(-1, 0), a2 = c(0, 1), a3 = 0.5, a4 = -0.5)
#' X <- kits_sim(seq(-3, 3, length.out = 400), betas, seed = 1)
#' cml <- kits_stats(X) |>
#'   algo_cml()
#' cml$coef
#' @export
algo_cml <- function(stats) {
  inp <- cml_inputs(stats)
  con <- cml_constraint(inp, stats$fixed)   # beta = beta0 + A b
  nll <- cml_nll(inp, con)                  # 條件對數概似（RTMB）
  est <- cml_optimize(nll, con)             # beta、V、-log L、收斂的紀錄
  coef <- data.frame(item_id = inp$step_item, item_score = inp$step_k, beta = est$beta,
                     se = sqrt(pmax(diag(est$V), 0)))
  list(coef = coef, betas = betas_by_item(coef), vcov = est$V, nll = est$value, optim = est$optim)
}

#' Thurstonian thresholds of the steps and their SEs
#'
#' For each step k of each item, the ability where P(X >= k | theta) = `rp` (by `uniroot()`), and
#' its SE by the delta method with the derivative from RTMB. It does not track scales or check its
#' input; [kits_thresholds()] and [kits_item()] do.
#'
#' @param cml An [algo_cml()] result (its `coef` and `vcov`).
#' @param rp The response probability (0.5: RP50).
#' @return A data.frame, one row per step: `item_id`, `item_score`, `beta`, `beta_se`, `thr`,
#'   `thr_se`.
#' @family algo
#' @examples
#' betas <- list(a1 = c(-1, 0), a2 = c(0, 1), a3 = 0.5, a4 = -0.5)
#' X <- kits_sim(seq(-3, 3, length.out = 400), betas, seed = 1)
#' kits_stats(X) |>
#'   algo_cml() |>
#'   algo_thresholds()
#' @export
algo_thresholds <- function(cml, rp = 0.5) {
  V <- cml$vcov
  cml$coef |>
    mutate(fit_row = row_number()) |>
    reframe(item_score = as.integer(item_score), beta, beta_se = se,
            step_thresholds(beta, V[fit_row, fit_row, drop = FALSE], rp), .by = item_id)
}

#' Jeffreys modal estimates of ability and their SEs
#'
#' Jeffreys modal estimate: theta-hat = argmax log L(theta) + 1/2 log I(theta). For the PCM it
#' equals Warm's (1989) WLE exactly (Firth, 1993; Magis, 2015), so [kits_subj()] reports it as the
#' WLE. One estimate per combination of items and raw score, by RTMB's Newton method; finite also
#' for a zero or a perfect score. SE = 1 / sqrt(I(theta-hat)), the asymptotic SE of the WLE
#' (Warm, 1989). With two maxima (few items with very different difficulties) the higher one is
#' taken (the value of the objective decides; it is found on a grid of theta, step 0.05); with two
#' equally high maxima, the smaller theta. It does not track scales or check its input;
#' [kits_subj()] and [kits_subj_table()] do.
#'
#' @param betas The step parameters: a named list, one numeric vector per item.
#' @param raw The raw scores, one per row (default: 0 to the top score of all items).
#' @param items The items answered, a list with one character vector per row (default: all items
#'   for every row).
#' @return A data.frame with `theta` and `se`, one row per raw score.
#' @references
#' Firth, D. (1993). Bias reduction of maximum likelihood estimates. *Biometrika, 80*(1), 27-38.
#' \doi{10.1093/biomet/80.1.27}
#'
#' Magis, D. (2015). A note on weighted likelihood and Jeffreys modal estimation of proficiency
#' levels in polytomous item response models. *Psychometrika, 80*(1), 200-204.
#' \doi{10.1007/s11336-013-9378-5}
#'
#' Warm, T. A. (1989). Weighted likelihood estimation of ability in item response theory.
#' *Psychometrika, 54*(3), 427-450. \doi{10.1007/BF02294627}
#' @family algo
#' @examples
#' betas <- list(a1 = c(-1, 0), a2 = c(0, 1), a3 = 0.5)
#' betas |>
#'   algo_jme()
#' algo_jme(betas, raw = c(1, 1), items = list(c("a1", "a2"), "a2"))
#' @export
algo_jme <- function(betas, raw = 0:sum(lengths(betas)),
                     items = rep(list(names(betas)), length(raw))) {
  # JME：theta = argmax f(theta)，f = log L(theta) + 1/2 log I(theta)（PCM 時等於 Warm 的 WLE）；原始分是充分統計量，
  # log L = r theta - sum_i log S_i(theta)，I = sum_i Var_i(theta)（algo_pcm()）。每列一個 theta，
  # 所有列的 -f 加起來是一個 RTMB tape，tape$newton() 一次解出來（wle_solve()，ability.R）。
  # f 有兩個高點時取高的（f 的值決定）；一樣高時取 theta 小的。
  # SE = 1 / sqrt(I(theta))：WLE 的漸近 SE（Warm, 1989）
  setup <- wle_setup(betas, items, raw)
  n <- length(raw)
  if (n == 0) {
    return(data.frame(theta = numeric(0), se = numeric(0)))
  }
  theta <- wle_solve(setup, raw, seq_len(n))
  grid <- wle_grid(setup, raw)   # 格點上的 f 高出 1e-8 以上（或解不出來）的組合：從格點再解
  up <- which(!(wle_f(theta, setup, raw) >= grid$f - 1e-8))
  if (length(up) > 0) theta[up] <- wle_solve(setup, raw, up, grid$theta[up])
  if (any(!is.finite(theta))) {
    kits_abort("The WLE did not converge for ", sum(!is.finite(theta)), " item pattern x raw score combination(s).",
               class = "kits_error_fit")
  }
  data.frame(theta = theta, se = 1 / sqrt(wle_parts(theta, setup)$info))   # SE = 1 / sqrt(I(theta))
}

#' The mean/sigma line
#'
#' The line new value = intercept + slope x old value with slope = SD(to) / SD(from) and the means
#' aligned. It does not track scales or check its input; [kits_link_line()] does.
#'
#' @param from,to The points on the two scales (numeric vectors of the same length). `from` can
#'   also be an [algo_screen()] result: then the kept points are used and `to` is not needed.
#' @return A list with `intercept` and `slope`.
#' @references
#' Kolen, M. J., & Brennan, R. L. (2014). *Test equating, scaling, and linking* (3rd ed.). Springer.
#' \doi{10.1007/978-1-4939-0317-7}
#' @family algo
#' @examples
#' from <- c(-2, -1, 0, 1, 2, 3)
#' to <- c(-1, 0, 1, 2, 3.2, 6)
#' algo_mean_sigma(from, to)
#' algo_screen(from, to) |>
#'   algo_mean_sigma()
#' @export
algo_mean_sigma <- function(from, to) {
  if (is.list(from)) {   # algo_screen() 的結果：只用留下的點
    kept <- !from$dropped
    to <- from$to[kept]
    from <- from$from[kept]
  }
  mean_sigma_moments(mean(from), stats::sd(from), mean(to), stats::sd(to))
}

#' One round of screening of the link points
#'
#' Fits the mean/sigma line to all points and flags the points whose residual exceeds `cut` times
#' the spread; with `unit = "item"` every point of a flagged item is dropped. It does not track
#' scales or check much; [kits_link_screen()] does.
#'
#' @param from,to The points on the two scales (numeric vectors of the same length).
#' @param item The item of each point (default: one item per point).
#' @param cut The cutoff, in units of the spread.
#' @param spread `"sd"` (the SD of the residuals), `"mad"` (1.4826 x their MAD) or `"se"` (the SE
#'   of each point from `se_from` and `se_to`).
#' @param unit `"step"` (drop the flagged points) or `"item"` (drop every point of a flagged item).
#' @param se_from,se_to The SE of each point (only for `spread = "se"`; `se_to` defaults to 0).
#' @return A list: `from`, `to`, `item`, `resid`, `z`, `flagged`, `dropped` (one value per point),
#'   `intercept0`, `slope0` (the line of all points), `cut` and `rule` (`spread`, `unit`, and the
#'   SEs). Give it to [algo_mean_sigma()] and [algo_link_error()].
#' @family algo
#' @examples
#' from <- c(-2, -1, 0, 1, 2, 3, 0.5, 1.5)
#' to <- c(-1, 0, 1, 2, 3, 4, 1.5, 5.5)
#' s <- algo_screen(from, to, spread = "mad")
#' s$dropped
#' s |>
#'   algo_mean_sigma()
#' @export
algo_screen <- function(from, to, item = seq_along(from), cut = 2, spread = "sd", unit = "step",
                        se_from = NULL, se_to = NULL) {
  # 會失敗的檢查都在這裡，重抽時失敗的就不算（algo_link_error()）
  rule <- screen_rule(spread, unit, cut, se_from, se_to, length(from))
  abort_if(length(from) < 3, "At least 3 points are needed; there are ", length(from), ".")
  abort_if(!is.finite(stats::sd(from)) || !is.finite(stats::sd(to)),
           "`from` and `to` are too large to compute their SD.")
  abort_if(stats::sd(from) == 0 || stats::sd(to) == 0,
           "All values of `from` (or of `to`) are the same, so no line can be fitted.")
  first <- algo_mean_sigma(from, to)
  resid <- line_resid(first, from, to)
  scale <- screen_scale(resid, rule, first$slope)
  z <- resid / scale
  flagged <- abs(resid) > cut * scale
  flagged[is.na(flagged)] <- FALSE
  if (stats::sd(resid) < 1e-8) {   # 點幾乎都在一條線上：不篩
    z[] <- 0
    flagged[] <- FALSE
  }
  dropped <- if (rule$unit == "item") item %in% item[flagged] else flagged
  abort_if(sum(!dropped) < 3,
           "The screening would keep only ", sum(!dropped), " point(s); at least 3 are needed. Use a larger `cut`.",
           class = "kits_error_fit")
  list(from = from, to = to, item = item, resid = resid, z = z, flagged = flagged, dropped = dropped,
       intercept0 = first$intercept, slope0 = first$slope, cut = cut, rule = rule)
}

#' The linking error by resampling whole items
#'
#' Resamples the items `n_boot` times and redoes the screening and the mean/sigma line on each
#' draw; the SDs of the intercept and slope over the draws are the linking error. It does not
#' track scales or check much; [kits_link_error()] does.
#'
#' @param screen An [algo_screen()] result (its points and screening rule are used again).
#' @param n_boot The number of item resamples.
#' @param seed The random seed (`NULL`: not repeatable); the current random state is left as it is.
#' @return A list: `intercept_error`, `slope_error`, `n_failed` (draws whose screening failed),
#'   `vcov` (2 x 2) and `draws` (one row per kept draw: `intercept`, `slope`).
#' @family algo
#' @examples
#' from <- c(-2, -1.1, 0, 0.9, 2, 3.1, 0.5, 1.4)
#' to <- c(-1, 0, 1.1, 2, 3, 4.2, 1.5, 2.4)
#' err <- algo_screen(from, to, item = rep(1:4, 2)) |>
#'   algo_link_error(n_boot = 100)
#' err$intercept_error
#' @export
algo_link_error <- function(screen, n_boot = 200, seed = 1) {
  # 重抽到兩次的題是同一題的兩份，unit = "item" 時一起丟掉（用 item[i] 分組）；每個點的 SE 跟著點走
  from <- screen$from
  to <- screen$to
  item <- screen$item
  rule <- screen$rule
  items <- unique(item)
  rows <- split(seq_along(item), factor(item, levels = items))
  one_draw <- function() {
    i <- unlist(rows[sample(items, replace = TRUE)], use.names = FALSE)
    link <- tryCatch(algo_screen(from[i], to[i], item[i], screen$cut, rule$spread, rule$unit, rule$se_from[i],
                                 rule$se_to[i]) |> algo_mean_sigma(),
                     error = function(e) list(intercept = NA, slope = NA))
    c(intercept = link$intercept, slope = link$slope)
  }
  draws <- with_seed(seed, t(replicate(n_boot, one_draw())))
  ok <- is.finite(draws[, "intercept"]) & is.finite(draws[, "slope"])
  kept <- draws[ok, , drop = FALSE]
  rownames(kept) <- NULL
  list(intercept_error = stats::sd(kept[, "intercept"]), slope_error = stats::sd(kept[, "slope"]),
       n_failed = sum(!ok), vcov = stats::cov(kept), draws = kept)   # vcov：2 x 2；少於 2 個重抽時 NA
}

#' Apply, invert or compose straight lines
#'
#' A line is a list with `intercept` and `slope` (new value = intercept + slope x old value); a
#' link is one too. A list of lines is composed in the order applied; `invert` gives the inverse;
#' with `x` the line is applied. It does not track scales; [kits_convert()] does.
#'
#' @param line A line, or a list of lines (applied in order).
#' @param x Optional values to convert.
#' @param invert Use the inverse line, x = (y - intercept) / slope.
#' @return Without `x`, the line (`intercept`, `slope`); with `x`, the converted values.
#' @family algo
#' @examples
#' l1 <- list(intercept = 1, slope = 2)
#' l2 <- list(intercept = -1, slope = 0.5)
#' algo_line(l1, x = c(0, 1))
#' algo_line(l1, invert = TRUE)
#' list(l1, l2) |>
#'   algo_line(x = c(0, 1))
#' @export
algo_line <- function(line, x = NULL, invert = FALSE) {
  # 一串線：照套用的順序合成（line_compose()，line.R），從恆等的線 0 + 1 x 開始，所以沒有線時也有
  if (is.null(line$slope)) line <- purrr::reduce(line, line_compose, .init = list(intercept = 0, slope = 1))
  if (invert) line <- list(intercept = -line$intercept / line$slope, slope = 1 / line$slope)   # x = (y - B) / A
  if (is.null(x)) line else line_apply(line, x)
}

#' A shash norm model by maximum likelihood
#'
#' Fits a sinh-arcsinh (shash) distribution of `y` whose four parameters (mu, sigma, epsilon,
#' delta) are polynomials in `x`, by nlminb with the gradient and Hessian from RTMB, from one start of delta.
#' It does not track scales or check its input; [kits_norm_fit()] does (several starts).
#'
#' @param x Transformed age (for example `sqrt(age_m)`); it is standardized inside.
#' @param y The scores.
#' @param degrees The degrees of the four polynomials.
#' @param delta0 The start value of delta.
#' @return A list: the coefficients `mu_est`, `sigma_est`, `epsilon_est`, `delta_est`, `degrees`,
#'   `age_mean`, `age_sd` (of `x`), `value` (-log L), and the convergence record (`convergence` and
#'   `message` of nlminb, `converged`, `grad_max`, `hessian_pd`, `at_bound`).
#' @references
#' Jones, M. C., & Pewsey, A. (2009). Sinh-arcsinh distributions. *Biometrika, 96*(4), 761-780.
#' \doi{10.1093/biomet/asp053}
#' @family algo
#' @examples
#' set.seed(1)
#' age <- runif(600, 3, 36)
#' D <- 20 + 15 * log(age) + rnorm(600, 0, 3)
#' fit <- algo_shash(sqrt(age), D, degrees = c(3, 2, 1, 1), delta0 = 0.5)
#' fit$value
#' @export
algo_shash <- function(x, y, degrees = c(5, 4, 3, 2), delta0 = 1) {
  # 模型、界限和起始值和 cNORM 的 cnorm.shash() 相同；最佳化用 nlminb（RTMB 的梯度和 Hessian），所以和 cNORM 差到約 1e-6。
  # 對數密度（u = (y - mu) / sigma，S = sinh(delta asinh(u) + epsilon)）：
  #   log f = log delta - log sigma - log(2 pi) / 2 - log(1 + u^2) / 2 - S^2 / 2 + log cosh(delta asinh(u) + epsilon)
  # 和 cNORM 一樣，log sigma、epsilon 限在 -10..10，log delta 限在 -2..2（每一筆資料），係數也有同樣的界。
  x_mean <- mean(x)
  x_sd <- stats::sd(x)
  xs <- (x - x_mean) / x_sd
  X_of <- function(degree) outer(xs, 0:degree, "^")
  X <- list(mu = X_of(degrees[1]), sigma = X_of(degrees[2]), epsilon = X_of(degrees[3]),
            delta = X_of(degrees[4]))
  logcosh <- function(a) abs(a) + log1p(exp(-2 * abs(a))) - log(2)
  nll <- function(par) {
    p <- shash_link(X, list(mu = par$c_mu, sigma = par$c_sigma, epsilon = par$c_epsilon,
                            delta = par$c_delta), ad_clamp)
    mu <- p$mu
    sigma <- p$sigma
    epsilon <- p$epsilon
    delta <- p$delta
    u <- (y - mu) / sigma
    arg <- delta * asinh(u) + epsilon
    -sum(log(delta) - log(sigma) - 0.5 * log(2 * pi) - 0.5 * log1p(u^2) - 0.5 * sinh(arg)^2 + logcosh(arg))
  }
  start <- list(c_mu = c(mean(y), rep(0, degrees[1])), c_sigma = c(log(stats::sd(y)), rep(0, degrees[2])),
                c_epsilon = rep(0, degrees[3] + 1), c_delta = c(log(delta0), rep(0, degrees[4])))
  bound <- SHASH_BOUNDS
  lower <- c(rep(-Inf, degrees[1] + 1), rep(-bound[["log_sigma"]], degrees[2] + 1),
             rep(-bound[["epsilon"]], degrees[3] + 1), rep(-bound[["log_delta"]], degrees[4] + 1))
  upper <- -lower
  obj <- RTMB::MakeADFun(nll, start, silent = TRUE)
  fit <- function(p) stats::nlminb(p, obj$fn, obj$gr, obj$he, lower = lower, upper = upper)
  check <- function(r) shash_convergence(r$par, obj$gr, obj$he, lower, upper)
  result <- fit(obj$par)
  conv <- check(result)
  for (restart in 1:3) {   # 停在夾住造成的平坦處（false convergence）時，從那裡再走；nlminb 只會讓 -log L 變小
    if (conv$converged) break
    result <- fit(result$par)
    conv <- check(result)
  }
  if (!conv$converged) {   # 還差一點：和 kits_item() 一樣，在自由的係數上走牛頓步
    result <- shash_newton(result, obj, lower, upper)
    conv <- check(result)
  }
  if (!conv$converged) {
    kits_warn("The norm fit stopped before converging (nlminb code ", result$convergence, "; max |gradient| ",
              signif(conv$grad_max, 3), if (!conv$hessian_pd) "; Hessian not positive definite", "). ",
              lower_degrees_advice(degrees), class = "kits_warning_convergence")
  }
  if (length(conv$at_bound) > 0) {
    kits_warn("Coefficient(s) ", paste(conv$at_bound, collapse = ", "), " of the norm model sit at their limits. ",
              lower_degrees_advice(degrees), class = "kits_warning_norms")
  }
  coef_of <- split(unname(result$par), rep(c("mu", "sigma", "epsilon", "delta"), degrees + 1))
  list(mu_est = coef_of$mu, sigma_est = coef_of$sigma, epsilon_est = coef_of$epsilon, delta_est = coef_of$delta,
       degrees = degrees, age_mean = x_mean, age_sd = x_sd, value = result$objective, convergence = result$convergence,
       message = result$message, converged = conv$converged, grad_max = conv$grad_max, hessian_pd = conv$hessian_pd,
       at_bound = conv$at_bound)
}

#' Same-age z-scores under a shash norm model
#'
#' The parameters of an [algo_shash()] model at `x` (clamped as in the fit), and the z-score
#' sinh(delta asinh((y - mu) / sigma) + epsilon) of each `y`. It does not track scales or check
#' its input; [kits_norm_z()] does (from a norm table).
#'
#' @param fit An [algo_shash()] result.
#' @param x Transformed age, as in the fit.
#' @param y The scores, one per `x`.
#' @return The z-scores (a numeric vector).
#' @family algo
#' @examples
#' set.seed(1)
#' age <- runif(600, 3, 36)
#' D <- 20 + 15 * log(age) + rnorm(600, 0, 3)
#' z <- algo_shash(sqrt(age), D, degrees = c(3, 2, 1, 1), delta0 = 0.5) |>
#'   algo_shash_z(sqrt(age), D)
#' c(mean(z), sd(z))
#' @export
algo_shash_z <- function(fit, x, y) {
  p <- shash_parameters(fit, x)
  shash_z(y, p$mu, p$sigma, p$epsilon, p$delta)
}

#' The raw-score distribution given theta (Lord-Wingersky)
#'
#' P(r | theta) for a set of PCM items by the Lord-Wingersky recursion: adding item j,
#' f_j(x) = sum_k f_(j-1)(x - k) P_jk. It does not track scales or check its input;
#' [kits_precision()] does.
#'
#' @param theta Abilities (a numeric vector).
#' @param betas The step parameters: a list, one numeric vector per item.
#' @return A matrix, one row per theta, one column per raw score 0..R.
#' @references
#' Thissen, D., Pommerich, M., Billeaud, K., & Williams, V. S. L. (1995). Item response theory for
#' scores on tests including polytomous items with ordered responses. *Applied Psychological
#' Measurement, 19*(1), 39-49. \doi{10.1177/014662169501900105}
#' @family algo
#' @examples
#' P <- list(a1 = c(-1, 0), a2 = 0.5) |>
#'   algo_lw(theta = c(-1, 0, 1))
#' P
#' rowSums(P)
#' @export
algo_lw <- function(theta, betas) {
  f <- matrix(1, length(theta), 1)
  for (beta in betas) {
    p <- pcm_probs(theta, beta)   # 這題的類別機率；相乘用 poly_times()（pcm_core.R；和 CML 的 gamma 同一個）
    f <- poly_times(f, map(seq_len(ncol(p)), \(k) p[, k]))
  }
  f
}
