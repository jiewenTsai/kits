# The calibration of kits_item() without kits: the CML of the partial credit model (PCM) in RTMB.
# Parts: "stats" (1, the sufficient statistics), "cml" (2-4), "thr" (5-6, the thresholds); after
# the section "compare", the other cases, each checked with kits: "fixed" (7, fixed item
# parameters) and "fixed_out" (8, fixed items left out of the calibration).
# Only -log L is written; RTMB gives its gradient and Hessian. The names follow kits_algorithm().

library(RTMB)
library(dplyr)
library(tidyr)
library(purrr)

# ==== inputs ====
# Replace with your own data. X: one row per person-wave (row names = visit_id), one column per
# item, scores 0, 1, 2, ...; NA = not given. Here: the language items of the demo data kitm03.
lan <- filter(kits::kitm03$responses, domain == "lan")
wide <- pivot_wider(lan, id_cols = visit_id, names_from = item_id, values_from = item_score)
X <- as.matrix(wide[-1])
rownames(X) <- wide$visit_id

# ==== part: stats ====
# 1 Long data: one row per answer (visit_id, item_id, x); the items in the order of sort()
long <- as_tibble(X, rownames = "visit_id") |>
  pivot_longer(-visit_id, names_to = "item_id", values_to = "x", values_drop_na = TRUE)
items <- sort(unique(long$item_id))
m <- map_dbl(items, \(i) max(long$x[long$item_id == i]))     # m_i: the top score of item i

# n_ig: the persons with a score of at least g on item i, one per step g = 1, ..., m_i
step_item <- rep(items, m)
step_k <- sequence(m)
n_ge <- map2_dbl(step_item, step_k, \(i, g) sum(long$x[long$item_id == i] >= g))

# n_br: the persons with raw score r = 0, 1, ... in each item pattern b (the items answered)
persons <- long |>
  summarise(pattern = paste(sort(match(item_id, items)), collapse = " "), raw = sum(x),
            .by = visit_id)
patterns <- persons |>
  summarise(raws = list(raw), .by = pattern) |>
  mutate(book_items = map(strsplit(pattern, " "), as.integer),
         n_raw = map2(raws, book_items, \(r, its) tabulate(r + 1, nbins = sum(m[its]) + 1)))

# ==== part: cml ====
# 2 gamma_r: the coefficient of z^r in prod_i (1 + exp(-beta_i1) z + exp(-beta_i1 - beta_i2) z^2
#   + ...), a product of polynomials; written with RTMB values, so RTMB can differentiate it
conv <- function(p1, p2) {                                     # coefficients of p1(z) p2(z)
  "[<-" <- ADoverload("[<-")
  out <- numeric(length(p1) + length(p2) - 1)
  for (k in seq_along(p2)) {
    idx <- seq_along(p1) + k - 1
    out[idx] <- out[idx] + p1 * p2[k]
  }
  out
}

# log gamma_r of the items its; gamma_r(beta) = exp(-r c) gamma_r(beta - c): no underflow
steps_of <- split(seq_along(step_item), factor(step_item, levels = items))
log_gamma <- function(beta, its) {
  centre <- mean(beta[unlist(steps_of[its])])
  out <- 1
  for (i in its) out <- conv(out, c(1, exp(-cumsum(beta[steps_of[[i]]] - centre))))
  log(out) - (seq_along(out) - 1) * centre
}

# 3 beta = beta0 + A eta: the free parameters eta, with the sum of all betas = 0 (the origin).
#   With fixed items (kits_item(fixed_params = )), beta0 holds their betas and A has no column
#   for their steps.
beta0 <- numeric(length(step_item))
A <- rbind(diag(length(step_item) - 1), -1)

# 4 Conditional log-likelihood: log L = - sum_ig beta_ig n_ig - sum_b sum_r n_br log gamma_br
nll <- function(par) {
  beta <- beta0 + drop(A %*% par$eta)
  loglik <- -sum(beta * n_ge)
  for (b in seq_len(nrow(patterns))) {
    loglik <- loglik - sum(patterns$n_raw[[b]] * log_gamma(beta, patterns$book_items[[b]]))
  }
  -loglik
}

# The optimizer of kits (fit_ad()): nlminb on the RTMB gradient, then Newton steps with the
#   RTMB Hessian while -log L goes down (at so strict a tolerance nlminb stops a little early)
fit_ad <- function(obj) {
  control <- list(rel.tol = 1e-14, eval.max = 2000, iter.max = 2000)
  par <- nlminb(obj$par, obj$fn, obj$gr, control = control)$par
  value <- obj$fn(par)
  for (step in 1:5) {
    delta <- solve(obj$he(par), drop(obj$gr(par)))
    new_value <- obj$fn(par - delta)
    if (!is.finite(new_value) || new_value > value + 1e-9 + 1e-12 * abs(value)) break
    par <- par - delta
    value <- new_value
    if (max(abs(delta)) < 1e-10) break
  }
  par
}

# The betas, and their covariance V = A H^-1 t(A) (H: the Hessian of -log L at eta)
obj <- MakeADFun(nll, list(eta = numeric(ncol(A))), silent = TRUE)
eta <- fit_ad(obj)
beta <- beta0 + drop(A %*% eta)
V <- A %*% solve(obj$he(eta)) %*% t(A)
cf <- tibble(item_id = step_item, item_score = step_k, beta = beta, beta_se = sqrt(diag(V)))
cf

# ==== part: thr ====
# 5 Thurstonian thresholds: P(X >= k | theta = thr_k) = 0.5, with the PCM
#   P(X = k | theta) = exp(k theta - sum_{g <= k} beta_g) / sum_h exp(h theta - sum_{g <= h} beta_g)
pcm <- function(theta, beta_i) {
  log_num <- cumsum(c(0 * theta, theta - beta_i))             # sum_{g <= k} (theta - beta_g)
  exp(log_num - max(log_num)) / sum(exp(log_num - max(log_num)))
}
p_ge <- function(theta, beta_i, k) sum(pcm(theta, beta_i)[-(1:k)])
thr <- map_dbl(seq_len(nrow(cf)), \(row) {
  beta_i <- cf$beta[cf$item_id == cf$item_id[row]]
  uniroot(\(t) p_ge(t, beta_i, cf$item_score[row]) - 0.5, c(-40, 40), tol = 1e-10)$root
})

# 6 SE of thr by the delta method, sqrt(J V t(J)), J = d thr / d beta from RTMB: the thresholds of
#   an item minimize sum_k (P(X >= k | theta_k) - 0.5)^2; tape$newton() turns this into the map
#   beta -> thr (started at the thresholds above), and its jacobian is J
thr_se <- numeric(nrow(cf))
for (i in items) {
  idx <- which(cf$item_id == i)
  m_i <- length(idx)
  gap <- function(p) {                                         # p = (theta_1..m_i, beta_1..m_i)
    total <- 0
    for (k in seq_len(m_i)) total <- total + (p_ge(p[k], p[m_i + seq_len(m_i)], k) - 0.5)^2
    total
  }
  solve_thr <- MakeTape(gap, c(thr[idx], cf$beta[idx]))$newton(seq_len(m_i), grad_tol = 1e-14)
  J <- matrix(solve_thr$jacobian(cf$beta[idx]), m_i, m_i)
  thr_se[idx] <- sqrt(diag(J %*% V[idx, idx] %*% t(J)))
}
thresholds <- mutate(cf, thr = thr, thr_se = thr_se)
thresholds

# ==== compare ====
# The same numbers as kits_item(), within the tolerance of the optimizer (below 1e-6)
item <- kits::kits_item(X, scale = "lan M03")
want <- coef(item)
diffs <- map_dbl(set_names(c("beta", "beta_se", "thr", "thr_se")),
                 \(v) max(abs(thresholds[[v]] - want[[v]])))
diffs

# ==== part: fixed ====
# 7 Fixed item parameters, kits_item(X, fixed_params = fixed): the fixed steps keep their betas
#   and set the origin (the sum of the betas is no longer 0). beta0 holds them and A has a column
#   for each free step only; nll() reads beta0 and A, so it is taped again. Here: two items fixed
#   at made-up values (in practice, the betas of an earlier calibration)
fixed <- cf |>
  filter(item_id %in% c("lana01", "lanb01")) |>
  transmute(item_id, item_score, beta = round(beta, 1) + 0.5)
key <- paste(step_item, step_k)
fixed_key <- paste(fixed$item_id, fixed$item_score)
beta0 <- if_else(key %in% fixed_key, fixed$beta[match(key, fixed_key)], 0)
A <- diag(length(key))[, !key %in% fixed_key, drop = FALSE]

# The fit of part cml; the fixed steps get SE 0 (their rows of A are 0)
obj_fixed <- MakeADFun(nll, list(eta = numeric(ncol(A))), silent = TRUE)
eta_fixed <- fit_ad(obj_fixed)
V_fixed <- A %*% solve(obj_fixed$he(eta_fixed)) %*% t(A)
cf_fixed <- tibble(item_id = step_item, item_score = step_k, beta = beta0 + drop(A %*% eta_fixed),
                   beta_se = sqrt(diag(V_fixed)))

# The thresholds of part thr on the new betas (their SE: the delta method of part thr, V_fixed)
thr_fixed <- map_dbl(seq_len(nrow(cf_fixed)), \(row) {
  beta_i <- cf_fixed$beta[cf_fixed$item_id == cf_fixed$item_id[row]]
  uniroot(\(t) p_ge(t, beta_i, cf_fixed$item_score[row]) - 0.5, c(-40, 40), tol = 1e-10)$root
})
thresholds_fixed <- mutate(cf_fixed, thr = thr_fixed)
thresholds_fixed

# With kits: the same calibration (below 1e-6)
want_fixed <- coef(kits::kits_item(X, scale = "lan M03", fixed_params = fixed))
diffs_fixed <- map_dbl(set_names(c("beta", "beta_se", "thr")),
                       \(v) max(abs(thresholds_fixed[[v]] - want_fixed[[v]])))
diffs_fixed

# ==== part: fixed_out ====
# 8 CML needs a 0 score and at least two scores on each item (else n_i1 counts all its persons
#   and beta_i1 goes to -Inf). A fixed item without them in these data is left out of the counts:
#   its betas are copied, its thresholds come from them as in part thr, and its SEs are 0. The
#   other items are calibrated as in part fixed. Here: lana00, a fixed item of an earlier wave
#   that every child passed with its top score
easy <- tibble(item_id = "lana00", item_score = 1:2, beta = c(-7, -6))
X_easy <- cbind(lana00 = 2, X)
thr_easy <- map_dbl(easy$item_score, \(k) {
  uniroot(\(t) p_ge(t, easy$beta, k) - 0.5, c(-40, 40), tol = 1e-10)$root
})
thresholds_easy <- bind_rows(mutate(easy, beta_se = 0, thr = thr_easy, thr_se = 0),
                             thresholds_fixed)
head(thresholds_easy)

# With kits: lana00 is left out (with a message) and keeps its fixed betas (below 1e-6)
item_easy <- kits::kits_item(X_easy, scale = "lan M03", fixed_params = bind_rows(easy, fixed))
item_easy$left_out
want_easy <- coef(item_easy)
diffs_easy <- map_dbl(set_names(c("beta", "beta_se", "thr")),
                      \(v) max(abs(thresholds_easy[[v]] - want_easy[[v]])))
diffs_easy
