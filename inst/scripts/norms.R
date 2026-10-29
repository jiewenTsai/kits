# The norm of kits_norms() without kits: the shash model by maximum likelihood in RTMB, with the
# polynomial model and bounds of cNORM, minimised by nlminb.
# Parts: "data" (1), "fit" (2-3), "table" (4-5, the parameters and percentiles), "z" (6); after
# the section "compare", the other case, checked with kits: "by" (7, a norm by group).
# Only -log L is written; RTMB gives its gradient and Hessian. The names follow kits_algorithm().

library(RTMB)
library(dplyr)
library(purrr)

# ==== inputs ====
# Replace with your own data. D: the scores, age_m: the age in months of each person-wave, ages:
# the ages of the norm (kits_norm_ages()), and the settings of kits_norms(). Here: the language
# WLE of kitm03 (on the lan M03 scale), in the waves that kits_norm_screen() keeps.
lan <- filter(kits::kitm03$responses, domain == "lan")
wide <- tidyr::pivot_wider(lan, id_cols = visit_id, names_from = item_id, values_from = item_score)
X <- as.matrix(wide[-1])
rownames(X) <- wide$visit_id
persons <- distinct(lan, visit_id, cohort, wave_m, age_m)
scores <- kits::kits_subj(X, kits::kits_item(X, scale = "lan M03"), persons = persons)
use <- filter(scores, kits::kits_norm_screen(scores)$normed)
D <- use$theta
age_m <- use$age_m
ages <- kits::kits_norm_ages(use)

degrees <- c(5, 4, 3, 2)                                        # of mu, sigma, epsilon, delta
pad <- 1
delta_starts <- c(0.5, 1, 1.5, 2)
age_x <- sqrt
pct <- c(3, 5, 10, 15, 25, 50, 75, 85, 90, 95, 97)

# ==== part: data ====
# 1 The data within the ages (plus pad); the table keeps the ages within the data
norm_data <- tibble(D = D, age_m = age_m) |>
  filter(is.finite(D), is.finite(age_m)) |>
  filter(age_m >= min(ages) - pad, age_m <= max(ages) + pad)
ages <- ages[ages >= floor(min(norm_data$age_m)) & ages <= ceiling(max(norm_data$age_m))]
norm_data

# ==== part: fit ====
# 2 Each parameter is a polynomial in x = (age_x(age) - m) / s. The shash log density, with
#   u = (D - mu) / sigma and a = delta asinh(u) + epsilon:
#   log f = log delta - log sigma - log(2 pi) / 2 - log(1 + u^2) / 2 - sinh(a)^2 / 2 + log cosh(a)
m <- mean(age_x(norm_data$age_m))
s <- sd(age_x(norm_data$age_m))
xd <- (age_x(norm_data$age_m) - m) / s
X_of <- function(degree) outer(xd, 0:degree, "^")

# As in cNORM, log sigma and epsilon stay within -10..10 and log delta within -2..2 at every point
clamp <- function(v, bound) (abs(v + bound) - abs(v - bound)) / 2   # pmin(pmax(v, -b), b)
logcosh <- function(a) abs(a) + log1p(exp(-2 * abs(a))) - log(2)
nll <- function(par) {
  mu <- drop(X_of(degrees[1]) %*% par$c_mu)
  sigma <- exp(clamp(drop(X_of(degrees[2]) %*% par$c_sigma), 10))
  epsilon <- clamp(drop(X_of(degrees[3]) %*% par$c_epsilon), 10)
  delta <- exp(clamp(drop(X_of(degrees[4]) %*% par$c_delta), 2))
  u <- (norm_data$D - mu) / sigma
  a <- delta * asinh(u) + epsilon
  -sum(log(delta) - log(sigma) - 0.5 * log(2 * pi) - 0.5 * log1p(u^2) - 0.5 * sinh(a)^2 +
         logcosh(a))
}

# 3 One fit per start of delta, by nlminb with the RTMB gradient and Hessian and the same
#   bounds; keep the largest likelihood
lower <- c(rep(-Inf, degrees[1] + 1), rep(-10, degrees[2] + 1), rep(-10, degrees[3] + 1),
           rep(-2, degrees[4] + 1))
fit_from <- function(d0) {
  start <- list(c_mu = c(mean(norm_data$D), rep(0, degrees[1])),
                c_sigma = c(log(sd(norm_data$D)), rep(0, degrees[2])),
                c_epsilon = rep(0, degrees[3] + 1), c_delta = c(log(d0), rep(0, degrees[4])))
  obj <- MakeADFun(nll, start, silent = TRUE)
  nlminb(obj$par, obj$fn, obj$gr, obj$he, lower = lower, upper = -lower)
}
fits <- map(delta_starts, fit_from)
best <- fits[[which.min(map_dbl(fits, "objective"))]]
coef_of <- split(unname(best$par), rep(c("c_mu", "c_sigma", "c_epsilon", "c_delta"), degrees + 1))
coef_of

# ==== part: table ====
# 4 mu, sigma, epsilon, delta at each age: the polynomials, kept within the limits of the fit
x <- (age_x(ages) - m) / s
poly_x <- function(coef) drop(outer(x, seq_along(coef) - 1, "^") %*% coef)
keep_in <- function(v, bound) pmin(pmax(v, -bound), bound)
norm <- tibble(age_m = ages,
               mu = poly_x(coef_of$c_mu),
               sigma = exp(keep_in(poly_x(coef_of$c_sigma), 10)),
               epsilon = keep_in(poly_x(coef_of$c_epsilon), 10),
               delta = exp(keep_in(poly_x(coef_of$c_delta), 2)))

# 5 Percentiles: P_q = mu + sigma sinh((asinh(qnorm(q / 100)) - epsilon) / delta)
for (q in pct) {
  norm[[paste0("P", q)]] <- with(norm, mu + sigma * sinh((asinh(qnorm(q / 100)) - epsilon) / delta))
}
norm

# ==== part: z ====
# 6 Same-age z = sinh(delta asinh((D - mu) / sigma) + epsilon), the parameters interpolated
#   between the ages of the table; PR = 100 pnorm(z); ages outside the norm get NA
at <- function(v) approx(norm$age_m, v, age_m)$y
z <- sinh(at(norm$delta) * asinh((D - at(norm$mu)) / at(norm$sigma)) + at(norm$epsilon))
summary(z)

# ==== compare ====
# The same norm as kits_norms() (below 1e-4, the convergence tolerance of kits)
want <- kits::kits_norms(D, age_m, ages = kits::kits_norm_ages(use), degrees = degrees, pad = pad,
                         delta_starts = delta_starts, pct = pct)
diffs <- c(map_dbl(set_names(c("mu", "P3", "P50", "P97")), \(v) max(abs(norm[[v]] - want[[v]]))),
           z = max(abs(z - kits::kits_norm_z(D, age_m, want)), na.rm = TRUE))
diffs

# ==== part: by ====
# 7 A norm by sex, kits_norms(D, age_m, ages, by = data.frame(sex = sex)): steps 1-4 (and 5) on
#   the person-waves of each group alone, with the same ages and settings; each person-wave gets
#   its z on the curves of its group. fit_from() and nll() read norm_data and xd, so the loop
#   sets them for each group (afterwards they hold the last group)
sex <- lan$sex[match(use$visit_id, lan$visit_id)]
norm_by <- list()
z_by <- rep(NA_real_, length(D))
for (g in sort(unique(sex))) {
  rows <- sex == g
  norm_data <- tibble(D = D[rows], age_m = age_m[rows]) |>                      # 1
    filter(is.finite(D), is.finite(age_m), age_m >= min(ages) - pad, age_m <= max(ages) + pad)
  ages_g <- ages[ages >= floor(min(norm_data$age_m)) & ages <= ceiling(max(norm_data$age_m))]
  m <- mean(age_x(norm_data$age_m))                                             # 2-3
  s <- sd(age_x(norm_data$age_m))
  xd <- (age_x(norm_data$age_m) - m) / s
  fits_g <- map(delta_starts, fit_from)
  best_g <- fits_g[[which.min(map_dbl(fits_g, "objective"))]]
  c_g <- split(unname(best_g$par), rep(c("mu", "sigma", "epsilon", "delta"), degrees + 1))
  x_g <- (age_x(ages_g) - m) / s                                                # 4
  poly_g <- function(coef) drop(outer(x_g, seq_along(coef) - 1, "^") %*% coef)
  norm_g <- tibble(sex = g, age_m = ages_g, mu = poly_g(c_g$mu),
                   sigma = exp(keep_in(poly_g(c_g$sigma), 10)),
                   epsilon = keep_in(poly_g(c_g$epsilon), 10),
                   delta = exp(keep_in(poly_g(c_g$delta), 2)))
  at <- function(v) approx(norm_g$age_m, v, age_m[rows])$y                     # 6
  z_by[rows] <- sinh(at(norm_g$delta) * asinh((D[rows] - at(norm_g$mu)) / at(norm_g$sigma)) +
                       at(norm_g$epsilon))
  norm_by[[g]] <- norm_g
}
norm_by <- bind_rows(norm_by)
norm_by

# With kits: the same norm by sex and the same z (below 1e-4)
want_by <- kits::kits_norms(D, age_m, ages = kits::kits_norm_ages(use), degrees = degrees,
                            pad = pad, delta_starts = delta_starts, pct = pct,
                            by = data.frame(sex = sex))
z_kits <- kits::kits_norm_z(D, age_m, want_by, by = data.frame(sex = sex))
diffs_by <- c(map_dbl(set_names(c("mu", "sigma", "epsilon", "delta")),
                      \(v) max(abs(norm_by[[v]] - want_by[[v]]))),
              z = max(abs(z_by - z_kits), na.rm = TRUE))
diffs_by
