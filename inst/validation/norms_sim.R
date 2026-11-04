# Simulated norm data and the cNORM reference, shared by engine.R and reference/make_norms_cnorm.R
# (no real data).
norm_sim <- function() {
  set.seed(5)
  age <- stats::runif(2000, 3, 36)
  D <- 20 + 15 * log(age) + stats::rnorm(2000, 0, 3)
  ages <- seq(4, 35, by = 0.5)
  keep <- age >= min(ages) & age <= max(ages)   # pad = 0: the fit uses the children within the ages
  list(age = age[keep], D = D[keep], ages = ages,
       age_row = round(age[keep] * 2) / 2,        # each child at the nearest row, so no interpolation
       degrees = c(5, 4, 3, 2), delta_starts = c(0.5, 1, 1.5, 2), pct = c(3, 10, 50, 90, 97))
}

# cNORM's shash fit of the same data: the best of the delta starts, as kits_norms() keeps
cnorm_ref <- function(sim) {
  fits <- lapply(sim$delta_starts, function(d0) suppressWarnings(
    cNORM::cnorm.shash(age = sqrt(sim$age), score = sim$D, mu_degree = sim$degrees[1],
                       sigma_degree = sim$degrees[2], epsilon_degree = sim$degrees[3],
                       delta_degree = sim$degrees[4], delta = d0, plot = FALSE)))
  best <- fits[[which.min(vapply(fits, function(f) f$result$value, numeric(1)))]]
  # a percentile of cNORM is the score whose norm score (T scale, 50 + 10 z) is 50 + 10 qnorm(p)
  pct_at <- function(p) vapply(sim$ages, function(a) stats::uniroot(
    function(s) stats::predict(best, sqrt(a), s) - (50 + 10 * stats::qnorm(p)), c(0, 200), tol = 1e-10)$root,
    numeric(1))
  pct <- sapply(sim$pct, function(p) pct_at(p / 100))
  colnames(pct) <- paste0("P", sim$pct)
  list(engine = paste("cNORM", utils::packageVersion("cNORM")), neg_loglik = best$result$value, pct = pct,
       z = (stats::predict(best, sqrt(sim$age_row), sim$D) - 50) / 10)
}
