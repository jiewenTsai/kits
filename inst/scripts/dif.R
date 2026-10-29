# Sex DIF of kits_dif() without kits: the CML of the PCM in each group (RTMB), then the item pairs
# of dexter's DIF() (Bechger & Maris, 2015).
# Parts: "cml" (1, the CML of one group, as cml.R), "dif" (2-3), "overall" (4, the overall test).
# The names follow kits_algorithm().

library(RTMB)
library(dplyr)
library(tidyr)
library(purrr)

# ==== inputs ====
# Replace with your own data. X: one row per person-wave (row names = visit_id), one column per
# item, scores 0, 1, 2, ...; NA = not given. sex: one per row of X; ref: the reference group.
# Here: the language items of kitm03.
lan <- filter(kits::kitm03$responses, domain == "lan")
wide <- pivot_wider(lan, id_cols = visit_id, names_from = item_id, values_from = item_score)
X <- as.matrix(wide[-1])
rownames(X) <- wide$visit_id
sex <- pull(distinct(lan, visit_id, sex), sex, name = visit_id)[rownames(X)]
ref <- "boy"

# ==== part: cml ====
# 1 The CML of the PCM in one group, with the sum of its betas = 0 (the steps of cml.R in one
#   function); it also keeps the observed scores of each item
fit_ad <- function(obj) {                                       # the optimizer of kits (cml.R)
  par <- nlminb(obj$par, obj$fn, obj$gr,
                control = list(rel.tol = 1e-14, eval.max = 2000, iter.max = 2000))$par
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
conv <- function(p1, p2) {                                      # coefficients of p1(z) p2(z)
  "[<-" <- ADoverload("[<-")
  out <- numeric(length(p1) + length(p2) - 1)
  for (k in seq_along(p2)) {
    idx <- seq_along(p1) + k - 1
    out[idx] <- out[idx] + p1 * p2[k]
  }
  out
}

cml <- function(X_g) {
  long <- as_tibble(X_g, rownames = "visit_id") |>
    pivot_longer(-visit_id, names_to = "item_id", values_to = "x", values_drop_na = TRUE)
  items <- sort(unique(long$item_id))
  m <- map_dbl(items, \(i) max(long$x[long$item_id == i]))
  scores <- map_chr(set_names(items), \(i) paste(sort(unique(long$x[long$item_id == i])),
                                                 collapse = ", "))
  step_item <- rep(items, m)
  n_ge <- map2_dbl(step_item, sequence(m), \(i, g) sum(long$x[long$item_id == i] >= g))
  patterns <- long |>
    summarise(pattern = paste(sort(match(item_id, items)), collapse = " "), raw = sum(x),
              .by = visit_id) |>
    summarise(raws = list(raw), .by = pattern) |>
    mutate(its = map(strsplit(pattern, " "), as.integer),
           n_raw = map2(raws, its, \(r, its) tabulate(r + 1, nbins = sum(m[its]) + 1)))

  steps_of <- split(seq_along(step_item), factor(step_item, levels = items))
  log_gamma <- function(beta, its) {
    centre <- mean(beta[unlist(steps_of[its])])
    out <- 1
    for (i in its) out <- conv(out, c(1, exp(-cumsum(beta[steps_of[[i]]] - centre))))
    log(out) - (seq_along(out) - 1) * centre
  }
  A <- rbind(diag(length(step_item) - 1), -1)                   # beta = A eta: sum of betas 0
  nll <- function(par) {
    beta <- drop(A %*% par$eta)
    loglik <- -sum(beta * n_ge)
    for (b in seq_len(nrow(patterns))) {
      loglik <- loglik - sum(patterns$n_raw[[b]] * log_gamma(beta, patterns$its[[b]]))
    }
    -loglik
  }
  obj <- MakeADFun(nll, list(eta = numeric(ncol(A))), silent = TRUE)
  eta <- fit_ad(obj)
  list(item_id = step_item, beta = drop(A %*% eta), V = A %*% solve(obj$he(eta)) %*% t(A),
       scores = scores)
}
fit_ref <- cml(X[which(sex == ref), , drop = FALSE])
head(fit_ref$beta)

# ==== part: dif ====
# 2 Each group on its own; the items compared: answered in both groups, with the same scores
groups <- sort(unique(sex[!is.na(sex)]))
focal <- setdiff(groups, ref)
fits <- map(set_names(groups), \(g) cml(X[which(sex == g), , drop = FALSE]))
common <- intersect(names(fits[[1]]$scores), names(fits[[2]]$scores))
common <- common[fits[[1]]$scores[common] == fits[[2]]$scores[common]]
steps <- map(fits, \(f) {
  keep <- f$item_id %in% common
  list(item_id = f$item_id[keep], beta = f$beta[keep], V = f$V[keep, keep, drop = FALSE])
})

# 3 dif_ig = (beta_ig^(r) - mean beta^(r)) - (beta_ig^(f) - mean beta^(f)), each group centred
#   on its own origin; dif_i = the mean over the steps of item i
rel <- map(steps, \(st) st$beta - mean(st$beta))
dif <- tibble(item_id = steps[[ref]]$item_id, dif = rel[[ref]] - rel[[focal]]) |>
  summarise(dif = mean(dif), .by = item_id)
dif

# ==== part: overall ====
# 4 The overall test on a common origin: centre the compared steps of each group (C = I - 11'/S),
#   chi^2 = d' (C (V1 + V2) C')^-1 d with d = C (beta^(1) - beta^(2)) without its first step
n_s <- length(steps[[1]]$beta)
C <- diag(n_s) - matrix(1 / n_s, n_s, n_s)
d <- drop(C %*% (steps[[1]]$beta - steps[[2]]$beta))
S <- C %*% (steps[[1]]$V + steps[[2]]$V) %*% t(C)
chi2 <- mahalanobis(d[-1], rep(0, n_s - 1), S[-1, -1, drop = FALSE])
overall <- c(stat = chi2, df = n_s - 1, p = pchisq(chi2, n_s - 1, lower.tail = FALSE))
overall

# ==== compare ====
# The same DIF and overall test as kits_dif() (below 1e-6)
want <- kits::kits_dif(X, sex, ref = ref)
diffs <- c(dif = max(abs(dif$dif - want$dif)),
           overall = abs(overall[["stat"]] - attr(want, "overall")$stat))
diffs
