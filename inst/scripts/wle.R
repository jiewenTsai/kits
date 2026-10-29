# The person scores of kits_subj() without kits: the WLE (Warm, 1989) in RTMB, its SE, person fit.
# Parts: "wle" (1-4, theta and its SE), "person_fit" (5-7, infit, outfit and flag); after the
# section "compare", the other case, checked with kits: "convert" (8, converted scores).
# Only the weighted likelihood is written; RTMB gives every derivative and the Newton method.
# The names follow kits_algorithm().

library(RTMB)
library(dplyr)
library(tidyr)
library(purrr)

# ==== inputs ====
# Replace with your own data. X: one row per person-wave (row names = visit_id), one column per
# item, scores 0, 1, 2, ...; NA = not given. cf: the item parameters (item_id, item_score, beta),
# from the calibration (cml.R, or coef(kits_item(X))). Here: the language items of kitm03.
lan <- filter(kits::kitm03$responses, domain == "lan")
wide <- pivot_wider(lan, id_cols = visit_id, names_from = item_id, values_from = item_score)
X <- as.matrix(wide[-1])
rownames(X) <- wide$visit_id
item <- kits::kits_item(X, scale = "lan M03")
cf <- select(coef(item), item_id, item_score, beta)

# ==== part: wle ====
# 1 The raw score is sufficient: log L(theta) = r theta - sum_i log S_i(theta), with
#   S_i = sum_k exp(k theta - sum_{g <= k} beta_ig). So there is one theta per item pattern x raw
#   score (a combination), each from its first person
beta_of <- cf |>
  arrange(item_score) |>
  summarise(beta = list(beta), .by = item_id) |>
  pull(beta, name = item_id)
long <- as_tibble(X, rownames = "visit_id") |>
  pivot_longer(-visit_id, names_to = "item_id", values_to = "x", values_drop_na = TRUE)
persons <- long |>
  summarise(set = paste(sort(item_id), collapse = " "), r = sum(x), .by = visit_id)
combos <- distinct(persons, set, r)
sets <- strsplit(combos$set, " ", fixed = TRUE)
r <- combos$r

# 2 f(theta) = log L(theta) + 1/2 log I(theta), the PCM information I = sum_i Var(X_i | theta);
#   every combination in one function of the vector theta, written with RTMB values
rows_of <- map(set_names(names(beta_of)), \(i) which(map_lgl(sets, \(set) i %in% set)))
item_terms <- function(theta, beta_i) {                         # log S_i(theta), Var(X_i | theta)
  k <- seq(0, length(beta_i))
  eta <- outer(theta, k) - rep(c(0, cumsum(beta_i)), each = length(theta))
  S <- drop(exp(eta) %*% rep(1, length(k)))
  P <- exp(eta) / S                                             # P(X_i = k | theta)
  E <- drop(P %*% k)
  list(log_S = log(S), var = drop((P * outer(-E, k, "+")^2) %*% rep(1, length(k))))
}
f_parts <- function(theta) {                                    # sum_i log S_i and I
  "[<-" <- ADoverload("[<-")
  log_S <- numeric(length(theta))
  info <- numeric(length(theta))
  for (i in names(rows_of)) {
    rows <- rows_of[[i]]
    terms <- item_terms(theta[rows], beta_of[[i]])
    log_S[rows] <- log_S[rows] + terms$log_S
    info[rows] <- info[rows] + terms$var
  }
  list(log_S = log_S, info = info)
}

# 3 Maximize f with RTMB's Newton method (tape$newton(); the Hessian is diagonal), from the start
#   value of TAM; with two maxima (few items far apart) it keeps the one it reaches (kits: the
#   higher one)
mean_beta <- map_dbl(sets, \(set) mean(unlist(beta_of[set])))
max_raw <- map_dbl(sets, \(set) sum(lengths(beta_of[set])))
theta0 <- mean_beta + log((r + 0.5) / (max_raw - r + 1))
neg_f <- function(d) {
  theta <- theta0 + d
  parts <- f_parts(theta)
  -sum(r * theta - parts$log_S + 0.5 * log(parts$info))
}
solve_f <- MakeTape(neg_f, numeric(length(r)))$newton(seq_along(r), grad_tol = 1e-12)
theta <- theta0 + solve_f(numeric(0))

# 4 SE = 1 / sqrt(I(theta)), the asymptotic SE of the WLE (Warm, 1989); then each person
combos <- mutate(combos, theta = theta, se = 1 / sqrt(f_parts(theta)$info))
subj <- persons |>
  left_join(combos, by = c("set", "r")) |>
  transmute(visit_id, raw = r, theta, se)
head(subj)

# ==== part: person_fit ====
# 5 raw_max: the top raw score of the items each person answered
raw_max <- long |>
  mutate(top = lengths(beta_of[item_id])) |>
  summarise(raw_max = sum(top), .by = visit_id)

# 6 Person fit at the WLE, one row per answer: E = sum_k k P_k, W = sum_k k^2 P_k - E^2;
#   infit = sum (x - E)^2 / sum W, outfit = mean of (x - E)^2 / W
moments <- function(theta, beta_i) {
  k <- seq(0, length(beta_i))
  eta <- outer(theta, k) - rep(c(0, cumsum(beta_i)), each = length(theta))
  P <- exp(eta - apply(eta, 1, max))
  P <- P / rowSums(P)
  tibble(E = drop(P %*% k), W = drop(P %*% k^2) - drop(P %*% k)^2)
}
person_fit <- long |>
  left_join(select(subj, visit_id, theta), by = "visit_id") |>
  mutate(moments(theta, beta_of[[item_id[1]]]), .by = item_id) |>
  summarise(infit = sum((x - E)^2) / sum(W), outfit = mean((x - E)^2 / W), .by = visit_id)

# 7 flag: zero, perfect, misfit (infit or outfit > 2), ok; no fit for a zero or a perfect score
subj <- subj |>
  left_join(raw_max, by = "visit_id") |>
  left_join(person_fit, by = "visit_id") |>
  mutate(inner = raw > 0 & raw < raw_max,
         infit = if_else(inner, infit, NA),
         outfit = if_else(inner, outfit, NA),
         flag = case_when(raw == 0 ~ "zero",
                          raw == raw_max ~ "perfect",
                          .default = if_else(pmax(infit, outfit) > 2, "misfit", "ok"))) |>
  select(-inner)
head(subj)

# ==== compare ====
# The same numbers as kits_subj() (below 1e-6), person by person
want <- kits::kits_subj(X, item)
want <- want[match(subj$visit_id, want$visit_id), ]
diffs <- map_dbl(set_names(c("theta", "se", "infit", "outfit")),
                 \(v) max(abs(subj[[v]] - want[[v]]), na.rm = TRUE))
diffs
all(subj$flag == want$flag)

# ==== part: convert ====
# 8 Converted scores, kits_convert(scores, link): the line of a link (link.R) on each score,
#   new theta = intercept + slope theta, new SE = slope SE. With the bootstrap draws of the link
#   (n_boot > 0), se_link is the SD of intercept + slope theta over the draws (at the old theta),
#   and se_total = sqrt(se^2 + se_link^2). Person fit and flag do not change.
convert <- function(subj, intercept, slope, draws) {
  V <- cov(draws)                                         # of the intercepts and slopes drawn
  subj |>
    mutate(se_link = sqrt(V[1, 1] + 2 * V[1, 2] * theta + V[2, 2] * theta^2),
           theta = intercept + slope * theta,
           se = slope * se,
           se_total = sqrt(se^2 + se_link^2))
}

# With kits: a link of these scores to the language scale of kitm36 (kits_linking(item, item36),
#   link.R), and the converted scores of kits_subj() (below 1e-6)
lan36 <- filter(kits::kitm36$responses, domain == "lan") |>
  pivot_wider(id_cols = visit_id, names_from = item_id, values_from = item_score)
X36 <- as.matrix(lan36[-1])
rownames(X36) <- lan36$visit_id
link <- kits::kits_linking(item, kits::kits_item(X36, scale = "lan M36"), n_boot = 200)
subj_36 <- convert(subj, link$intercept, link$slope, link$draws)
want_36 <- kits::kits_convert(kits::kits_subj(X, item), link)            # on lan M36
want_36 <- want_36[match(subj_36$visit_id, want_36$visit_id), ]
diffs_36 <- map_dbl(set_names(c("theta", "se", "se_link", "se_total")),
                   \(v) max(abs(subj_36[[v]] - want_36[[v]])))
diffs_36

# ==== part: warm ====
# 9 Warm's (1989) own route, for comparison only: the root of the score equation
#   d log L / d theta + J / (2 I) = 0, I = sum_ik P'_ik^2 / P_ik, J = sum_ik P'_ik P''_ik / P_ik,
#   with uniroot(). P', P'' and d log S_i / d theta come from RTMB tapes (jacfun()); then
#   d log L / d theta = r - sum_i d log S_i / d theta. For the PCM J = dI / dtheta, so the root
#   is the maximum of f in step 3 (the Jeffreys modal estimate): the same theta
item_tapes <- map(beta_of, \(beta_i) {
  k <- seq(0, length(beta_i))
  eta <- function(t) t * k - c(0, cumsum(beta_i))
  P <- MakeTape(\(t) exp(eta(t)) / sum(exp(eta(t))), 0)            # P(X_i = k | theta)
  P1 <- P$jacfun()
  log_S <- MakeTape(\(t) log(sum(exp(eta(t)))), 0)
  list(P = P, P1 = P1, P2 = P1$jacfun(), dlog_S = log_S$jacfun())
})
warm_score <- function(theta, set, r) {
  score <- r
  I <- 0
  J <- 0
  for (tp in item_tapes[set]) {
    p <- drop(tp$P(theta))
    p1 <- drop(tp$P1(theta))
    score <- score - drop(tp$dlog_S(theta))
    I <- I + sum(p1^2 / p)
    J <- J + sum(p1 * drop(tp$P2(theta)) / p)
  }
  score + J / (2 * I)
}
warm <- map2_dbl(sets, r, \(set, r_j) {
  uniroot(warm_score, c(-15, 15), set = set, r = r_j, tol = 1e-12)$root
})
max(abs(warm - combos$theta))                                     # the theta of step 3

# With kits: the same root as kits::algo_jme(), the WLE of kits_subj() (about 1e-8)
jme <- kits::algo_jme(beta_of, raw = r, items = sets)
max(abs(warm - jme$theta))
