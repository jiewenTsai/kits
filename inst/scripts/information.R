# The two functions of theta packed in a calibration, without kits: kits_information() (the test
# information of the PCM) and kits_expected_score() (the expected raw score, the TCC).
# Parts: "information" (1), "expected_score" (2-3; RTMB checks that its slope is the information).
# The names follow kits_algorithm().

library(RTMB)
library(dplyr)
library(purrr)

# ==== inputs ====
# Replace with your own item parameters: cf, the betas (item_id, item_score, beta) of a
# calibration (cml.R, or coef(kits_item(X))). Here: the language items of kitm03.
lan <- filter(kits::kitm03$responses, domain == "lan")
wide <- tidyr::pivot_wider(lan, id_cols = visit_id, names_from = item_id, values_from = item_score)
X <- as.matrix(wide[-1])
rownames(X) <- wide$visit_id
item <- kits::kits_item(X, scale = "lan M03")
cf <- select(coef(item), item_id, item_score, beta)
theta <- seq(-6, 6, by = 2)

# ==== part: information ====
# 1 I(theta) = sum_i Var(X_i | theta), Var(X_i | theta) = sum_k k^2 P_ik - (sum_k k P_ik)^2, with
#   P_ik = P(X_i = k | theta) of the PCM
beta_of <- cf |>
  arrange(item_score) |>
  summarise(beta = list(beta), .by = item_id) |>
  pull(beta, name = item_id)
pcm <- function(theta, beta_i) {
  log_num <- cumsum(c(0 * theta, theta - beta_i))             # sum_{g <= k} (theta - beta_g)
  exp(log_num - max(log_num)) / sum(exp(log_num - max(log_num)))
}
item_variance <- function(t, beta_i) {
  P <- pcm(t, beta_i)
  k <- seq_along(P) - 1
  sum(k^2 * P) - sum(k * P)^2
}
information <- function(theta) {
  map_dbl(theta, \(t) sum(map_dbl(beta_of, \(beta_i) item_variance(t, beta_i))))
}
information(theta)

# ==== part: expected_score ====
# 2 The expected raw score: sum_i E(X_i | theta), E(X_i | theta) = sum_k k P_ik
expected_one <- function(t) {
  total <- 0
  for (beta_i in beta_of) total <- total + sum(seq_along(beta_i) * pcm(t, beta_i)[-1])
  total
}
expected_score <- function(theta) map_dbl(theta, expected_one)
expected_score(theta)

# 3 Its slope is the information: d E(X | theta) / d theta = I(theta), the slope from RTMB
slope <- map_dbl(theta, \(t) MakeTape(expected_one, t)$jacobian(t))
max(abs(slope - information(theta)))

# ==== compare ====
# The same functions as kits_information() and kits_expected_score() (below 1e-10)
diffs <- c(information = max(abs(information(theta) - kits::kits_information(item)(theta))),
           expected = max(abs(expected_score(theta) - kits::kits_expected_score(item)(theta))))
diffs
