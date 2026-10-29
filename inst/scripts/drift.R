# The drift check of kits_drift() without kits: which items of a new wave behave differently from
# the earlier calibration ref (a mean/sigma line with one round of screening; base R and dplyr).
# kits_drift() is a free kits_item() of the new data, then kits_linking() of the judged steps.
# Parts: "points" (1-2, the common steps that can be judged), "line" (3), "items" (4-5).
# The names follow kits_algorithm().

library(dplyr)
library(tidyr)
library(purrr)

# ==== inputs ====
# Replace with your own data. X_new: the new responses (one row per person-wave, one column per
# item); ref: the thresholds of the earlier calibration (coef(kits_item())); free: the thresholds
# of the new data calibrated freely (kits_item(); only items with a 0 score, two or more scores and
# no gap). Here: the language items of kitm36 against the calibration of kitm03.
wide_of <- function(responses) {
  w <- pivot_wider(filter(responses, domain == "lan"), id_cols = visit_id,
                   names_from = item_id, values_from = item_score)
  X <- as.matrix(w[-1])
  rownames(X) <- w$visit_id
  X
}
item_ref <- kits::kits_item(wide_of(kits::kitm03$responses), scale = "lan M03")
ref <- as.data.frame(coef(item_ref))
X_new <- wide_of(kits::kitm36$responses)
free <- as.data.frame(coef(kits::kits_item(X_new, scale = "lan M36")))

# ==== part: points ====
# 1 Why an item in both cannot be judged ("" when it can): a higher top score than in ref, no 0
#   score or one score only, or a gap in its scores
long <- as_tibble(X_new, rownames = "visit_id") |>
  pivot_longer(-visit_id, names_to = "item_id", values_to = "x", values_drop_na = TRUE)
items <- long |>
  filter(item_id %in% ref$item_id) |>
  summarise(max_new = max(x), can_fit = any(x == 0) && n_distinct(x) > 1,
            gap = !all(seq(min(x), max(x)) %in% x), .by = item_id) |>
  left_join(summarise(ref, max_ref = max(item_score), .by = item_id), by = "item_id") |>
  mutate(reason = case_when(max_new > max_ref ~ "score range differs",
                            !can_fit ~ "no 0 score or one score only",
                            gap ~ "score gap",
                            .default = "")) |>
  arrange(item_id)

# 2 The common steps: the free threshold in the new data (thr_new) and in ref (thr_ref), and the
#   persons below (score < k) and above (score >= k) the step; judged with min_n of each
min_n <- 30
points <- free |>
  select(item_id, item_score, thr_new = thr) |>
  inner_join(select(ref, item_id, item_score, thr_ref = thr), by = c("item_id", "item_score")) |>
  filter(item_id %in% items$item_id[items$reason == ""]) |>
  arrange(item_id, item_score) |>
  mutate(n_below = map2_dbl(item_id, item_score, \(i, k) sum(long$x[long$item_id == i] < k)),
         n_above = map2_dbl(item_id, item_score, \(i, k) sum(long$x[long$item_id == i] >= k)),
         judged = pmin(n_below, n_above) >= min_n)
items <- mutate(items, reason = if_else(reason == "" & !item_id %in% points$item_id[points$judged],
                                        "too little information", reason))

# ==== part: line ====
# 3 mean/sigma through the judged steps, new data -> ref (thr_ref = intercept + slope thr_new),
#   with one round of screening: drop the steps with |e| > cut SD(e) (as kits_linking())
cut <- 2
mean_sigma <- function(from, to) {
  slope <- sd(to) / sd(from)
  c(intercept = mean(to) - slope * mean(from), slope = slope)
}
use <- filter(points, judged)
line0 <- mean_sigma(use$thr_new, use$thr_ref)
e <- use$thr_ref - (line0[["intercept"]] + line0[["slope"]] * use$thr_new)
dropped <- if (sd(e) < 1e-8) rep(FALSE, length(e)) else abs(e) > cut * sd(e)
line <- mean_sigma(use$thr_new[!dropped], use$thr_ref[!dropped])
line

# ==== part: items ====
# 4 The residuals from the final line; a judged step that the screening dropped has drifted
points <- mutate(points,
                 resid = thr_ref - (line[["intercept"]] + line[["slope"]] * thr_new),
                 dropped = replace(rep(NA, n()), judged, dropped))

# 5 Each item: drift (one of its judged steps dropped); suggested for fixing when judged and not
#   drifted, or not judged for too little information, for no 0 score or one score only, or for a
#   higher top score than in ref (kits_item(fixed_params = ) keeps its old steps, estimates the new)
by_item <- points |>
  summarise(drift = any(dropped, na.rm = TRUE), max_abs_resid = max(abs(resid)), .by = item_id)
items <- items |>
  left_join(by_item, by = "item_id") |>
  mutate(judged = reason == "",
         drift = if_else(judged, drift, NA),
         suggest_fix = (judged & !drift) |
           reason %in% c("too little information", "no 0 score or one score only",
                         "score range differs")) |>
  select(item_id, max_abs_resid, drift, suggest_fix, reason)
items

# ==== compare ====
# The same line and the same items as kits_drift() (below 1e-6)
want <- kits::kits_drift(X_new, item_ref)
diffs <- c(line = max(abs(line - c(want$intercept, want$slope))),
           resid = max(abs(items$max_abs_resid - want$items$max_abs_resid), na.rm = TRUE))
diffs
identical(items$drift, want$items$drift) && identical(items$suggest_fix, want$items$suggest_fix)
