# A link of kits_linking() without kits: mean/sigma with one round of screening (base R).
# Parts: "points" (1), "line" (2-3, mean/sigma and the screening), "error" (4-5, the linking
# error), "domains" (6, the domain check). The demo links two domains to D through anchor items.
# After the section "compare", the other cases, each checked with kits: "slope" (7, a link
# between two calibrations and its slope against 1), "rules" (8, the other screening rules).
# The names follow kits_algorithm().

library(dplyr)
library(purrr)

# ==== inputs ====
# Replace with your own data. thresholds: the thresholds (item_id, item_score, thr) of the
# calibrations, with a domain column; anchors: the anchor items (item_id, tau, the difficulty on
# D). Here: the language and motor items of kitm03, calibrated with kits (cml.R).
wide_of <- function(d) {
  r <- filter(kits::kitm03$responses, domain == d)
  w <- tidyr::pivot_wider(r, id_cols = visit_id, names_from = item_id, values_from = item_score)
  X <- as.matrix(w[-1])
  rownames(X) <- w$visit_id
  X
}
item_lan <- kits::kits_item(wide_of("lan"), scale = "lan M03")
item_motor <- kits::kits_item(wide_of("motor"), scale = "motor M03")
thresholds <- bind_rows(lan = as.data.frame(coef(item_lan)),
                        motor = as.data.frame(coef(item_motor)), .id = "domain")
anchors <- kits::kitm03$anchors

# ==== part: points ====
# 1 One point per anchor item: the threshold of its pass step (score >= 1) on the KIT scale
#   (from) and its tau on D (to); the items in the order of sort()
points <- thresholds |>
  filter(item_score == 1) |>
  inner_join(select(anchors, item_id, tau), by = "item_id") |>
  arrange(item_id) |>
  transmute(id = item_id, item = item_id, from = thr, to = tau, domain)
points
# Between two calibrations (kits_linking(item36, item03)) each common step is a point instead:
#   inner_join() the two threshold tables by item_id and item_score; from = thr of M36, to = M03

# ==== part: line ====
# 2 mean/sigma: slope = SD(to) / SD(from), intercept = mean(to) - slope mean(from)
mean_sigma <- function(from, to) {
  slope <- sd(to) / sd(from)
  c(intercept = mean(to) - slope * mean(from), slope = slope)
}

# 3 Screening once: the line through all points, e = to - (intercept0 + slope0 from); drop the
#   points with |e| > cut SD(e), then the line through the kept points
#   (kits_link_screen(spread = "mad") uses mad(e); unit = "item" drops every step of an item)
cut <- 2
screen <- function(from, to) {
  line0 <- mean_sigma(from, to)
  e <- to - (line0[["intercept"]] + line0[["slope"]] * from)
  if (sd(e) < 1e-8) return(rep(FALSE, length(e)))              # all on a line: nothing dropped
  abs(e) > cut * sd(e)
}
dropped <- screen(points$from, points$to)
line <- mean_sigma(points$from[!dropped], points$to[!dropped])
line
points$id[dropped]

# Convert (kits_convert()): new = intercept + slope x old, new SE = slope x old SE
line[["intercept"]] + line[["slope"]] * c(-2, 0, 2)

# ==== part: error ====
# 4 Linking error: resample whole items n_boot times and redo the screening and the line; a draw
#   that cannot be fitted (fewer than 3 points, or one value only) is left out and counted
n_boot <- 200
seed <- 1
link_once <- function(from, to) {
  if (length(from) < 3 || sd(from) == 0 || sd(to) == 0) return(c(intercept = NA, slope = NA))
  drop <- screen(from, to)
  if (sum(!drop) < 3) return(c(intercept = NA, slope = NA))
  mean_sigma(from[!drop], to[!drop])
}

items <- unique(points$item)
rows <- split(seq_len(nrow(points)), factor(points$item, levels = items))
set.seed(seed)
draws <- t(replicate(n_boot, {
  i <- unlist(rows[sample(items, replace = TRUE)], use.names = FALSE)
  link_once(points$from[i], points$to[i])
}))
ok <- is.finite(draws[, "intercept"]) & is.finite(draws[, "slope"])
c(intercept_error = sd(draws[ok, "intercept"]), slope_error = sd(draws[ok, "slope"]),
  n_failed = sum(!ok))
# Between two calibrations kits also compares the slope with 1 (a PCM link needs only a shift;
#   Kolen & Brennan, 2014, p. 229): a note when |(slope - 1) / slope_error| > 2

# 5 The linking error at theta: the SD of intercept + slope theta over the draws (Kolen &
#   Brennan, 2014, pp. 249-252)
V_draws <- cov(draws[ok, ])
se_link <- function(theta) {
  sqrt(V_draws[1, 1] + theta^2 * V_draws[2, 2] + 2 * theta * V_draws[1, 2])
}
se_link(c(-2, 0, 2))

# ==== part: domains ====
# 6 Domain check: on the kept points, one slope slope_s shared by the domains with an intercept
#   intercept_d per domain (it is not the slope of the line above)
kept <- points |>
  filter(!dropped) |>
  mutate(dx = from - mean(from), dy = to - mean(to), .by = domain)
slope_s <- sqrt(sum(kept$dy^2) / sum(kept$dx^2))
domains <- kept |>
  summarise(n = n(), intercept = mean(to) - slope_s * mean(from), .by = domain) |>
  arrange(domain)
domains

# ==== compare ====
# The same link as kits_linking() (below 1e-6)
link <- kits::kits_linking(list(lan = item_lan, motor = item_motor), anchors, n_boot = n_boot)
diffs <- c(line = max(abs(line - c(link$intercept, link$slope))),
           error = abs(sd(draws[ok, "intercept"]) - link$intercept_error),
           domains = max(abs(domains$intercept - link$domains$intercept)),
           slope_s = abs(slope_s - link$domains$slope[1]))
diffs
identical(points$id[dropped], link$dropped)

# ==== part: slope ====
# 7 Between two calibrations, kits_linking(item36, item03): one point per common step, its
#   threshold on M36 (from) and on M03 (to), with the SE of both. Steps 2-4 are the same; the
#   slope is also compared with 1, z = (slope - 1) / slope_error, since a PCM link needs only a
#   shift (Kolen & Brennan, 2014, p. 229); kits notes |z| > 2. Here: the language items of kitm36

# With kits: the calibration of the language items of kitm36 (cml.R)
lan36 <- filter(kits::kitm36$responses, domain == "lan") |>
  tidyr::pivot_wider(id_cols = visit_id, names_from = item_id, values_from = item_score)
X36 <- as.matrix(lan36[-1])
rownames(X36) <- lan36$visit_id
item36 <- kits::kits_item(X36, scale = "lan M36")

points36 <- as.data.frame(coef(item36)) |>
  inner_join(filter(thresholds, domain == "lan"), by = c("item_id", "item_score"),
             suffix = c("_from", "_to")) |>
  arrange(item_id, item_score) |>
  transmute(id = paste0(item_id, ":", item_score), item = item_id, from = thr_from, to = thr_to,
            se_from = thr_se_from, se_to = thr_se_to)
dropped36 <- screen(points36$from, points36$to)
line36 <- mean_sigma(points36$from[!dropped36], points36$to[!dropped36])

# Resample whole items as in part error, then z
rows36 <- split(seq_len(nrow(points36)), factor(points36$item, levels = unique(points36$item)))
set.seed(seed)
draws36 <- t(replicate(n_boot, {
  i <- unlist(rows36[sample(names(rows36), replace = TRUE)], use.names = FALSE)
  link_once(points36$from[i], points36$to[i])
}))
slope_error36 <- sd(draws36[, "slope"], na.rm = TRUE)
z_slope <- (line36[["slope"]] - 1) / slope_error36
c(line36, slope_error = slope_error36, z = z_slope)

# With kits: the same link and its row "slope vs 1" (below 1e-6)
link36 <- kits::kits_linking(item36, item_lan, n_boot = n_boot)
z_kits <- filter(as.data.frame(link36$diagnostics), check == "slope vs 1")$value
diffs_slope <- c(line = max(abs(line36 - c(link36$intercept, link36$slope))),
                 z = abs(z_slope - z_kits))
diffs_slope

# ==== part: rules ====
# 8 The other rules of the screening (kits_link_screen(), kits_linking(spread = , unit = )):
#   spread = "mad": |e| > cut x 1.4826 MAD(e), a spread that the outliers do not inflate;
#   spread = "se": |e| / sqrt(se_to^2 + slope0^2 se_from^2) > cut, each point by its own SE
#   (an anchor's tau has no SE: se_to = 0); unit = "item": drop every step of an item with a
#   flagged step. On the points of part slope (several steps per item, SEs on both sides)
screen_rule <- function(pts, spread, unit) {
  line0 <- mean_sigma(pts$from, pts$to)
  e <- pts$to - (line0[["intercept"]] + line0[["slope"]] * pts$from)
  s <- switch(spread, sd = sd(e), mad = mad(e),                 # mad() has the 1.4826
              se = sqrt(pts$se_to^2 + line0[["slope"]]^2 * pts$se_from^2))
  flagged <- abs(e) > cut * s
  if (unit == "item") pts$item %in% pts$item[flagged] else flagged
}
rules <- tibble(spread = c("mad", "se", "sd", "mad"), unit = c("step", "step", "item", "item"))
rules <- rules |>
  mutate(dropped = map2(spread, unit, \(s, u) screen_rule(points36, s, u)),
         n_dropped = map_int(dropped, sum),
         line = map(dropped, \(d) mean_sigma(points36$from[!d], points36$to[!d])),
         intercept = map_dbl(line, "intercept"),
         slope = map_dbl(line, "slope"))
select(rules, spread, unit, n_dropped, intercept, slope)

# With kits: the same points dropped and the same lines (below 1e-6)
links <- map2(rules$spread, rules$unit, \(s, u) {
  kits::kits_linking(item36, item_lan, spread = s, unit = u)
})
diffs_rules <- abs(c(rules$intercept, rules$slope) -
                     c(map_dbl(links, "intercept"), map_dbl(links, "slope")))
max(diffs_rules)
all(map2_lgl(rules$dropped, links, \(d, k) identical(d, k$points$dropped)))
