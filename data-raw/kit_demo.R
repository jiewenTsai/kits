# Builds data/kitm03.rda and data/kitm36.rda: simulated KIT responses for teaching
# (language and motor, 500 children per cohort).
#
# No row comes from a real child. The script reads the real KIT data only to take aggregates:
# - the PCM item parameters of each cohort x domain calibration;
# - the standard form (most common item pattern) of each wave x domain;
# - the share of children seen at each wave, and the share who answered each item;
# - the mean and the true-score SD of theta at each wave;
# - the correlations between domains (same wave) and between adjacent waves, corrected for
#   measurement error;
# - the spread of age around each wave (percentiles), and the share of boys.
# The real cohort M03 stops at 48 months. The demo M03 cohort goes on to 60 and 72 months:
# there its children answer like the real M36 children of that age (M36 item parameters,
# booklets, wave means and SDs, age spread). Each child keeps its rank from its own earlier
# waves (the M36 correlation between adjacent waves). M03 item parameters would not do:
# in the real M03 data these items were answered only up to 48 months.
#
# Run from the kits/ folder with kits installed: Rscript data-raw/kit_demo.R
# Needs the real data in ../../kit_five_20261024 (not part of the package).

library(kits)
library(dplyr)
library(tidyr)
library(purrr)

N_CHILDREN <- 500
SEED <- 20261028
KIT <- "../../kit_five_20261024"
WAVES <- list(M03 = c(3, 6, 12, 18, 24, 36, 48, 60, 72), M36 = c(36, 48, 60, 72, 84, 96))
DOMAINS <- c("lan", "motor")

# ---- 1 Real data: only aggregates leave this section ----

old <- setwd(KIT)
kit <- new.env()
sys.source("R/functions.R", envir = kit)
resp <- suppressMessages(kit$read_responses())
setwd(old)

persons <- resp |>
  distinct(person_id, pid, cohort = sample, wave_m, age_m, sex)
# kits names (0.23.0): visit_id = KIT person_id, subj_id = KIT pid
kits_persons <- persons |>
  rename(visit_id = person_id, subj_id = pid)

calibrate <- function(dom, coh) {
  X <- resp |>
    filter(domain == dom, sample == coh) |>
    select(visit_id = person_id, item_id, item_score) |>
    pivot_wider(names_from = item_id, values_from = item_score, names_sort = TRUE)
  item <- suppressMessages(kits_item(X, scale = paste(dom, coh)))
  subj <- suppressMessages(kits_subj(X, item, persons = kits_persons)) |>
    as_tibble() |>
    rename(person_id = visit_id, pid = subj_id)
  list(item = item, subj = subj)
}
fits <- expand_grid(domain = DOMAINS, cohort = c("M03", "M36")) |>
  mutate(fit = map2(domain, cohort, calibrate))

# Item parameters (beta of each step), one table
betas <- fits |>
  mutate(tab = map(fit, \(f) as_tibble(f$item$thresholds))) |>
  select(domain, cohort, tab) |>
  unnest(tab) |>
  select(cohort, domain, item_id, item_score, beta)

# Standard booklet of each cohort x wave x domain, and the share who answered each of its items
sets <- resp |>
  group_by(sample, wave_m, domain, person_id) |>
  summarise(set = paste(sort(item_id), collapse = " "), .groups = "drop")
standard <- sets |>
  count(sample, wave_m, domain, set) |>
  slice_max(n, n = 1, by = c(sample, wave_m, domain), with_ties = FALSE) |>
  select(cohort = sample, wave_m, domain, set)
answered <- resp |>
  count(cohort = sample, wave_m, domain, item_id, name = "n_item") |>
  left_join(count(sets, cohort = sample, wave_m, domain, name = "n_domain"),
            by = c("cohort", "wave_m", "domain")) |>
  mutate(p_item = n_item / n_domain)
booklets <- standard |>
  mutate(item_id = strsplit(set, " ")) |>
  unnest(item_id) |>
  select(-set) |>
  left_join(select(answered, cohort, wave_m, domain, item_id, p_item),
            by = c("cohort", "wave_m", "domain", "item_id"))

# Attendance: share of the cohort's children seen at each wave, and in each domain
children <- persons |>
  summarise(n_children = n_distinct(pid), .by = cohort)
attend <- persons |>
  count(cohort, wave_m, name = "n_wave") |>
  left_join(children, by = "cohort") |>
  mutate(p_wave = n_wave / n_children)
attend_domain <- sets |>
  count(cohort = sample, wave_m, domain, name = "n_domain") |>
  left_join(select(attend, cohort, wave_m, n_wave), by = c("cohort", "wave_m")) |>
  mutate(p_domain = n_domain / n_wave)

# Wave means and true-score SDs of theta (WLE variance minus mean error variance)
scores <- fits |>
  mutate(subj = map(fit, "subj")) |>
  select(domain, subj) |>
  unnest(subj)
waves <- scores |>
  summarise(mean = mean(theta), var_obs = var(theta), var_err = mean(se^2),
            .by = c(cohort, wave_m, domain)) |>
  mutate(sd = sqrt(pmax(var_obs - var_err, 0.25 * var_obs)), rel = sd^2 / var_obs)

# Correlations corrected for measurement error, capped at 0.95
wide_by <- function(coh, cols) {
  scores |>
    filter(cohort == coh) |>
    select(pid, wave_m, domain, theta) |>
    pivot_wider(names_from = all_of(cols), values_from = theta)
}
true_cor <- function(x, y, rel_x, rel_y) {
  ok <- !is.na(x) & !is.na(y)
  min(cor(x[ok], y[ok]) / sqrt(rel_x * rel_y), 0.95)
}
rel_of <- function(coh, w, dom) {
  waves$rel[waves$cohort == coh & waves$wave_m == w & waves$domain == dom]
}
domain_cor <- function(coh) {
  by_wave <- wide_by(coh, "domain")
  pairs <- expand_grid(d1 = DOMAINS, d2 = DOMAINS) |> filter(d1 < d2)
  pairs |>
    mutate(r = map2_dbl(d1, d2, \(a, b) {
      ws <- intersect(waves$wave_m[waves$cohort == coh & waves$domain == a],
                      waves$wave_m[waves$cohort == coh & waves$domain == b])
      mean(map_dbl(ws, \(w) {
        rows <- by_wave$wave_m == w
        true_cor(by_wave[[a]][rows], by_wave[[b]][rows], rel_of(coh, w, a), rel_of(coh, w, b))
      }))
    }))
}
wave_cor <- function(coh) {
  ws <- sort(unique(waves$wave_m[waves$cohort == coh]))
  by_wave <- scores |>
    filter(cohort == coh) |>
    select(pid, wave_m, domain, theta) |>
    pivot_wider(names_from = wave_m, values_from = theta, names_prefix = "w")
  tibble(w1 = head(ws, -1), w2 = tail(ws, -1)) |>
    mutate(r = map2_dbl(w1, w2, \(a, b) {
      doms <- intersect(waves$domain[waves$cohort == coh & waves$wave_m == a],
                        waves$domain[waves$cohort == coh & waves$wave_m == b])
      mean(map_dbl(doms, \(d) {
        rows <- by_wave$domain == d
        true_cor(by_wave[[paste0("w", a)]][rows], by_wave[[paste0("w", b)]][rows],
                 rel_of(coh, a, d), rel_of(coh, b, d))
      }))
    }))
}
cor_domain <- map(c(M03 = "M03", M36 = "M36"), domain_cor)
cor_wave <- map(c(M03 = "M03", M36 = "M36"), wave_cor)

# Age around each wave: percentiles 1..99 of (age_m - wave_m); share of boys
age_spread <- persons |>
  mutate(dev = age_m - wave_m) |>
  reframe(p = 1:99 / 100, q = quantile(dev, 1:99 / 100, na.rm = TRUE), .by = c(cohort, wave_m))
p_boy <- persons |>
  distinct(cohort, pid, sex) |>
  summarise(p_boy = mean(sex == "boy"), .by = cohort)

rm(resp, persons, sets, scores, fits, kit)

# ---- 2 Extend M03 to 60 and 72 months with the M36 aggregates ----

late <- c(60, 72)
copy_late <- function(tab) {
  tab |>
    filter(cohort == "M36", wave_m %in% late) |>
    mutate(cohort = "M03") |>
    bind_rows(tab)
}
# On the M36 scale: these waves are simulated with the M36 item parameters (see source_of())
waves <- copy_late(waves)
booklets <- copy_late(booklets)
# M36 children also join at 84 months, so its share seen at 36-72 is low.
# M03 keeps its 48 month share.
attend <- attend |>
  filter(cohort == "M03", wave_m == 48) |>
  select(-wave_m) |>
  cross_join(tibble(wave_m = late)) |>
  bind_rows(attend)
attend_domain <- copy_late(attend_domain)
age_spread <- copy_late(age_spread)
cor_wave$M03 <- bind_rows(cor_wave$M03, filter(cor_wave$M36, w1 >= 48, w2 <= 72))

# ---- 3 Simulate ----

set.seed(SEED)

# Correlation of all domain x wave thetas of one child: domains (x) waves
wave_matrix <- function(coh) {
  ws <- WAVES[[coh]]
  adj <- cor_wave[[coh]]$r[match(head(ws, -1), cor_wave[[coh]]$w1)]
  R <- diag(length(ws))
  for (j in seq_along(ws)[-1]) {
    for (i in seq_len(j - 1)) R[i, j] <- R[j, i] <- prod(adj[i:(j - 1)])
  }
  R
}
domain_matrix <- function(coh) {
  C <- diag(length(DOMAINS))
  dimnames(C) <- list(DOMAINS, DOMAINS)
  tab <- cor_domain[[coh]]
  for (k in seq_len(nrow(tab))) C[tab$d1[k], tab$d2[k]] <- C[tab$d2[k], tab$d1[k]] <- tab$r[k]
  C
}

draw_children <- function(coh, first_id) {
  ws <- WAVES[[coh]]
  S <- kronecker(domain_matrix(coh), wave_matrix(coh))
  z <- matrix(rnorm(N_CHILDREN * nrow(S)), N_CHILDREN) %*% chol(S)
  key <- expand_grid(domain = DOMAINS, wave_m = ws)
  ids <- sprintf("K%04d", first_id + seq_len(N_CHILDREN) - 1)
  z_long <- tibble(pid = rep(ids, times = nrow(key)),
                   domain = rep(key$domain, each = N_CHILDREN),
                   wave_m = rep(key$wave_m, each = N_CHILDREN),
                   z = as.vector(z))
  boy <- runif(N_CHILDREN) < p_boy$p_boy[p_boy$cohort == coh]
  sex <- tibble(pid = ids, sex = if_else(boy, "boy", "girl"))
  z_long |>
    left_join(filter(waves, cohort == coh), by = c("domain", "wave_m")) |>
    filter(!is.na(mean)) |>
    mutate(theta = mean + sd * z, cohort = coh) |>
    left_join(sex, by = "pid") |>
    select(cohort, pid, sex, wave_m, domain, theta)
}

# Which child is seen at which wave and in which domain; age at the visit
visits <- function(kids) {
  seen <- kids |>
    distinct(cohort, pid, sex, wave_m) |>
    left_join(select(attend, cohort, wave_m, p_wave), by = c("cohort", "wave_m")) |>
    filter(runif(n()) < p_wave)
  ages <- seen |>
    mutate(u = runif(n())) |>
    rowwise() |>
    mutate(age_m = wave_m + {
      sp <- age_spread[age_spread$cohort == cohort & age_spread$wave_m == wave_m, ]
      approx(sp$p, sp$q, xout = min(max(u, 0.01), 0.99))$y
    }) |>
    ungroup() |>
    mutate(age_m = round(age_m, 2))
  kids |>
    inner_join(select(ages, cohort, pid, wave_m, age_m), by = c("cohort", "pid", "wave_m")) |>
    left_join(select(attend_domain, cohort, wave_m, domain, p_domain),
              by = c("cohort", "wave_m", "domain")) |>
    filter(!is.na(p_domain), runif(n()) < p_domain) |>
    select(-p_domain)
}

# Which calibration simulates a wave: its own, except M03 at 60 and 72 months (M36)
source_of <- function(cohort, wave_m) if_else(cohort == "M03" & wave_m %in% late, "M36", cohort)

# PCM responses to the standard booklet, then drop items at the real rate of not answering
respond <- function(seen) {
  cells <- seen |> nest(kids = c(pid, sex, age_m, theta), .by = c(cohort, wave_m, domain))
  pmap(cells, \(cohort, wave_m, domain, kids) {
    book <- booklets |>
      filter(.data$cohort == .env$cohort, .data$wave_m == .env$wave_m, .data$domain == .env$domain)
    b <- betas[betas$cohort == source_of(cohort, wave_m) & betas$domain == domain, ]
    item_betas <- map(set_names(book$item_id), \(i) b$beta[b$item_id == i])
    X <- kits_sim(kids$theta, item_betas, prefix = "x")
    as_tibble(X) |>
      mutate(pid = kids$pid, sex = kids$sex, age_m = kids$age_m, theta = kids$theta) |>
      pivot_longer(all_of(book$item_id), names_to = "item_id", values_to = "item_score") |>
      left_join(select(book, item_id, p_item), by = "item_id") |>
      filter(runif(n()) < p_item) |>
      transmute(pid, sample = cohort, sex, wave_m, age_m, item_id, item_score, domain, theta,
                source = source_of(cohort, wave_m))
  }) |>
    list_rbind()
}

# A rare score (a 0 on an easy item) can be missing by chance in a few hundred children;
# CML then cannot estimate the item (scores must start at 0 without gaps, over all waves of a
# cohort). Redraw only that item, for the same children, until its scores have no gap.
fill_scores <- function(rows, sample, domain, item_id) {
  beta_of <- function(src) {
    betas$beta[betas$cohort == src & betas$domain == domain & betas$item_id == item_id]
  }
  no_gap <- function(x) all(0:max(x) %in% x)
  tries <- 0
  while (!no_gap(rows$item_score) && tries < 1000) {
    rows <- rows |>
      mutate(item_score = as.vector(kits_sim(theta, list(i = beta_of(source[1])), prefix = "x")),
             .by = source)
    tries <- tries + 1
  }
  if (tries > 0) message(sample, " ", item_id, ": redrawn ", tries, " times")
  rows
}

responses <- list(M03 = 1, M36 = 1 + N_CHILDREN) |>
  imap(\(first, coh) draw_children(coh, first)) |>
  list_rbind() |>
  visits() |>
  respond()
kit_demo <- responses |>
  nest(rows = -c(sample, domain, item_id)) |>
  mutate(rows = pmap(list(rows, sample, domain, item_id), fill_scores)) |>
  unnest(rows) |>
  select(pid, sample, sex, wave_m, age_m, item_id, item_score, domain) |>
  mutate(person_id = sprintf("%s_%s_M%02d", sample, pid, wave_m), .before = 1) |>
  arrange(sample, pid, wave_m, domain, item_id) |>
  mutate(item_score = as.integer(item_score))

# Child ids: random 7-hex strings, like short git hashes, with at least one letter (so no
# reader turns them into numbers). They are drawn last with their own seed, so the responses
# are the same as with the working ids K0001, ...; they are not derived from any real id.
set.seed(SEED + 1)
work_ids <- unique(kit_demo$pid)
hex <- sprintf("%07x", sample.int(16^7, 3 * length(work_ids)) - 1L)
hex <- unique(hex[grepl("[a-f]", hex)])
stopifnot(length(hex) >= length(work_ids))
id_map <- tibble(pid = work_ids, hex_id = hex[seq_along(work_ids)])
kit_demo <- kit_demo |>
  left_join(id_map, by = "pid") |>
  mutate(pid = hex_id, person_id = sprintf("%s_%s_M%02d", sample, pid, wave_m)) |>
  select(-hex_id)

# One list per cohort: the responses, and the anchor key (gsed2406) for the items it has.
# kits names: visit_id (KIT person_id), subj_id (KIT pid), cohort (KIT sample).
cohort_list <- function(coh) {
  responses <- kit_demo |>
    filter(sample == coh) |>
    rename(visit_id = person_id, subj_id = pid, cohort = sample)
  anchors <- kits_anchors() |> filter(item_id %in% responses$item_id)
  list(responses = responses, anchors = anchors)
}
OUT <- Sys.getenv("KIT_DEMO_DIR", "data")
kitm03 <- cohort_list("M03")
kitm36 <- cohort_list("M36")
save(kitm03, file = file.path(OUT, "kitm03.rda"), compress = "xz")
save(kitm36, file = file.path(OUT, "kitm36.rda"), compress = "xz")

# ---- 4 Report (aggregates only) ----

kit_demo |> distinct(person_id, sample, wave_m) |> count(sample, wave_m) |> print(n = 20)
sizes <- file.size(file.path(OUT, c("kitm03.rda", "kitm36.rda")))
cat("rows:", nrow(kitm03$responses), nrow(kitm36$responses), " anchors:", nrow(kitm03$anchors),
    nrow(kitm36$anchors), " bytes:", sizes, "\n")

# Calibrate the demo data like the real data and compare (item betas, wave means)
demo_fit <- expand_grid(domain = DOMAINS, cohort = c("M03", "M36")) |>
  mutate(item = map2(domain, cohort, \(dom, coh) {
    X <- kit_demo |>
      filter(domain == dom, sample == coh) |>
      select(visit_id = person_id, item_id, item_score) |>
      pivot_wider(names_from = item_id, values_from = item_score, names_sort = TRUE)
    suppressMessages(kits_item(X))
  }))
demo_betas <- demo_fit |>
  mutate(tab = map(item, \(it) as_tibble(it$thresholds))) |>
  select(domain, cohort, tab) |>
  unnest(tab) |>
  select(cohort, domain, item_id, item_score, beta_demo = beta)
betas |>
  inner_join(demo_betas, by = c("cohort", "domain", "item_id", "item_score")) |>
  summarise(steps = n(), r = cor(beta, beta_demo), max_abs_diff = max(abs(beta - beta_demo)),
            .by = c(cohort, domain)) |>
  print()
