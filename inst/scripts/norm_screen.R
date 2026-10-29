# The screening of the waves and the ages of the norm without kits: kits_norm_screen() and
# kits_norm_ages() (closed form; dplyr for the data).
# Parts: "screen" (1-2, the waves left out), "ages" (3, the ages of the norm).
# The names follow kits_algorithm().

library(dplyr)

# ==== inputs ====
# Replace with your own data. D: the score of each person-wave, wave_m: its wave (months), cohort:
# its cohort (or NA), age_m: its age in months. Here: the language WLE of kitm03.
lan <- filter(kits::kitm03$responses, domain == "lan")
wide <- tidyr::pivot_wider(lan, id_cols = visit_id, names_from = item_id, values_from = item_score)
X <- as.matrix(wide[-1])
rownames(X) <- wide$visit_id
persons <- distinct(lan, visit_id, cohort, wave_m, age_m)
scores <- kits::kits_subj(X, kits::kits_item(X, scale = "lan M03"), persons = persons)
D <- scores$theta
wave_m <- scores$wave_m
cohort <- scores$cohort
age_m <- scores$age_m

# ==== part: screen ====
# 1 The percentage of children at the lowest score of each (cohort x) wave, with D rounded to 3
#   decimals and the percentage to 1 decimal; a wave with cut % or more is left out
cut <- 20
waves <- tibble(D = round(D, 3), wave_m = wave_m, cohort = cohort) |>
  filter(is.finite(D), !is.na(wave_m)) |>
  summarise(n = n(), floor_pct = round(100 * mean(D == min(D)), 1), .by = c(cohort, wave_m)) |>
  mutate(normed = floor_pct < cut) |>
  arrange(cohort, wave_m)
waves

# 2 Each person-wave goes into the norm when its (cohort x) wave does
normed <- tibble(D = D, wave_m = wave_m, cohort = cohort) |>
  left_join(select(waves, cohort, wave_m, wave_normed = normed), by = c("cohort", "wave_m")) |>
  mutate(normed = is.finite(D) & wave_normed %in% TRUE) |>
  pull(normed)
table(normed)

# ==== part: ages ====
# 3 The ages of the norm: from the 1st percentile of age in the first wave to the 99th in the
#   last wave, both rounded outward to half a month
in_norm <- tibble(age_m = age_m, wave_m = wave_m)[normed, ] |>
  filter(is.finite(age_m))
first_wave <- filter(in_norm, wave_m == min(wave_m))$age_m
last_wave <- filter(in_norm, wave_m == max(wave_m))$age_m
lo <- floor(quantile(first_wave, 0.01, names = FALSE) / 0.5) * 0.5
hi <- ceiling(quantile(last_wave, 0.99, names = FALSE) / 0.5) * 0.5
ages <- round(seq(lo, hi, by = 0.5), 10)
range(ages)

# ==== compare ====
# The same waves, person-waves and ages as kits_norm_screen() and kits_norm_ages()
screen <- kits::kits_norm_screen(D, wave_m, cohort = cohort, cut = cut)
c(waves = identical(waves$floor_pct, screen$waves$floor_pct),
  normed = identical(normed, screen$normed),
  ages = identical(ages, kits::kits_norm_ages(age_m[normed], wave_m[normed])))
