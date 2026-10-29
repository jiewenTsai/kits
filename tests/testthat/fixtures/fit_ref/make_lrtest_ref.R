# Reference values for kits/tests/testthat/test-lrtest.R: eRm::LRtest() (Andersen LR test, PCM by CML).
# Run from the root of kits_20261027/ with eRm installed (eRm is not a dependency of kits):
#   Rscript make_lrtest_ref.R
# Writes kits/tests/testthat/fixtures/fit_ref/lrtest_ref.rds: the data and eRm's LR, df, p (and group sizes).
pkgload::load_all("kits", quiet = TRUE)   # only for kits_sim()
library(eRm)
lr_of <- function(e) list(LR = e$LR, df = e$df, p = e$pvalue, n = vapply(e$X.list, nrow, 1L))

# 1 Complete data: 1000 persons, 6 items with 2 or 3 categories
betas <- list(a1 = c(-1, 0), a2 = c(0, 1), a3 = 0.5, a4 = -0.5, a5 = c(-0.5, 0.8), a6 = 1)
X <- kits_sim(stats::qnorm(seq(0.0005, 0.9995, length.out = 1000)), betas, seed = 3)
set.seed(9)
split3 <- sample(rep(1:3, length.out = 1000))   # a user split vector: 3 random groups
median_ref <- lr_of(LRtest(PCM(X), splitcr = "median"))
user_ref <- lr_of(LRtest(PCM(X), splitcr = split3))

# (No eRm reference for an incomplete design: with NA, eRm 1.0.10's PCM()$loglik is not the conditional
# log-likelihood; on a two-booklet design it exceeded the CML maximum found by kits and checked by
# enumerating every response pattern, and LRtest() gave LR = 281 on data simulated from the model.)

# 2 An unobservable step: nobody in the low group has the top score 2 of a5
Xu <- X
rv <- rowSums(Xu)
low <- rv <= stats::median(rv)
Xu[low & Xu[, "a5"] == 2, "a5"] <- 1
rv_u <- rowSums(Xu)
low_u <- rv_u <= stats::median(rv_u)
stopifnot(!any(Xu[low_u, "a5"] == 2), any(Xu[!low_u, "a5"] == 2))
# unobserved = "drop": eRm leaves a5 out of both groups and of the joint model (with a warning)
drop_all_ref <- lr_of(suppressWarnings(LRtest(PCM(Xu), splitcr = "median")))
# unobserved = "steps": the low group estimates only step 1 of a5 (eRm's PCM() of a group uses the
# categories observed in it); the joint model uses every response
m_all <- PCM(Xu)
m_low <- PCM(Xu[low_u, ])
m_high <- PCM(Xu[!low_u, ])
steps_ref <- list(LR = 2 * (m_low$loglik + m_high$loglik - m_all$loglik),
                  df = length(m_low$etapar) + length(m_high$etapar) - length(m_all$etapar),
                  loglik = c(all = m_all$loglik, low = m_low$loglik, high = m_high$loglik))
steps_ref$p <- stats::pchisq(steps_ref$LR, steps_ref$df, lower.tail = FALSE)

ref <- list(X = X, split3 = split3, median = median_ref, user = user_ref,
            Xu = Xu, drop = drop_all_ref, steps = steps_ref,
            made_with = sprintf("eRm %s, %s", utils::packageVersion("eRm"), R.version.string))
dir.create("kits/tests/testthat/fixtures/fit_ref", showWarnings = FALSE, recursive = TRUE)
saveRDS(ref, "kits/tests/testthat/fixtures/fit_ref/lrtest_ref.rds", version = 3)
str(ref[c("median", "user", "drop", "steps", "made_with")])
