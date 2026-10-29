# Reference values for the validation tests of kits 0.23.0 (conditional item fit, item-restscore).
# kits does not depend on eRm, iarm or vcdExtra (not even Suggests): this script was run once,
# by hand, with those packages in a separate library, and the tests compare with the saved rds.
# It is not run by the tests (and is left out of the built package by .Rbuildignore).
#
# Made on 2026-10-28 with R 4.5, eRm 1.0.10, iarm 0.4.3, vcdExtra (the version iarm 0.4.3 loads),
# and kits 0.23.0 only for kits_sim() (the simulated data are saved too).
#
# Run from kits/tests/testthat/fixtures/fit_ref:
#   R_LIBS=<library with kits, eRm, iarm, vcdExtra> Rscript make_fit_ref.R

library(kits)
versions <- c(eRm = as.character(utils::packageVersion("eRm")),
              iarm = as.character(utils::packageVersion("iarm")),
              vcdExtra = as.character(utils::packageVersion("vcdExtra")))

# eRm's PCM thresholds as a kits parameter table (item_id, item_score, beta): the same PCM steps
erm_parms <- function(fit, items) {
  th <- eRm::thresholds(fit)$threshtable[[1]][, -1, drop = FALSE]
  do.call(rbind, lapply(seq_len(nrow(th)), function(i) {
    v <- as.numeric(stats::na.omit(th[i, ]))
    data.frame(item_id = items[i], item_score = seq_along(v), beta = v)
  }))
}

# iarm::item_restscore() returns rounded, formatted text: its own code, returning the unrounded
# numbers (observed and expected gamma, SE, p, BH-adjusted p)
iarm_restscore <- function(object) {
  f <- utils::getFromNamespace("item_restscore", "iarm")
  b <- as.list(body(f))
  text <- vapply(b, function(s) paste(deparse(s), collapse = ""), "")
  at_p <- which(grepl("^pkorr <-", text))
  out <- quote(return(data.frame(observed = mm[, 1], expected = expected, se = mm[, 2],
                                 p = pvalue, p_adj = pkorr)))
  body(f) <- as.call(c(b[1:at_p], out))
  suppressMessages(f(object))
}

# 1 Conditional infit and outfit: iarm::out_infit(eRm::PCM(Z)) on complete PCM data
b <- list(a1 = c(-1, 0), a2 = c(0, 1), a3 = 0.5, a4 = -0.5, a5 = c(-2, -1, 1), a6 = -1.5)
Z <- kits_sim(seq(-2.5, 2.5, length.out = 500), b, seed = 4)
fit <- eRm::PCM(as.data.frame(Z))
io <- suppressMessages(iarm::out_infit(fit))
saveRDS(list(data = Z, parms = erm_parms(fit, colnames(Z)), infit = unname(io$Infit),
             outfit = unname(io$Outfit), versions = versions,
             call = "iarm::out_infit(eRm::PCM(as.data.frame(data)))"),
        "cond_fit.rds")

# 2 Item-restscore, polytomous: iarm::item_restscore(eRm::PCM(X))
betas <- list(a1 = c(-1, 0.2), a2 = c(-0.5, 0, 1), a3 = 0.5, a4 = c(-0.3, 0.6),
              a5 = c(0, 0.4, 1.1), a6 = -0.7)
set.seed(3)
X <- kits_sim(rnorm(800, 0, 1.3), betas, seed = 11)
fit <- eRm::PCM(X)
saveRDS(list(data = X, parms = erm_parms(fit, colnames(X)), ref = iarm_restscore(fit),
             versions = versions, call = "iarm::item_restscore(eRm::PCM(data)), unrounded"),
        "restscore_pcm.rds")

# 3 Item-restscore, dichotomous: iarm::item_restscore(eRm::RM(X2)); eRm's item difficulties
b2 <- list(d1 = -1.2, d2 = -0.5, d3 = 0, d4 = 0.3, d5 = 0.8, d6 = 1.4)
set.seed(4)
X2 <- kits_sim(rnorm(600), b2, seed = 5)
fit <- eRm::RM(X2)
parms <- data.frame(item_id = colnames(X2), item_score = 1L, beta = -unname(stats::coef(fit)))
saveRDS(list(data = X2, parms = parms, ref = iarm_restscore(fit), versions = versions,
             call = "iarm::item_restscore(eRm::RM(data)), unrounded"),
        "restscore_rm.rds")

# 4 Goodman-Kruskal gamma and its asymptotic SE: vcdExtra::GKgamma()
set.seed(6)
tabs <- lapply(1:3, function(k) matrix(stats::rpois(4 * 7, 6), 7, 4))
gk <- lapply(tabs, vcdExtra::GKgamma)
saveRDS(list(tables = tabs, gamma = vapply(gk, function(g) g$gamma, 0),
             sigma = vapply(gk, function(g) g$sigma, 0), versions = versions,
             call = "vcdExtra::GKgamma(table)"),
        "gkgamma.rds")
