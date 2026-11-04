# Engine-equivalence layer of the kits validation report.
#
# Rule: each in-scope numeric function of kits is compared ONCE, on FIXED data, against an established
# engine (dexter, psychotools, TAM, mirt; eRm, iarm, cNORM, kitsPreview are not installed).
# Data: only kits::kits_sim() with fixed seeds and the demo data kits::kitm03. No real data.
#
# Run:  source(system.file("validation", "engine.R", package = "kits"))
#       needs dexter, psychotools, TAM and mirt; writes engine_results.csv to KIT_VALIDATION_OUT
#
# Parametrisation notes (read from kits 0.26.0 R/ sources):
#   * kits beta = PCM step location (Masters 1982): eta_k = sum_{g<=k} (theta - beta_g).
#     dexter coef()$beta, psychotools threshpar(type = "mode") and TAM's PCM xsi use the same step
#     parametrisation; mirt's Rasch/gpcm uses intercepts d_k = -sum_{g<=k} beta_g with a = 1.
#   * Identification: kits CML fixes sum(all step betas) = 0 (cml_constraint()). dexter and psychotools
#     threshpar() use the same sum-zero origin; we still centre both sides on the mean of the steps
#     (a no-op for CML, an origin alignment for TAM's MML, whose origin is the person mean = 0).
#   * kits WLE = Jeffreys modal estimate (= Warm's WLE for the PCM); SE = 1 / sqrt(I(theta)).
#     TAM::tam.wle reports the same SE; dexter reports Warm's sqrt(I + (I'/2I)^2) / I (not compared).

t_start <- proc.time()
suppressPackageStartupMessages({
  library(kits)
  library(dexter)
  library(psychotools)
  library(TAM)
  library(mirt)
})

# Results go to KIT_VALIDATION_OUT (default: a folder in tempdir()); the installed package is not written to
out_dir <- Sys.getenv("KIT_VALIDATION_OUT", file.path(tempdir(), "kit_validation"))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

results <- list()
add <- function(fun, engine, ref_fun, data, quantity, a, b, tol, note = "", pass = NULL, rel = FALSE) {
  a <- as.numeric(a)
  b <- as.numeric(b)
  ok <- is.finite(a) & is.finite(b)
  d <- if (rel) abs(a[ok] - b[ok]) / pmax(abs(b[ok]), 1e-300) else abs(a[ok] - b[ok])
  mad <- if (any(ok)) max(d) else NA_real_
  r <- if (sum(ok) > 2 && stats::sd(a[ok]) > 0 && stats::sd(b[ok]) > 0) stats::cor(a[ok], b[ok]) else NA_real_
  if (is.null(pass)) pass <- isTRUE(mad <= tol)
  results[[length(results) + 1]] <<- data.frame(
    `function` = fun, reference_engine = engine, reference_function = ref_fun, data = data,
    quantity = quantity, n_values = sum(ok), max_abs_diff = signif(mad, 4), correlation = signif(r, 8),
    tolerance = tol, pass = pass, note = note, check.names = FALSE)
  invisible(NULL)
}
add_none <- function(fun, engine, note) {
  results[[length(results) + 1]] <<- data.frame(
    `function` = fun, reference_engine = engine, reference_function = NA_character_, data = NA_character_,
    quantity = NA_character_, n_values = 0L, max_abs_diff = NA_real_, correlation = NA_real_,
    tolerance = NA_real_, pass = NA, note = note, check.names = FALSE)
}
ver <- function(p) paste(p, as.character(utils::packageVersion(p)))

# ---- helpers ----------------------------------------------------------------------------------------
to_long <- function(X) {
  d <- data.frame(person_id = rep(rownames(X), ncol(X)), item_id = rep(colnames(X), each = nrow(X)),
                  item_score = as.vector(X), stringsAsFactors = FALSE)
  d[!is.na(d$item_score), ]
}
kits_coef <- function(item) {
  ck <- as.data.frame(coef(item))
  ck[order(ck$item_id, ck$item_score), ]
}
centre <- function(x) x - mean(x)
key <- function(d) paste(d$item_id, d$item_score)

# mirt model with every parameter fixed at kits' betas (Rasch / gpcm with a = 1, d_k = -cumsum(beta))
mirt_fixed <- function(X, parms) {
  X <- as.data.frame(X)
  vals <- mirt(X, 1, itemtype = "Rasch", pars = "values")
  bo <- split(parms$beta[order(parms$item_score)], parms$item_id[order(parms$item_score)])
  for (i in names(bo)) {
    d <- -cumsum(bo[[i]])
    nm <- if (length(d) == 1 && any(vals$item == i & vals$name == "d")) "d" else paste0("d", seq_along(d))
    for (k in seq_along(d)) vals$value[vals$item == i & vals$name == nm[k]] <- d[k]
  }
  vals$value[vals$name == "MEAN_1"] <- 0
  vals$value[vals$name == "COV_11"] <- 1
  vals$est <- FALSE
  mirt(X, 1, itemtype = "Rasch", pars = vals, TOL = NaN, verbose = FALSE)
}

# TAM PCM with xsi fixed at kits' betas (TAM's xsi order: by item column, then category)
tam_fixed <- function(X, parms) {
  X <- X[, sort(colnames(X)), drop = FALSE]
  tm0 <- tam.mml(X, irtmodel = "PCM", verbose = FALSE, control = list(maxiter = 1))
  nm <- rownames(tm0$xsi)                                  # "a1_Cat1", ... or "d1" for dichotomous
  k <- ifelse(grepl("_Cat", nm), paste(sub("_Cat.*$", "", nm), sub("^.*_Cat", "", nm)), paste(nm, 1))
  xsi <- parms$beta[match(k, key(parms))]
  stopifnot(!anyNA(xsi))
  list(model = tam.mml(X, irtmodel = "PCM", xsi.fixed = cbind(seq_along(xsi), xsi), verbose = FALSE), X = X)
}

# ---- fixed data -------------------------------------------------------------------------------------
b_dich <- stats::setNames(as.list(seq(-2, 2, length.out = 8)), paste0("d", 1:8))
b_poly <- list(a1 = c(-1, 0), a2 = c(0, 1), a3 = c(-2, -1, 1), a4 = c(-0.5, 0.5), a5 = c(0.5, 1.5, 2.5))
b_mix <- c(stats::setNames(as.list(c(-1.5, -0.5, 0.5, 1.5)), paste0("d", 1:4)),
           list(p1 = c(-1, 0), p2 = c(0, 1), p3 = c(-1, 0, 1), p4 = c(0.5, 1.5)))

set.seed(101); X_dich <- kits_sim(rnorm(1000, 0, 1.5), b_dich, seed = 101)
set.seed(102); X_poly <- kits_sim(rnorm(1000, 0, 1.5), b_poly, seed = 102)
set.seed(103); X_mix <- kits_sim(rnorm(1500, 0, 1.5), b_mix, seed = 103)
set.seed(104); X_inc <- X_mix
X_inc[sample(length(X_inc), round(0.35 * length(X_inc)))] <- NA     # incomplete design, many item patterns
X_inc <- X_inc[rowSums(!is.na(X_inc)) > 1, ]

# demo data: kitm03 language domain (booklet design by wave), long format
kitm03_lan <- kitm03$responses[kitm03$responses$domain == "lan", ]

D_dich <- sprintf("kits_sim seed 101: 1000 x 8 dichotomous, complete")
D_poly <- sprintf("kits_sim seed 102: 1000 x 5 polytomous (12 steps), complete")
D_inc <- sprintf("kits_sim seed 103/104: %d x 8 mixed, 35%% MCAR missing (incomplete design)", nrow(X_inc))
D_m03 <- sprintf("kitm03 demo, domain lan: %d responses, booklet design by wave", nrow(kitm03_lan))

it_dich <- kits_item(X_dich)
it_poly <- kits_item(X_poly)
it_inc <- kits_item(X_inc)
it_m03 <- kits_item(kits_resp(kitm03_lan), scale = "lan M03")

# ---- 1-2. kits_item (CML) vs dexter::fit_enorm (CML): betas and SEs ----------------------------------
cmp_dexter <- function(item, long, data_label) {
  f <- dexter::fit_enorm(long[, c("person_id", "item_id", "item_score")])
  cd <- coef(f)
  ck <- kits_coef(item)
  cd <- cd[match(key(ck), key(cd)), ]
  stopifnot(!anyNA(cd$beta))
  add("kits_item", ver("dexter"), "fit_enorm(method = 'CML') + coef()", data_label, "beta (step location)",
      centre(ck$beta), centre(cd$beta), 1e-5,
      "both sum-zero over steps; centred anyway (no-op). CML vs CML: expect agreement to optimiser tolerance")
  add("kits_item", ver("dexter"), "fit_enorm(method = 'CML') + coef()$SE_beta", data_label, "beta_se",
      ck$beta_se, cd$SE_beta, 1e-5, "same sum-zero constraint, so the SEs are comparable")
}
cmp_dexter(it_dich, to_long(X_dich), D_dich)
cmp_dexter(it_poly, to_long(X_poly), D_poly)
cmp_dexter(it_inc, to_long(X_inc), D_inc)
m03_long <- data.frame(person_id = kitm03_lan$visit_id, item_id = kitm03_lan$item_id,
                       item_score = kitm03_lan$item_score)
cmp_dexter(it_m03, m03_long[!is.na(m03_long$item_score), ], D_m03)

# 1b. kits_item vs psychotools::pcmodel (CML)
cmp_psy <- function(item, X, data_label) {
  pc <- pcmodel(X, hessian = TRUE)
  tp <- threshpar(pc, type = "mode")
  v <- unlist(tp)
  nm <- rownames(vcov(tp))          # "a1-C1", same order as unlist(tp)
  it_id <- sub("-C[0-9]+$", "", nm)
  sc <- as.integer(sub("^.*-C", "", nm))
  pd <- data.frame(item_id = it_id, item_score = sc, beta = unname(v),
                   se = sqrt(diag(vcov(tp))))
  ck <- kits_coef(item)
  pd <- pd[match(key(ck), key(pd)), ]
  stopifnot(!anyNA(pd$beta))
  add("kits_item", ver("psychotools"), "pcmodel() + threshpar(type = 'mode')", data_label,
      "beta (step location)", centre(ck$beta), centre(pd$beta), 1e-5,
      "threshpar() default origin = mean of thresholds 0 (same as kits)")
  add("kits_item", ver("psychotools"), "vcov(threshpar(pcmodel()))", data_label, "beta_se",
      ck$beta_se, pd$se, 1e-5, "sum-zero vcov on both sides")
}
cmp_psy(it_poly, X_poly, D_poly)
cmp_psy(it_inc, X_inc, D_inc)

# 1c. kits_item (CML) vs TAM::tam.mml (MML, normal ability distribution): close, not exact
tm_poly <- tam.mml(X_poly, irtmodel = "PCM", verbose = FALSE, control = list(conv = 1e-8, convD = 1e-8, maxiter = 2000))
ck <- kits_coef(it_poly)
xs <- tm_poly$xsi
xk <- paste(sub("_Cat.*$", "", rownames(xs)), sub("^.*_Cat", "", rownames(xs)))
add("kits_item", ver("TAM"), "tam.mml(irtmodel = 'PCM')$xsi", D_poly, "beta (step location)",
    centre(ck$beta), centre(xs$xsi[match(key(ck), xk)]), 0.15,
    paste("MML (N(mu, sigma) ability) vs CML: not exact by design. Origin aligned by centring both on the",
          "mean step (TAM fixes the person mean at 0). Pass = max diff <= 0.15 and r > 0.995"),
    pass = NULL)
results[[length(results)]]$pass <- results[[length(results)]]$max_abs_diff <= 0.15 &&
  results[[length(results)]]$correlation > 0.995

# ---- 3. kits_subj (WLE) vs TAM::tam.wle with item parameters fixed at kits' betas --------------------
cmp_wle <- function(item, X, data_label) {
  tf <- tam_fixed(X, kits_coef(item))
  w <- tam.wle(tf$model, progress = FALSE, convM = 1e-10, Msteps = 200)
  s <- kits_subj(X, item)
  w <- w[match(rownames(X), rownames(tf$X)), ]
  inner <- s$raw > 0 & s$raw < s$raw_max
  add("kits_subj", ver("TAM"), "tam.wle() on tam.mml(xsi.fixed = kits beta)", data_label,
      "theta (WLE), all persons incl. zero/perfect", s$theta, w$theta, 1e-6,
      sprintf("%d zero/perfect scores included; both finite WLEs for extremes", sum(!inner)))
  add("kits_subj", ver("TAM"), "tam.wle()$error", data_label, "se = 1/sqrt(I(theta))",
      s$se, w$error, 1e-6, "kits SE is 1/sqrt(I) (asymptotic WLE SE); TAM's error is the same definition")
}
cmp_wle(it_poly, X_poly, D_poly)
cmp_wle(it_inc, X_inc, D_inc)

# ---- 4. kits_itemfit ----------------------------------------------------------------------------------
# 4a. type = "wle": mean squares at each person's WLE; mirt::itemfit(fit_stats = "infit", Theta = WLE).
#     kits leaves out zero/perfect persons, so mirt gets the same persons.
s_poly <- kits_subj(X_poly, it_poly)
inner <- s_poly$raw > 0 & s_poly$raw < s_poly$raw_max
fw <- kits_itemfit(X_poly, it_poly, type = "wle")
mm_in <- mirt_fixed(X_poly[inner, ], kits_coef(it_poly))
mf <- itemfit(mm_in, fit_stats = "infit", Theta = matrix(s_poly$theta[inner]))
mf <- mf[match(fw$item_id, mf$item), ]
add("kits_itemfit", ver("mirt"), "itemfit(fit_stats = 'infit', Theta = kits WLE) on fixed-parameter model",
    D_poly, "infit and outfit (type = 'wle')", c(fw$infit, fw$outfit), c(mf$infit, mf$outfit), 1e-8,
    "same WLE theta and betas given to mirt; zero/perfect persons removed on both sides")
# 4b. type = "conditional" (default; Christensen-Kreiner-Mesbah, given the raw score): the engine is iarm,
#     not installed. kits ships values computed with iarm 0.4.3 + eRm 1.0.10 (tests/testthat/fixtures/fit_ref).
fx <- system.file("validation", "reference", package = "kits")   # copies of tests/testthat/fixtures/fit_ref
cf <- readRDS(file.path(fx, "cond_fit.rds"))
fc <- kits_itemfit(cf$data, cf$parms)
add("kits_itemfit", "iarm 0.4.3 + eRm 1.0.10 (stored values; NOT installed)", "out_infit(eRm::PCM()) via fixture cond_fit.rds",
    "kits_sim seed 4: 500 x 6 mixed (fixture data)", "infit and outfit (type = 'conditional'), eRm's betas",
    c(fc$infit, fc$outfit), c(cf$infit, cf$outfit), 1e-10,
    "not a live engine run: values stored in the kits repo; conditional fit has no installed engine")

# ---- 5. kits_information / kits_expected_score vs mirt from the same parameters -------------------
theta_grid <- seq(-5, 5, by = 0.25)
mm_poly <- mirt_fixed(X_poly, kits_coef(it_poly))
add("kits_information", ver("mirt"), "testinfo() on fixed-parameter model", D_poly,
    "test information I(theta), theta = -5..5", kits_information(it_poly)(theta_grid),
    testinfo(mm_poly, matrix(theta_grid)), 1e-10, "d_k = -cumsum(beta), a = 1")
add("kits_expected_score", ver("mirt"), "expected.test()", D_poly, "TCC E(raw | theta)",
    kits_expected_score(it_poly)(theta_grid), expected.test(mm_poly, matrix(theta_grid)), 1e-10, "")

# ---- 6. kits_pcm_prob vs mirt::probtrace -----------------------------------------------------------
pp <- kits_pcm_prob(it_poly, theta_grid)$prob
mp <- do.call(rbind, lapply(sort(colnames(X_poly)), function(i) {
  pt <- probtrace(extract.item(mm_poly, i), matrix(theta_grid))
  data.frame(item_id = i, theta = rep(theta_grid, each = ncol(pt)), x = rep(0:(ncol(pt) - 1), length(theta_grid)),
             p = as.vector(t(pt)))
}))
mg <- merge(as.data.frame(pp), mp, by = c("item_id", "theta", "x"))
stopifnot(nrow(mg) == nrow(pp))
add("kits_pcm_prob", ver("mirt"), "probtrace(extract.item())", D_poly, "P(X_i = x | theta)",
    mg$prob, mg$p, 1e-12, "")

# ---- kits_thresholds (Thurstonian RP50 thresholds) vs TAM::tam.threshold ----------------------------
tf_poly <- tam_fixed(X_poly, kits_coef(it_poly))
tt <- tam.threshold(tf_poly$model)
ck <- kits_coef(it_poly)
thr_tam <- tt[cbind(match(ck$item_id, rownames(tt)), ck$item_score)]
thr_k <- as.data.frame(kits_thresholds(it_poly$fit))
thr_k <- thr_k[match(key(ck), key(thr_k)), ]
add("kits_thresholds", ver("TAM"), "tam.threshold(prob.lvl = 0.5) on fixed-parameter model", D_poly,
    "Thurstonian threshold P(X >= k) = 0.5", thr_k$thr, thr_tam, 5e-4,
    "TAM finds thresholds by bisection to ~1e-4; kits by uniroot. Tolerance reflects TAM's grid")

# ---- 7. kits_esf vs psychotools::elementary_symmetric_functions --------------------------------------
bo <- split(ck$beta, ck$item_id)
e_k <- kits_esf(it_poly)
e_p <- elementary_symmetric_functions(lapply(bo, cumsum), order = 0)[[1]]
add("kits_esf", ver("psychotools"), "elementary_symmetric_functions(cumulative item-category pars)", D_poly,
    "gamma_r (relative difference)", e_k$gamma$gamma, e_p, 1e-10,
    "psychotools takes cumulative pars sum_{g<=x} beta_g; compared as relative difference", rel = TRUE)

# ---- kits_item_restscore vs iarm (stored values; iarm not installed) ---------------------------------
rr <- readRDS(file.path(fx, "restscore_pcm.rds"))
rk <- kits_item_restscore(rr$data, rr$parms)
add("kits_item_restscore", "iarm 0.4.3 + eRm 1.0.10 (stored values; NOT installed)",
    "item_restscore(eRm::PCM()) unrounded, fixture restscore_pcm.rds", "kits_sim seed 11: 800 x 6 polytomous (fixture data)",
    "gamma_obs, gamma_exp, se, p (eRm's betas)",
    c(rk$gamma_obs, rk$gamma_exp, rk$se, rk$p), c(rr$ref$observed, rr$ref$expected, rr$ref$se, rr$ref$p), 1e-10,
    "not a live engine run: values stored in the kits repo")

# ---- 8. kits_linking (mean/sigma with one screening round) vs hand formula ---------------------------
b_link <- list(a1 = -2, a2 = -1, a3 = c(-0.5, 0.5), a4 = 0, a5 = c(0.5, 1.5), a6 = 2, a7 = 3, a8 = c(-1, 1))
set.seed(201)
young <- kits_sim(rnorm(800, 0), b_link[c(1:5, 8)], seed = 201)
b_old <- modifyList(b_link, list(a4 = 1.5))      # a4 drifts in the older cohort: screening should drop it
old <- kits_sim(rnorm(800, 1.5), b_old[c(3:8)], prefix = "q", seed = 202)
iy <- kits_item(young, scale = "young")
io <- kits_item(old, scale = "old")
lk <- kits_linking(io, iy)
hand_screen_ms <- function(from, to, cut = 2) {
  ms <- function(f, t) { s <- stats::sd(t) / stats::sd(f); c(mean(t) - s * mean(f), s) }
  l0 <- ms(from, to)
  res <- to - (l0[1] + l0[2] * from)
  drop <- abs(res / stats::sd(res)) > cut
  c(ms(from[!drop], to[!drop]), sum(drop))
}
ty <- as.data.frame(iy$thresholds); to_ <- as.data.frame(io$thresholds)
common <- intersect(key(ty), key(to_))
h <- hand_screen_ms(to_$thr[match(common, key(to_))], ty$thr[match(common, key(ty))])
add("kits_linking", "hand formula (no installed engine: equateIRT/plink not installed)",
    "mean/sigma on common-step RP50 thresholds after one 2-SD residual screen",
    "kits_sim seeds 201/202: 2 cohorts x 800, 4 common items (6 steps), a4 drifted",
    "intercept, slope", c(lk$intercept, lk$slope), h[1:2], 1e-12,
    sprintf("slope = SD(to)/SD(from), intercept = mean(to) - slope*mean(from); hand drops %d point(s), kits drops %d",
            h[3], length(lk$dropped)))
# algo_mean_sigma directly (pure formula)
f0 <- c(-2, -1, 0, 1, 2, 3); t0 <- c(-1, 0, 1, 2, 3.2, 6)
am <- algo_mean_sigma(f0, t0)
add("algo_mean_sigma", "hand formula", "slope = sd(to)/sd(from); intercept = mean(to) - slope*mean(from)",
    "fixed vectors from ?algo_mean_sigma", "intercept, slope", c(am$intercept, am$slope),
    c(mean(t0) - sd(t0) / sd(f0) * mean(f0), sd(t0) / sd(f0)), 1e-14, "")

# ---- 9a. kits_drift: no reference engine; hand mean/sigma on its judged points -----------------------
b_dr <- list(a1 = c(-1, 0), a2 = c(0, 1), a3 = 0.5, a4 = -0.5, a5 = c(-0.5, 0.5), a6 = 1.5, a7 = -3.5)
old_d <- kits_sim(seq(-3, 3, length.out = 400), b_dr, seed = 301)
new_d <- kits_sim(seq(-2, 4, length.out = 400), modifyList(b_dr, list(a3 = 2)), prefix = "n", seed = 302)
dr <- kits_drift(new_d, kits_item(old_d))
pts <- as.data.frame(dr$points)
pj <- pts[pts$judged, ]
hd <- hand_screen_ms(pj$thr_new, pj$thr_ref)
add("kits_drift", "no reference engine (hand formula)", "mean/sigma + 2-SD screen on judged thresholds",
    "kits_sim seeds 301/302: 400 + 400, a3 drifted by +1.5", "intercept, slope",
    c(dr$intercept, dr$slope), hd[1:2], 1e-12,
    sprintf("no package implements kits' drift rule; only the linking arithmetic is checked. Flagged: %s (true: a3)",
            paste(dr$items$item_id[dr$items$drift %in% TRUE], collapse = ",")))

# ---- 9b. kits_dif: Andersen LR vs psychotools; overall DIF vs dexter::DIF -----------------------------
b_dif <- list(a1 = c(-1, 0), a2 = c(0, 1), a3 = 0.5, a4 = -0.5, a5 = 0, a6 = c(-0.5, 0.5, 1.5), a7 = 1)
gb <- modifyList(b_dif, list(a2 = c(-0.6, 0.4), a5 = 0.8))
set.seed(401)
X_dif <- rbind(kits_sim(rnorm(600), b_dif, prefix = "b", seed = 401),
               kits_sim(rnorm(600, 0.3), gb, prefix = "g", seed = 402))
sex <- ifelse(startsWith(rownames(X_dif), "b"), "boy", "girl")
dk <- kits_dif(X_dif, sex)
pa <- pcmodel(X_dif); pb <- pcmodel(X_dif[sex == "boy", ]); pg <- pcmodel(X_dif[sex == "girl", ])
LR_psy <- 2 * (as.numeric(logLik(pb)) + as.numeric(logLik(pg)) - as.numeric(logLik(pa)))
D_dif <- "kits_sim seeds 401/402: 600 boys + 600 girls, 7 items, DIF on a2, a5"
add("kits_dif", ver("psychotools"), "2*(logLik(pcmodel(boys)) + logLik(pcmodel(girls)) - logLik(pcmodel(all)))",
    D_dif, "Andersen LR statistic", attr(dk, "lrtest")$LR, LR_psy, 1e-5,
    "conditional log-likelihoods; complete data, all steps observed in both groups")
ld <- to_long(X_dif); ld$sex <- sex[match(ld$person_id, rownames(X_dif))]
dx <- dexter::DIF(ld, person_property = "sex")   # mirt also exports DIF()
add("kits_dif", ver("dexter"), "DIF(person_property = 'sex')$DIF_overall$stat", D_dif,
    "overall DIF chi-square (Bechger & Maris)", attr(dk, "overall")$stat, dx$DIF_overall$stat, 1e-5,
    sprintf("df kits %d vs dexter %d", attr(dk, "overall")$df, dx$DIF_overall$df))
# per-item dif = mean over steps of (beta_boy - beta_girl), each group's own CML on the joint origin;
# psychotools per-group threshpar() (sum-zero per group), aligned on the common origin by centring
tb <- unlist(threshpar(pb, type = "mode")); tg <- unlist(threshpar(pg, type = "mode"))
itn <- sub("-C[0-9]+$", "", rownames(vcov(threshpar(pb, type = "mode"))))
dpsy <- tapply(tb - tg, itn, mean)
dkit <- stats::setNames(dk$dif, dk$item_id)[names(dpsy)]
add("kits_dif", ver("psychotools"), "mean over steps of threshpar(boys) - threshpar(girls)", D_dif,
    "item dif, after removing the common origin shift", dkit - mean(dkit), dpsy - mean(dpsy), 1e-5,
    "per-item DIF is identified only up to a shift; compared after centring over items")

# ---- kits_precision vs brute-force enumeration with mirt probtrace + TAM WLE -------------------------
form <- sort(colnames(X_poly))
pr <- as.data.frame(kits_precision(it_poly, list(all = form), theta = seq(-3, 3, by = 0.5)))
grid <- as.matrix(expand.grid(lapply(bo[form], function(b) 0:length(b))))
colnames(grid) <- form
rownames(grid) <- paste0("pat", seq_len(nrow(grid)))
tfg <- tam_fixed(grid, ck)
wg <- tam.wle(tfg$model, progress = FALSE, convM = 1e-10, Msteps = 200)$theta
Pg <- sapply(pr$theta, function(t) {
  pi <- sapply(form, function(i) probtrace(extract.item(mm_poly, i), matrix(t))[1, grid[, i] + 1])
  apply(pi, 1, prod)
})
e_ref <- colSums(Pg * wg)
rmse_ref <- sqrt(colSums(Pg * outer(wg, pr$theta, "-")^2))
add("kits_precision", paste(ver("mirt"), "+", ver("TAM")), "full enumeration of 432 patterns: probtrace() x tam.wle()",
    D_poly, "E[WLE | theta] and RMSE, theta = -3..3", c(pr$expected, pr$rmse), c(e_ref, rmse_ref), 1e-6,
    "replaces kits' Lord-Wingersky recursion by explicit enumeration; WLEs from TAM")

# ---- 10. norms and descriptive functions: no installed engine ----------------------------------------
add_none("kits_norms / kits_norm_fit / kits_norm_z / kits_norm_table / algo_shash",
         "reference engine not installed (cNORM, gamlss)",
         "shash continuous norming; kits' own test-against-kitsPreview.R compares with cNORM via kitsPreview but is SKIPPED here (kitsPreview, cNORM not installed)")
add_none("kits_norm_screen / kits_norm_ages", "reference engine not installed (kitsPreview)",
         "wave screening rules are kits/kitsPreview-specific; no external engine")
add_none("kits_link_error / algo_link_error", "no reference engine",
         "bootstrap linking error (kits-specific resampling rule); no package implements the same scheme")
add_none("kits_vertical", "no reference engine",
         "descriptive growth statistics (means, SDs, Yen effect size); hand-formula check belongs to another layer")
add_none("kits_dif_combine", "reference engine not installed (kitsPreview)",
         "combination rule of two cohorts' DIF is kits-specific")

res <- do.call(rbind, results)
utils::write.csv(res, file.path(out_dir, "engine_results.csv"), row.names = FALSE)
el <- (proc.time() - t_start)[["elapsed"]]
cat(sprintf("kits %s; %d comparisons, %d pass, %d fail, %d without engine; %.1f s\n",
            utils::packageVersion("kits"), sum(!is.na(res$pass)), sum(res$pass %in% TRUE),
            sum(res$pass %in% FALSE), sum(is.na(res$pass)), el))
print(res[, c("function", "reference_engine", "quantity", "max_abs_diff", "correlation", "pass")], right = FALSE)
