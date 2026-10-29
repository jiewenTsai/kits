# Kolen & Brennan (2014) 的第二批功能：連結的獨立檢查（同齡的兩群、合併校準）、篩點的規則（spread、unit）、
# 反應機率（RP）的門檻和通過月齡、品質管制的註記。每一項都和手算比。

set.seed(1)
betas <- list(a1 = -2, a2 = -1, a3 = c(-0.5, 0.5), a4 = 0, a5 = c(0.5, 1.5), a6 = 2, a7 = 3)
X03 <- kits_sim(rnorm(500, 0), betas[1:5], seed = 1)
X36 <- kits_sim(rnorm(500, 2), betas[3:7], prefix = "q", seed = 2)
item03 <- kits_item(X03, scale = "lan M03")
item36 <- kits_item(X36, scale = "lan M36")

# mean/sigma 和一次篩點（手算）
hand_link <- function(from, to, keep_rule) {
  b0 <- sd(to) / sd(from)
  a0 <- mean(to) - b0 * mean(from)
  e <- to - (a0 + b0 * from)
  drop <- keep_rule(e, b0)
  b <- sd(to[!drop]) / sd(from[!drop])
  list(intercept = mean(to[!drop]) - b * mean(from[!drop]), slope = b, drop = drop)
}

test_that("篩點：預設和 kitsPreview 的規則相同，並記下規則；mad 和 se 照手算", {
  from <- c(-2, -1, 0, 1, 2, 3, 0.5, 1.5, -0.5, 2.5)
  to <- 1 + 1.1 * from + c(0.05, -0.1, 0.08, -0.02, 0.1, -0.06, 0.03, 1.2, -0.04, 0.02)
  id <- paste0("p", 1:10)
  d <- link_of_points(from, to, id = id)
  h <- hand_link(from, to, function(e, b0) abs(e) > 2 * sd(e))
  expect_equal(c(d$intercept, d$slope), c(h$intercept, h$slope))
  expect_equal(d$points$dropped, h$drop)
  expect_equal(d$settings[c("spread", "unit")], list(spread = "sd", unit = "step"))
  expect_equal(d$points$z, d$points$resid / sd(d$points$resid))

  m <- link_of_points(from, to, id = id, spread = "mad", cut = 2.5)
  h <- hand_link(from, to, function(e, b0) abs(e) > 2.5 * 1.4826 * median(abs(e - median(e))))
  expect_equal(c(m$intercept, m$slope), c(h$intercept, h$slope))
  expect_equal(m$points$dropped, h$drop)

  se_f <- rep(c(0.05, 0.1), 5)
  se_t <- rep(0.08, 10)
  s <- link_of_points(from, to, id = id, spread = "se", se_from = se_f, se_to = se_t)
  h <- hand_link(from, to, function(e, b0) abs(e) / sqrt(se_t^2 + b0^2 * se_f^2) > 2)
  expect_equal(c(s$intercept, s$slope), c(h$intercept, h$slope))
  expect_equal(s$points$dropped, h$drop)
  expect_equal(s$points$se_from, se_f)
  # 只有 from 的 SE（錨題）：se_to 是 0
  s1 <- link_of_points(from, to, id = id, spread = "se", se_from = se_f)
  h <- hand_link(from, to, function(e, b0) abs(e) / (b0 * se_f) > 2)
  expect_equal(s1$points$dropped, h$drop)

  expect_error(link_of_points(from, to, id = id, spread = "se"), "se_from",
               class = "kits_error_input")
  expect_error(link_of_points(from, to, id = id, spread = "range"), "spread",
               class = "kits_error_input")
  expect_error(link_of_points(from, to, id = id, unit = "items"), "unit",
               class = "kits_error_input")
  expect_error(link_of_points(from, to, id = id, spread = "se", se_from = rep(0, 10)), "SE 0",
               class = "kits_error_input")
})

test_that("unit = \"item\"：一個步驟被標記時，整題的步驟都丟掉；重抽時也一樣", {
  from <- c(-2, -1.2, -1, 0, 0.4, 1, 2, 2.6, 3)
  to <- 0.5 + from
  to[5] <- to[5] + 1.5                     # i3 的第 2 步偏了
  item <- c("i1", "i2", "i3", "i3", "i3", "i4", "i5", "i5", "i6")
  id <- paste0(item, ":", c(1, 1, 1, 2, 3, 1, 1, 2, 1))
  step <- link_of_points(from, to, id = id, item = item)
  whole <- link_of_points(from, to, id = id, item = item, unit = "item")
  expect_equal(step$dropped, "i3:3")
  expect_equal(whole$dropped, c("i3:1", "i3:2", "i3:3"))
  expect_equal(whole$points$flagged, step$points$dropped)
  keep <- !item %in% "i3"
  expect_equal(whole$slope, sd(to[keep]) / sd(from[keep]))
  # 重抽：第一次重抽照手算重做（同樣的亂數）
  wb <- link_of_points(from, to, id = id, item = item, unit = "item", n_boot = 50)
  items <- unique(item)
  rows <- split(seq_along(item), factor(item, levels = items))
  set.seed(1)
  i <- unlist(rows[sample(items, replace = TRUE)], use.names = FALSE)
  first <- algo_screen(from[i], to[i], item[i], cut = 2, unit = "item") |>
    algo_mean_sigma()
  expect_identical(unname(wb$draws[1, ]), c(first$intercept, first$slope))
})

test_that("kits_linking() 的 spread = \"se\"：兩個校準用兩邊的 thr_se；錨題只用 from 的 thr_se", {
  co <- kits_linking(item36, item03, spread = "se")
  common <- merge(coef(item36), coef(item03), by = c("item_id", "item_score"), suffixes = c("_f", "_t"))
  h <- hand_link(common$thr_f, common$thr_t, function(e, b0) abs(e) / sqrt(common$thr_se_t^2 + b0^2 * common$thr_se_f^2) > 2)
  expect_equal(c(co$intercept, co$slope), c(h$intercept, h$slope))
  anchors <- data.frame(item_id = c("a1", "a2", "a3", "a4", "a5"), tau = c(40, 47, 51, 56, 60))
  dl <- kits_linking(item03, anchors, spread = "se", cut = 3)
  pass <- coef(item03)[coef(item03)$item_score == 1, ]
  pass <- pass[match(anchors$item_id, pass$item_id), ]
  h <- hand_link(pass$thr, anchors$tau, function(e, b0) abs(e) / (b0 * pass$thr_se) > 3)
  expect_equal(c(dl$intercept, dl$slope), c(h$intercept, h$slope))
  expect_equal(dl$points$se_to, rep(0, 5))
  expect_error(kits_linking(as.data.frame(coef(item36))[c("item_id", "item_score", "thr")], item03, spread = "se"),
               "thr_se", class = "kits_error_input")
})

test_that("連結照它的規則印出篩點；kits_algorithm(co) 說是哪個規則", {
  for (rule in list(c("mad", "item"), c("se", "item"), c("sd", "step"))) {
    co <- kits_linking(item36, item03, spread = rule[1], unit = rule[2], n_boot = 50, cut = 1.5)
    capture.output(algo <- kits_algorithm(co)$facts)
    expect_true(any(grepl(switch(rule[1], sd = "SD(e)", mad = "MAD(e)", se = "se_from"), algo, fixed = TRUE)))
    expect_equal(any(grepl("whole items", algo, fixed = TRUE)), rule[2] == "item")
    expect_output(print(co), "screening: ")
  }
  expect_s3_class(plot(kits_linking(item36, item03, spread = "mad")), "plotly")
  expect_s3_class(plot(kits_linking(item36, item03, spread = "se")), "plotly")
})

test_that("summary(link, theta = )：篩點前後的線在 theta 的值和差；有重抽時除以 se_link", {
  from <- c(-2, -1, 0, 1, 2, 3, 0.5, 1.5)
  to <- 1 + from
  to[8] <- to[8] + 3
  link <- link_of_points(from, to, id = paste0("p", 1:8), n_boot = 100, from_scale = "x",
                         to_scale = "y")
  th <- c(-1, 0, 2)
  sm <- summary(link, theta = th)
  sc <- sm$screening
  expect_equal(sc$before, link$intercept0 + link$slope0 * th)
  expect_equal(sc$after, link$intercept + link$slope * th)
  expect_equal(sc$shift, sc$after - sc$before)
  expect_equal(sc$se_link, link_se_of(link, theta = th, from = "x")$se_link)
  expect_equal(sc$z, sc$shift / sc$se_link)
  expect_equal(attr(sc, "scale"), "y")
  expect_output(print(sm), "Screening: the line before")
  expect_identical(sm$diagnostics, summary(link)$diagnostics)
  expect_null(summary(link)$screening)
  plain <- link_of_points(from, to, id = paste0("p", 1:8))
  expect_false("se_link" %in% names(summary(plain, theta = th)$screening))
  expect_error(summary(link, theta = "a"), "theta", class = "kits_error_input")
  expect_error(summary(link, theta = NA_real_), "theta", class = "kits_error_input")
})

test_that("領域的平衡：留下的點各領域幾個；一個領域超過 2/3 時註記（p. 239）；印出來也有", {
  set.seed(1)
  steps <- data.frame(item_id = paste0("a", 1:12), item_score = 1, thr = rnorm(12), domain = rep(c("motor", "lan", "cog"), c(9, 2, 1)))
  anchors <- data.frame(item_id = steps$item_id, tau = 60 + 7 * steps$thr + rnorm(12, 0, 0.3))
  dl <- kits_linking(steps, anchors)
  d <- dl$diagnostics
  row <- d[d$check == "domain balance of kept points", ]
  kept <- table(dl$points$domain[!dl$points$dropped])
  expect_equal(row$value, max(kept) / sum(kept))
  expect_equal(row$status, "note")
  expect_match(row$message, "p. 239")
  expect_output(print(dl), "kept points by domain: motor")
  steps$domain <- rep(c("motor", "lan", "cog"), each = 4)
  d2 <- kits_linking(steps, anchors)$diagnostics
  expect_equal(d2$status[d2$check == "domain balance of kept points"], "ok")
})

test_that("共同題佔各校準的比例：低於 20% 時註記（p. 288）", {
  co <- kits_linking(item36, item03)
  row <- co$diagnostics[co$diagnostics$check == "common items / items", ]
  expect_equal(co$items_total, c(from = 5, to = 5))
  expect_equal(row$value, 100 * 3 / 5)
  expect_equal(row$status, "ok")
  # 很多題、只有 3 題共同
  many <- c(betas[3:5], setNames(as.list(seq(-2, 2, length.out = 14)), paste0("b", 1:14)))
  Xm <- kits_sim(rnorm(800), many, prefix = "m", seed = 3)
  im <- kits_item(Xm, scale = "big")
  cm <- kits_linking(im, item03)
  row <- cm$diagnostics[cm$diagnostics$check == "common items / items", ]
  expect_equal(row$status, "note")
  expect_equal(row$value, 100 * 3 / 17)
  expect_null(kits_linking(item03, data.frame(item_id = c("a1", "a2", "a4"), tau = c(40, 47, 55)))$items_total)
})

test_that("RP：coef(item, rp = .5) 不變；其他 RP 的門檻是 P(X >= k) = rp 的 theta；二點計分時 = beta + logit(rp)，SE 相同", {
  expect_identical(coef(item03, rp = 0.5), coef(item03))
  t80 <- coef(item03, rp = 0.8)
  expect_equal(attr(t80, "rp"), 0.8)
  one <- t80[t80$item_id == "a1", ]
  expect_equal(one$thr, one$beta + qlogis(0.8), tolerance = 1e-8)
  expect_equal(one$thr_se, coef(item03)$thr_se[coef(item03)$item_id == "a1"], tolerance = 1e-6)
  # 多點計分：照 PCM 手算 P(X >= k | thr) = .8
  b <- t80$beta[t80$item_id == "a5"]
  for (k in 1:2) {
    t <- t80$thr[t80$item_id == "a5" & t80$item_score == k]
    eta <- c(0, cumsum(t - b))
    expect_equal(sum(exp(eta)[(k + 1):3]) / sum(exp(eta)), 0.8, tolerance = 1e-8)
  }
  expect_true(all(t80$thr > coef(item03)$thr))
  # 選列保留 RP；不同 RP 不能疊；連結只收 RP50
  expect_equal(attr(t80[1:2, ], "rp"), 0.8)
  expect_error(rbind(t80, coef(item03)), "RP", class = "kits_error_input")
  expect_equal(attr(rbind(t80[1:2, ], t80[3:4, ]), "rp"), 0.8)
  expect_error(kits_linking(item36, t80), "RP80", class = "kits_error_input")
  expect_error(kits_linking(coef(item03, rp = 0.8), data.frame(item_id = c("a1", "a2", "a4"), tau = c(40, 47, 55))),
               "RP50", class = "kits_error_input")
  # 換到 D：RP 跟著表走
  dl <- kits_linking(item03, data.frame(item_id = c("a1", "a2", "a3", "a4", "a5"), tau = c(40, 47, 51, 55, 60)))
  d80 <- kits_convert(t80, dl)
  expect_equal(d80$thr, dl$intercept + dl$slope * t80$thr)
  expect_equal(attr(d80, "rp"), 0.8)
  expect_error(coef(kits_convert(item03, dl), rp = 0.8), class = "kits_error_scale")
  expect_error(coef(item03, rp = 1), "rp", class = "kits_error_input")
})

test_that("age_pass()：p = .5 和以前的 kits_age50() 相同；其他 p 是 (logit(p) - b0) / b1；沒跨過 p 時 NA", {
  set.seed(4)
  age <- runif(900, 3, 30)
  wave <- round(age / 3) * 3
  P <- cbind(a1 = rbinom(900, 1, plogis((age - 10) / 3)), a2 = rbinom(900, 1, plogis((age - 18) / 2)),
             a3 = rbinom(900, 1, 0.97))
  a50 <- age50_of(P, age, wave)
  ap <- kits:::age_pass(P, age, wave, p = c(0.25, 0.5, 0.9))
  expect_equal(names(ap), c("item_id", "p25", "p50", "p90", "separated", "age_last_fail", "age_first_pass"))
  expect_equal(ap$item_id, c("a1", "a2", "a3"))
  expect_equal(stats::setNames(ap$p50, ap$item_id), a50)
  cf <- coef(glm(P[, "a1"] ~ age, family = binomial))
  expect_equal(unlist(ap[1, c("p25", "p50", "p90")], use.names = FALSE), unname((qlogis(c(0.25, 0.5, 0.9)) - cf[1]) / cf[2]))
  expect_true(all(is.na(unlist(ap[3, c("p25", "p50")]))))               # 每一波都超過 .5
  one <- kits:::age_pass(P[, "a1"], age, wave, p = c(0.5, 0.9))
  expect_equal(unlist(one[c("p50", "p90")], use.names = FALSE), unlist(ap[1, c("p50", "p90")], use.names = FALSE))
  expect_equal(age50_of(P[, "a1"], age, wave), unname(a50["a1"]))
  expect_equal(ap$separated, c(FALSE, FALSE, NA))
  expect_equal(ap$age_last_fail[1], max(age[P[, "a1"] == 0]))
  expect_equal(ap$age_first_pass[1], min(age[P[, "a1"] == 1]))
  # 題目表的 rp = 0.9：age90 是 p90
  norm <- data.frame(domain = NA_character_, age_m = 3:30, P10 = -4, P50 = 0, P90 = 4)
  thr <- data.frame(item_id = c("a1", "a2", "a3"), item_score = 1, thr = 0)
  attr(thr, "rp") <- 0.9
  data <- data.frame(visit_id = seq_len(900), age_m = age, wave_m = wave, P)
  expect_identical(kits_item_table(thr, norm, data = data, rp = 0.9)$age90, ap$p90)
})

test_that("kits_item_positions(rp = )：用 age80 欄、記下 RP；成長曲線和階梯的圖跟著標", {
  norm <- data.frame(score = "lan", age_m = 3:24)
  norm$P50 <- 20 + 15 * log(norm$age_m)
  norm$P10 <- norm$P50 - 5
  norm$P90 <- norm$P50 + 5
  items <- data.frame(item_id = c("a1", "a2"), score = "lan", column = "lan", D = c(45, 60), age80 = c(6, NA),
                      DIF = "none", stem = c("x", "y"))
  pos <- old_item_positions(items, norm, rp = 0.8)
  expect_equal(attr(pos, "rp"), 0.8)
  expect_equal(pos$x, c(6, pos$age_reach[2]))
  expect_error(old_item_positions(items, norm), "age50", class = "kits_error_input")
  p <- old_plot_growth(norm, pos, scores = "lan")
  ttl <- plotly::plotly_build(p)$x$layout$title
  expect_match(if (is.list(ttl)) ttl$text else ttl, "80% passing")
  expect_s3_class(old_plot_ladder(pos, columns = "lan", marks = data.frame(age_m = 12, D = 57)),
                  "plotly")
})

test_that("kits_item 的題本人數：主要的題本少於 400 人次時註記（p. 304）", {
  r <- kits:::booklet_size_row(c(390, 600, 5))
  expect_equal(r$status, "note")
  expect_equal(r$value, 390)
  expect_match(r$message, "p. 304")
  expect_equal(kits:::booklet_size_row(c(400, 600, 5))$status, "ok")
  expect_equal(kits:::booklet_size_row(rep(10, 30))$value, 10)     # 沒有 5% 的題本：看最大的
  d <- item03$diagnostics
  expect_equal(d$status[d$check == "smallest main item pattern"], "ok")
  small <- kits_item(X03[1:300, ])$diagnostics
  expect_equal(small$status[small$check == "smallest main item pattern"], "note")
})

test_that("kits_vertical()：共同波次的平均差除以連結誤差（z_link）、SD 比；超過時註記", {
  set.seed(2)
  s <- data.frame(cohort = rep(c("M03", "M36"), each = 300), wave_m = rep(c(24, 36, 36, 48), each = 150))
  s$theta <- 40 + 8 * log(s$wave_m) + rnorm(600, 0, 4) + ifelse(s$cohort == "M36", 2, 0)
  s$se <- 1.5
  s$se_link <- 0.5
  v <- kits_vertical(s)
  co <- attr(v, "cohorts")
  a <- s$theta[s$cohort == "M03" & s$wave_m == 36]
  b <- s$theta[s$cohort == "M36" & s$wave_m == 36]
  expect_equal(co$z_link, (mean(b) - mean(a)) / sqrt(0.5^2 + 0.5^2))
  expect_equal(co$sd_ratio, sd(b) / sd(a))
  d <- attr(v, "diagnostics")
  expect_equal(d$status[d$check == "cohorts at wave 36 vs linking error"], if (abs(co$z_link) > 2) "note" else "ok")
  s$theta[s$cohort == "M36" & s$wave_m == 36] <- mean(b) + (b - mean(b)) * 1.5
  d2 <- attr(kits_vertical(s), "diagnostics")
  expect_equal(d2$status[d2$check == "cohorts at wave 36, SD ratio"], "note")
  expect_null(attr(kits_vertical(s[names(s) != "se_link"]), "cohorts")$z_link)
  # 參照世代沒有 se_link（沒有換算）：當成 0
  s$se_link[s$cohort == "M03"] <- NA
  co2 <- attr(kits_vertical(s), "cohorts")
  m03 <- s$theta[s$cohort == "M03" & s$wave_m == 36]
  m36 <- s$theta[s$cohort == "M36" & s$wave_m == 36]
  expect_equal(co2$se_link, 0.5)
  expect_equal(co2$z_link, (mean(m36) - mean(m03)) / 0.5)
  expect_equal(co2$z, (mean(m36) - mean(m03)) / sqrt(var(m03) / 150 + var(m36) / 150 + 0.25))
})


# ---- kits_link_check() --------------------------------------------------------------------------------------------

lc_betas <- lapply(seq(-3, 3, length.out = 10), function(b) b + c(-0.5, 0.5))
names(lc_betas) <- paste0("a", 1:10)
set.seed(1)
lc_wave <- rep(c(24, 36, 36, 48), each = 400)
lc_theta <- 2 * log(lc_wave / 24) + rnorm(1600)
young <- kits_sim(lc_theta[1:800], lc_betas[1:7], prefix = "y", seed = 1)
old <- kits_sim(lc_theta[801:1600], lc_betas[4:10], prefix = "o", seed = 2)
lc_persons <- data.frame(visit_id = c(rownames(young), rownames(old)), wave_m = lc_wave)
item_young <- kits_item(young, scale = "young")
item_old <- kits_item(old, scale = "old")
lc_link <- kits_linking(item_old, item_young, n_boot = 100)
chk <- kits_link_check(lc_link, item_old, item_young, old, young, persons = lc_persons)

test_that("合併校準：疊起來的資料用 kits_item() 校準；各自到合併量尺的線（所有門檻的 mean/sigma）合成 from -> to", {
  both <- matrix(NA, 1600, 10, dimnames = list(c(rownames(young), rownames(old)), names(lc_betas)))
  both[rownames(young), colnames(young)] <- young
  both[rownames(old), colnames(old)] <- old
  conc <- coef(kits_item(both))
  line <- function(item) {
    m <- merge(coef(item), conc, by = c("item_id", "item_score"))
    slope <- sd(m$thr.y) / sd(m$thr.x)
    c(intercept = mean(m$thr.y) - slope * mean(m$thr.x), slope = slope)
  }
  f <- line(item_old)
  t <- line(item_young)
  row <- chk$methods[chk$methods$method == "concurrent CML", ]
  expect_equal(row$slope, f[["slope"]] / t[["slope"]], tolerance = 1e-8)
  expect_equal(row$intercept, (f[["intercept"]] - t[["intercept"]]) / t[["slope"]],
               tolerance = 1e-8)
  expect_equal(coef(chk$concurrent)$thr, conc$thr, tolerance = 1e-8)
  expect_equal(chk$lines$intercept, c(f[["intercept"]], t[["intercept"]]), tolerance = 1e-8)
})

test_that("同齡的兩群：常態能力分佈的邊際最大概似和暴力積分（integrate）相同；線照平均和 SD 對齊", {
  # 暴力法：一個波次的一個世代，對每個原始分用 integrate() 算邊際概似，再用 optim 找最大
  r <- kits:::resp_long(old[1:400, ])
  parms <- kits:::param_frame(item_old)
  raw <- tapply(r$item_score, r$visit_id, sum)
  cb <- lapply(split(parms$beta, parms$item_id), function(b) c(0, cumsum(b)))
  log_s <- function(th) sum(vapply(cb, function(c) log(sum(exp((seq_along(c) - 1) * th - c))), numeric(1)))
  nll <- function(par) {
    tab <- table(raw)
    -sum(vapply(names(tab), function(x) {
      f <- function(th) vapply(th, function(t) exp(as.numeric(x) * t - log_s(t)), numeric(1)) * dnorm(th, par[1], exp(par[2]))
      tab[[x]] * log(integrate(f, -Inf, Inf, rel.tol = 1e-10)$value)
    }, numeric(1)))
  }
  brute <- optim(c(0, 0), nll, method = "BFGS", control = list(reltol = 1e-12))
  ln <- kits:::latent_normal(r, parms)
  expect_equal(unname(ln[c("mean", "sd")]), c(brute$par[1], exp(brute$par[2])), tolerance = 1e-4)
  w <- chk$waves
  expect_equal(w$mean_from, unname(ln[["mean"]]), tolerance = 1e-8)
  sa <- chk$methods[chk$methods$method == "same-age groups, wave 36", ]
  expect_equal(sa$slope, w$sd_to / w$sd_from)
  expect_equal(sa$intercept, w$mean_to - sa$slope * w$mean_from)
})

test_that("在 theta 的比較：值、差、se_link、z、同齡的線的抽樣誤差（delta 法，和數值微分相同）", {
  at <- chk$at
  th <- unique(at$theta)
  s_old <- kits_subj(old, item_old, persons = lc_persons)
  expect_equal(th, unique(as.vector(tapply(s_old$theta, s_old$wave_m, median))))
  main <- at[at$method == "common items (the link)", ]
  expect_equal(main$value, lc_link$intercept + lc_link$slope * th)
  expect_equal(main$diff, rep(0, length(th)))
  expect_equal(main$se_link, link_se_of(lc_link, theta = th, from = lc_link$from_scale)$se_link)
  sa <- at[at$method == "same-age groups, wave 36", ]
  expect_equal(sa$z, sa$diff / sa$se_link)
  expect_equal(sa$z_both, sa$diff / sqrt(sa$se_link^2 + sa$se_method^2))
  # delta 法：對 4 個估計值做數值微分
  w <- chk$waves
  val <- function(p, t) {
    b <- p[4] / p[2]
    p[3] + b * (t - p[1])
  }
  p0 <- c(w$mean_from, w$sd_from, w$mean_to, w$sd_to)
  se <- c(w$se_mean_from, w$se_sd_from, w$se_mean_to, w$se_sd_to)
  for (k in seq_along(th)) {
    g <- vapply(1:4, function(j) {
      h <- replace(numeric(4), j, 1e-6)
      (val(p0 + h, th[k]) - val(p0 - h, th[k])) / 2e-6
    }, numeric(1))
    expect_equal(sa$se_method[k], sqrt(sum(g^2 * se^2)), tolerance = 1e-6)
  }
  expect_true(all(is.na(at$se_method[at$method == "concurrent CML"])))
})

test_that("kits_link_check()：沒有 persons 時跳過同齡的檢查；錯的輸入說清楚；print、summary、公式、kits_check", {
  nop <- kits_link_check(lc_link, item_old, item_young, old, young)
  expect_null(nop$waves)
  expect_equal(nop$methods$method, c("common items (the link)", "concurrent CML"))
  expect_equal(unique(nop$at$theta), unname(quantile(kits_subj(old, item_old)$theta, c(0.25, 0.5, 0.75))))
  expect_match(nop$diagnostics$message[nop$diagnostics$check == "same-age groups"], "skipped")
  # 沒有重抽：用 tol
  nob <- kits_link_check(kits_linking(item_old, item_young), item_old, item_young, old, young, theta = 0, tol = 1e-6)
  expect_false("z" %in% names(nob$at))
  expect_equal(nob$diagnostics$status[nob$diagnostics$check == "concurrent CML"], "note")
  expect_error(kits_link_check(lc_link, item_young, item_old, young, old), class = "kits_error_scale")
  anchors <- data.frame(item_id = c("a1", "a2", "a3"), tau = c(40, 45, 50))
  expect_error(kits_link_check(kits_linking(item_young, anchors), item_young, item_old, young, old), class = "kits_error_input")
  expect_error(kits_link_check(lc_link, item_old, item_young, old, young, tol = -1), class = "kits_error_input")
  expect_output(print(chk), "concurrent CML")
  expect_output(print(summary(chk)), "Diagnostics")
  expect_s3_class(kits_check(chk, all = TRUE), "kits_check")
  expect_equal(kits_info(chk)$settings$tol, 0.2)
})

test_that("同齡線的收斂用梯度和 Hessian 判斷（nlminb 的 singular convergence 不算沒收斂）", {
  betas <- lapply(seq(-3, 3, length.out = 10), function(b) b + c(-0.5, 0.5))
  names(betas) <- paste0("a", 1:10)
  set.seed(1)
  wave <- rep(c(24, 36, 36, 48), each = 400)
  theta <- 2 * log(wave / 24) + rnorm(1600)
  young <- kits_sim(theta[1:800], betas[1:7], prefix = "y", seed = 1)
  old <- kits_sim(theta[801:1600], betas[4:10], prefix = "o", seed = 2)
  persons <- data.frame(visit_id = c(rownames(young), rownames(old)), wave_m = wave)
  item_young <- kits_item(young, scale = "young")
  item_old <- kits_item(old, scale = "old")
  link <- kits_linking(item_old, item_young, n_boot = 20)
  chk <- kits_link_check(link, item_old, item_young, old, young, persons = persons)
  expect_true(all(chk$waves$converged))
  expect_false(any(grepl("convergence", kits_check(chk)$check)))
})
