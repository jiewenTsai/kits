# R 內建的泛型（stats）：confint()、vcov()、update()、logLik()、nobs()、simulate()。
# 使用者已經會的名字，不用記新的函數；AIC()、BIC() 跟著 logLik() 和 nobs() 就有。

#' Confidence intervals, covariances and other base methods for kits results
#'
#' The methods of R's own generics, so the usual names work:
#' * `confint()`: confidence intervals.
#'   For a [kits_item()], of each threshold: `thr` +/- z `thr_se` (normal; the SE is in
#'   `coef(item)`).
#'   For a link ([kits_linking()]), of the `intercept` and the `slope`: with bootstrap draws
#'   (`n_boot > 0`) the quantiles of the draws, otherwise none (NA). For a [kits_dif()], of each
#'   item's `dif`: `dif` +/- z SE, with the SE from the covariances of the two calibrations (`attr(,
#'   "se")`).
#' * `vcov()`: the covariance matrix of the step parameters `beta` of a [kits_item()] (rows and
#'   columns `item:score`; 0 for fixed steps).
#' * `update()`: fit again with some settings changed, for example `update(norm, degrees = c(4, 3,
#'   2, 1))` or `update(item, scale = "lan M03")`. The call recorded in [kits_info()] is run again
#'   in the calling environment, so the data must still be there under the same name (results keep
#'   no data); data given by value (recorded as `<data>`) cannot be refitted.
#' * `logLik()`, `nobs()`, and so `AIC()` and `BIC()`: for a [kits_norms()], the log-likelihood of
#'   the shash fit (summed over the groups of a norm by group), its number of coefficients and of
#'   person-waves; for a [kits_item()], the conditional log-likelihood of the CML, the number of
#'   free step parameters and of person-waves. Compare only fits of the same data (for example norms
#'   of different `degrees`).
#' * `simulate()`: responses from the step parameters of a [kits_item()] at given `theta`, as
#'   [kits_sim()] does;
#'   [kits_sim()] needs no calibration.
#'
#' @param object A kits result (see above).
#' @param parm Which rows: item ids (a [kits_item()] or a [kits_dif()]), or `"intercept"`, `"slope"`
#'   (a link);
#'   by default all.
#' @param level The confidence level.
#' @param nsim Number of simulated response matrices.
#' @param seed Random seed, as for [kits_sim()].
#' @param theta For `simulate()`: the abilities of the simulated person-waves, on the calibration
#'   scale.
#' @param ... For `update()`: the settings to change. Otherwise not used.
#' @return `confint()`: a data.frame with the estimate, `lower` and `upper`. `vcov()`: a matrix.
#'   `update()`: the new result. `logLik()`: a `logLik` object. `nobs()`: a number. `simulate()`: a
#'   response matrix (a list of them when `nsim > 1`).
#' @examples
#' betas <- list(a1 = c(-1, 0), a2 = c(0, 1), a3 = 0.5, a4 = -0.5)
#' X <- kits_sim(rnorm(400, 0, 1.5), betas, seed = 1)
#' item <- kits_item(X, scale = "lan M03")
#' confint(item)
#' vcov(item)[1:2, 1:2]
#' logLik(item)
#' AIC(item)
#' X_new <- simulate(item, theta = rnorm(5), seed = 2)
#' X_new
#'
#' set.seed(1)
#' age <- runif(600, 3, 36)
#' D <- 20 + 15 * log(age) + rnorm(600, 0, 3)
#' norm <- kits_norms(D, age, ages = 3:36, degrees = c(3, 2, 1, 1))
#' lower <- update(norm, degrees = c(2, 1, 1, 1))
#' AIC(norm, lower)
#' nobs(norm)
#' @name kits_methods
#' @importFrom stats coef confint vcov logLik nobs simulate update
NULL

# 常態近似的信賴區間的 z
level_z <- function(level) {
  need_number(level, "level", "one number between 0 and 1, e.g. 0.95", above = 0, below = 1)
  stats::qnorm(1 - (1 - level) / 2)
}

#' @rdname kits_methods
#' @export
confint.kits_item <- function(object, parm = NULL, level = 0.95, ...) {
  z <- level_z(level)
  thr <- as.data.frame(coef(object))
  if (!is.null(parm)) thr <- thr[thr$item_id %in% pick_items(unique(thr$item_id), parm), ]
  out <- transmute(thr, item_id, item_score, thr, thr_se, lower = thr - z * thr_se,
                   upper = thr + z * thr_se)
  kits_table(out)
}

#' @rdname kits_methods
#' @export
confint.kits_link <- function(object, parm = c("intercept", "slope"), level = 0.95, ...) {
  level_z(level)
  probs <- c((1 - level) / 2, 1 - (1 - level) / 2)
  draws <- object$draws
  ci <- function(term) {
    if (!link_has_draws(object)) {
      return(c(NA_real_, NA_real_))
    }
    unname(stats::quantile(draws[, term], probs))
  }
  terms <- intersect(c("intercept", "slope"), parm)
  rows <- map(terms, \(term) data.frame(term = term, estimate = object[[term]], lower = ci(term)[1],
                                         upper = ci(term)[2]))
  out <- list_rbind(rows)
  out$method <- if (link_has_draws(object)) "bootstrap quantiles" else "no draws (n_boot = 0)"
  kits_table(out)
}

#' @rdname kits_methods
#' @export
confint.kits_dif <- function(object, parm = NULL, level = 0.95, ...) {
  z <- level_z(level)
  se <- attr(object, "se")
  abort_if(is.null(se), "This DIF result has no item DIF (it needs two groups) or keeps no SE (an older kits).")
  d <- as.data.frame(object)
  if (!is.null(parm)) d <- d[d$item_id %in% pick_items(d$item_id, parm), ]
  out <- transmute(d, item_id, dif, se = unname(se[item_id]), lower = dif - z * se,
                   upper = dif + z * se)
  kits_table(out)
}

#' @rdname kits_methods
#' @export
vcov.kits_item <- function(object, ...) {
  cf <- object$fit$coef
  V <- object$fit$est$acov.beta
  dimnames(V) <- rep(list(paste0(cf$item_id, ":", cf$item_score)), 2)
  V
}

#' @rdname kits_methods
#' @export
logLik.kits_item <- function(object, ...) {
  V <- object$fit$est$acov.beta
  free <- qr(V)$rank   # 自由的步驟參數（CML 有一個限制，固定的步驟不算）
  structure(-object$fit$est$nll, df = free, nobs = object$n_persons, class = "logLik")
}

#' @rdname kits_methods
#' @export
nobs.kits_item <- function(object, ...) object$n_persons

# 常模的每一組：-log L、係數的個數、人次（一個常模時只有一組）
norm_fits <- function(object) {
  groups <- attr(object, "groups")
  parts <- if (is.null(groups)) list(object) else groups
  degrees <- attr(object, "kits_args")$degrees
  map(parts, \(g) list(value = g$neg_loglik[1], df = sum(degrees + 1), n = g$n[1])) |>
    map(as.data.frame) |>
    list_rbind()
}

#' @rdname kits_methods
#' @export
logLik.kits_norms <- function(object, ...) {
  fits <- norm_fits(object)
  structure(-sum(fits$value), df = sum(fits$df), nobs = sum(fits$n), class = "logLik")
}

#' @rdname kits_methods
#' @export
nobs.kits_norms <- function(object, ...) sum(norm_fits(object)$n)

#' @rdname kits_methods
#' @export
simulate.kits_item <- function(object, nsim = 1, seed = NULL, theta, ...) {
  if (missing(theta)) {
    kits_abort("Give the abilities of the persons to simulate: simulate(item, theta = ).",
               class = "kits_error_input")
  }
  parms <- param_frame(object)
  betas <- betas_by_item(parms)
  one <- function(i) kits_sim(theta, betas, prefix = if (nsim > 1) sprintf("sim%d_", i),
                              seed = seed_of(seed, i))
  if (nsim == 1) {
    return(one(1))
  }
  map(seq_len(nsim), one)
}

# 第 i 次模擬的種子
seed_of <- function(seed, i) if (is.null(seed)) NULL else seed + i - 1

#' @rdname kits_methods
#' @export
update.kits_norms <- function(object, ...) update_kits(object, ..., envir = parent.frame())

#' @rdname kits_methods
#' @export
update.kits_item <- function(object, ...) update_kits(object, ..., envir = parent.frame())

# 照 kits_info() 記下的呼叫，改幾個設定，在呼叫的環境裡再跑一次
update_kits <- function(object, ..., envir) {
  text <- info_of(object)$call
  abort_if(is.null(text) || grepl("<data>", text, fixed = TRUE),
           "The call of this result was not recorded with the names of its data, ",
           "so it cannot be refitted; call the function again.")
  call <- str2lang(text)
  changes <- list(...)
  abort_if(length(changes) > 0 && (is.null(names(changes)) || any(names(changes) == "")),
           "Give the settings to change by name, e.g. update(norm, degrees = c(4, 3, 2, 1)).")
  for (name in names(changes)) call[[name]] <- changes[[name]]
  eval(call, envir)
}
