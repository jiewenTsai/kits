# kits_item() 的三步，一步一個匯出的函數（讓學生一步一步跑；kits_algorithm(item) 列出它們）：
#   kits_stats()：充分統計量 -> kits_cml()：CML -> kits_thresholds()：門檻
# 和 kits_item() 是同一段程式（calibration_stats()、cml_estimate()、item_thresholds()），所以數字相同。
# 再加上兩個 theta 的函數（像 dexter 的 information()、expected_score()）：
#   kits_information()：測驗訊息；kits_expected_score()：期望原始分（test characteristic curve，tcc()）

#' The calibration step by step: counts, CML, thresholds
#'
#' [kits_item()] calibrates in three steps, and each step is also a function:
#' `kits_stats()` counts the sufficient statistics (`n_ig`, persons with a score of at least g on
#' item i; `n_br`, persons with raw score r in item pattern b), `kits_cml()` maximizes the
#' conditional likelihood of the PCM from these counts, and `kits_thresholds()` gives the threshold
#' of every step with its SE. They give the numbers of `item$stats`, `item$fit` and `coef(item)`.
#'
#' Details: `kits_algorithm("cml")` and tutorial_algorithms/01_calibration.html;
#' the models tutorial (tutorial_models/kits_models.html), sections 2, 3 and 5.
#'
#' @inheritParams kits_item
#' @param fixed_params Optional fixed step parameters (`item_id`, `item_score`, `beta`), as in
#'   [kits_item()]. Fixed items with no 0 score or one score are left out (`left_out`).
#' @param stats A `kits_stats()` result (or `item$stats` of a [kits_item()]).
#' @param fit A `kits_cml()` result (or `item$fit`).
#' @param x A `kits_stats`.
#' @param n Rows of each table to print (`Inf` for all).
#' @param ... Not used (an argument the function does not have stops with an error).
#' @return
#' * `kits_stats()`: a `kits_stats`, a list that prints briefly. Its elements:
#'   `items` (the items, in the order of `sort()`);
#'   `steps` (one row per step, `item_id`, `item_score` (the step g) and `n` (`n_ig`));
#'   `item patterns` (one row per item pattern, `pattern_id`, `items` (a list of item ids), `n`
#'   (persons) and `max_raw`);
#'   `raw` (one row per item pattern and raw score 0 to `max_raw`, `pattern_id`, `raw` and `n`
#'   (`n_br`));
#'   `fixed` (the fixed parameters of the items in the counts, or `NULL`) and `left_out` (the fixed
#'   parameters of the fixed items left out, or `NULL`);
#'   `n_persons` (the number of person-waves) and `info` (see [kits_info()]).
#' * `kits_cml()`: a `kits_fit`, with `coef` (`item_id`, `item_score`, `beta`, `SE_beta`), `est`
#'   (`beta`, `acov.beta`, `nll`), `inputs` (`design`, `method = "CML"`, `n_persons`, `left_out`)
#'   and `optim` (the convergence record).
#' * `kits_thresholds()`: a table with one row per step: `item_id`, `item_score`, `beta`, `beta_se`,
#'   `thr`, `thr_se`.
#' @seealso [kits_item()] (all three steps in one call), [kits_algorithm()].
#' @examples
#' betas <- list(a1 = c(-1, 0), a2 = c(0, 1), a3 = 0.5, a4 = -0.5)
#' X <- kits_sim(seq(-3, 3, length.out = 400), betas, seed = 1)
#' X[1:100, "a4"] <- NA   # two item patterns
#'
#' # 1 the sufficient statistics
#' stats <- kits_stats(X)
#' stats
#' stats$steps
#'
#' # 2 the CML
#' fit <- kits_cml(stats)
#' fit
#' coef(fit)
#'
#' # 3 the thresholds
#' thresholds <- kits_thresholds(fit)
#' thresholds
#'
#' # The same numbers are in the calibration
#' identical(stats$steps, kits_item(X)$stats$steps)
#' @name kits_item_steps
NULL

#' @rdname kits_item_steps
#' @export
kits_stats <- function(data, fixed_params = NULL, ...) {
  no_extra_args(..., fun = "kits_stats")
  call <- call_text(sys.call())
  data_expr <- substitute(data)
  resp <- resp_long(data, gaps = is.null(fixed_params))
  fixed <- if (!is.null(fixed_params)) param_frame(fixed_params)
  info <- make_info(call, list(fixed_params = fixed_text(fixed)), data_fingerprint(resp),
                    input = data_input(data, data_expr))
  calibration_stats(resp, fixed, info)
}

#' @rdname kits_item_steps
#' @export
print.kits_stats <- function(x, n = 6, ...) {
  need_rows(n)
  cat(sprintf(paste0("<kits_stats> sufficient statistics of the CML: %d person-waves, %d items, ",
                     "%d steps, %d item pattern(s)\n"),
              x$n_persons, length(x$items), nrow(x$steps), nrow(x$patterns)))
  cat("n_ig, persons with a score of at least g on item i ($steps):\n")
  print_rows(x$steps, 2, n, row.names = FALSE, more = "steps")
  cat("Item patterns ($patterns); n_br, persons per raw score r of each item pattern, ",
      "is in $raw:\n", sep = "")
  print_rows(x$patterns, 2, n, row.names = FALSE, more = "item patterns")
  if (!is.null(x$fixed) && nrow(x$fixed) > 0) {
    cat(sprintf("Fixed: %d step(s) of %d item(s)\n", nrow(x$fixed), n_distinct(x$fixed$item_id)))
  }
  if (!is.null(x$left_out)) cat("Fixed, left out (no 0 score or one score):", paste(unique(x$left_out$item_id), collapse = ", "),
                                "\n")
  cat("Next: kits_cml(x). The algorithm: kits_algorithm(x).\n")
  invisible(x)
}

#' @rdname kits_item_steps
#' @export
kits_cml <- function(stats) {
  abort_if(!inherits(stats, "kits_stats"),
           "`stats` must be a kits_stats() result, e.g. kits_cml(kits_stats(X)) or kits_cml(item$stats).")
  cml_estimate(stats)
}

#' @rdname kits_item_steps
#' @export
kits_thresholds <- function(fit) {
  abort_if(!inherits(fit, "kits_fit"),
           "`fit` must be a kits_cml() result (or item$fit); for a kits_item, use coef(item).")
  kits_table(item_thresholds(fit))
}


#' Test information and expected score as functions of theta
#'
#' Each returns a function of theta for a set of items:
#' `kits_information()` the test information I(theta) = sum_i Var(X_i | theta) under the PCM, and
#' `kits_expected_score()` the expected raw score sum_i E(X_i | theta) (the test characteristic
#' curve). theta is on the calibration scale of the step parameters (`beta`).
#'
#' Details: `kits_algorithm("information")` and tutorial_algorithms/06_scale.html;
#' the models tutorial (tutorial_models/kits_models.html), section 1.
#'
#' @param item The item parameters: a [kits_item()] result, its `fit` ([kits_cml()]), a [kits()], or
#'   a data.frame with `item_id`, `item_score`, `beta`.
#' @param items Optional item ids: the items of the test (default: all items).
#' @param pattern_id Optional: one item pattern of `item$design` (a [kits_item()] or a [kits_cml()]
#'   result); its items.
#' @param ... Not used: an argument the function does not have stops with an error (class
#'   `kits_error_input`).
#' @return A function of theta (a numeric vector), of class `kits_information` or
#'   `kits_expected_score`, that returns one value per theta. It keeps `attr(, "parms")`
#'   (the betas), `attr(, "items")` and `attr(, "scale")`.
#' @examples
#' betas <- list(a1 = c(-1, 0), a2 = c(0, 1), a3 = 0.5, a4 = -0.5)
#' X <- kits_sim(seq(-3, 3, length.out = 400), betas, seed = 1)
#' item <- kits_item(X)
#' info <- kits_information(item)
#' info
#' info(c(-2, 0, 2))
#' 1 / sqrt(info(0))   # the SE of theta at 0
#'
#' expected <- kits_expected_score(item, items = c("a1", "a2"))
#' expected(c(-2, 0, 2))
#' @export
kits_information <- function(item, ...) UseMethod("kits_information")

#' @rdname kits_information
#' @export
kits_information.default <- function(item, items = NULL, pattern_id = NULL, ...) {
  no_extra_args(..., fun = "kits_information")
  pcm_function(item, items, pattern_id, "information")
}

#' @rdname kits_information
#' @export
kits_information.kits_item <- kits_information.default

#' @rdname kits_information
#' @export
kits_information.kits <- function(item, ...) {
  kits_information(item$item, ...)   # kits() 的校準
}

#' @rdname kits_information
#' @export
kits_information.kits_fit <- kits_information.default

#' @rdname kits_information
#' @export
kits_expected_score <- function(item, ...) UseMethod("kits_expected_score")

#' @rdname kits_information
#' @export
kits_expected_score.default <- function(item, items = NULL, pattern_id = NULL, ...) {
  no_extra_args(..., fun = "kits_expected_score")
  pcm_function(item, items, pattern_id, "expected_score")
}

#' @rdname kits_information
#' @export
kits_expected_score.kits_item <- kits_expected_score.default

#' @rdname kits_information
#' @export
kits_expected_score.kits <- function(item, ...) {
  kits_expected_score(item$item, ...)   # kits() 的校準
}

#' @rdname kits_information
#' @export
kits_expected_score.kits_fit <- kits_expected_score.default

# 一組題的訊息或期望原始分，寫成 theta 的函數（tcc() 的 slope、tau）
pcm_function <- function(item, items, pattern_id, what) {
  parms <- param_frame(item)
  its <- function_items(item, parms, items, pattern_id)
  parms <- filter(parms, item_id %in% its)
  f <- theta_function(parms, its, if (what == "information") "slope" else "tau")
  structure(f, class = c(paste0("kits_", what), "function"), parms = parms, items = its, scale = beta_scale_of(item))
}

# 函數的環境裡只有 parms、its、part（不帶整個 kits_item）
theta_function <- function(parms, its, part) {
  force(parms)
  force(its)
  force(part)
  function(theta) {
    abort_if(!is.numeric(theta) || anyNA(theta), "`theta` must be numbers (no NA).")
    tcc(theta, parms, its)[[part]]
  }
}

# 一組題的期望總分（test characteristic curve，tau）和它的導數（slope）：PCM 時導數是各題分數變異數的和（測驗訊息）。
# 每題的 E、Var 來自 algo_pcm()（algo.R），從 0 照題的順序加。theta 可以是 RTMB 的 advector（information_slope()）
tcc <- function(theta, parms, items) {
  m <- map(betas_by_item(parms)[items], \(beta) algo_pcm(theta, beta))
  list(tau = Reduce(`+`, map(m, "E"), 0), slope = Reduce(`+`, map(m, "Var"), 0))
}

# 函數用的題：items、pattern_id 的題，或全部
function_items <- function(item, parms, items, pattern_id) {
  abort_if(!is.null(items) && !is.null(pattern_id), "Give `items` or `pattern_id`, not both.")
  if (!is.null(pattern_id)) {
    items <- booklet_items(item, pattern_id)
  }
  if (is.null(items)) {
    return(sort(unique(parms$item_id)))
  }
  need_params(items, parms$item_id, "Item(s) without parameters: ")
  sort(unique(items))
}

# 一個題本的題（kits_item 或 kits_fit 的 design）
booklet_items <- function(item, pattern_id) {
  if (inherits(item, "kits")) item <- item$item
  design <- if (inherits(item, "kits_item")) item$design else if (inherits(item, "kits_fit")) item$inputs$design
  abort_if(is.null(design),
           paste0("`pattern_id` needs a kits_item() or kits_cml() result, which knows its ", "item patterns."))
  abort_if(length(pattern_id) != 1 || !pattern_id %in% design$pattern_id,
           "`pattern_id` must be one of the item patterns: ",
           paste(utils::head(unique(design$pattern_id), 10), collapse = ", "), ".")
  design$item_id[design$pattern_id == pattern_id]
}

#' @rdname kits_information
#' @param x A `kits_information` or `kits_expected_score` function.
#' @param ... Not used.
#' @export
print.kits_information <- function(x, ...) {
  print_theta_function(x, "test information I(theta) = sum_i Var(X_i | theta)")
}

#' @rdname kits_information
#' @export
print.kits_expected_score <- function(x, ...) {
  print_theta_function(x, "expected raw score sum_i E(X_i | theta)")
}

print_theta_function <- function(x, what) {
  its <- attr(x, "items")
  cat(sprintf("<%s> %s, a function of theta\n", class(x)[1], what))
  cat(sprintf("  %d item(s): %s\n", length(its), some_ids(its, 8)))
  cat(sprintf("  theta on scale: %s\n", scale_label(attr(x, "scale"))))
  cat("  Call it with theta, e.g. f(c(-2, 0, 2)). The algorithm: kits_algorithm(f).\n")
  invisible(x)
}
