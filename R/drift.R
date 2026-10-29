#' Item drift: do the items behave the same in a new wave of data?
#'
#' Calibrates the new data freely ([kits_item()]), links its step thresholds to those of `ref`
#' ([kits_linking()]: a mean/sigma line with one round of screening), and flags the items with a
#' dropped step as drifted. Steps need at least `min_n` persons below and above to be judged.
#' Items that did not drift, or carry too little information to judge, are suggested for fixing
#' (`drift$fix`), ready for `kits_item(X_new, fixed_params = drift$fix)` on the scale of `ref`.
#' An item with a higher score in the new data than in `ref` (`"score range differs"`) is not
#' judged, since its thresholds depend on all its steps; it is suggested for fixing, and
#' [kits_item()] then keeps its old steps and estimates the new one.
#' Details: `vignette("kits-advanced", package = "kits")`, sections 8.1 and 8.2.
#'
#' @param data New data: a matrix (one row per person-wave, one column per item, row names are
#'   person ids), a wide data.frame of the same shape, or a long data.frame with `visit_id`,
#'   `item_id`, `item_score` (see [kits_item()]). Items with a score gap are not judged.
#' @param ref The [kits_item()] result of the earlier calibration, or its `thresholds` table, or
#'   a table of step parameters (`item_id`, `item_score`, `beta`; the thresholds are computed).
#' @param cut The cutoff for dropping, in residual SDs.
#' @param min_n The minimum number of persons needed both below and above a step.
#' @return A `kits_drift`: a list that prints briefly. Its elements:
#'   * `fix`: the `ref` parameters (`item_id`, `item_score`, `beta`) of the items suggested for
#'     fixing, with `attr(, "scale")`; NULL when `ref` has no `beta` column.
#'   * `intercept`, `slope`: the line from the free calibration of the new data to the `ref` scale.
#'   * `points`: one row per step: the thresholds, the numbers below and above, `judged`, the
#'     residual and whether it was dropped.
#'   * `items`: one row per item: `item_id`, `max_abs_resid`, `drift` (NA when not judged),
#'     `suggest_fix`, `reason` (why the item is not judged).
#'   * `steps_free`: the step table of the free calibration of the new data.
#'   * `scale`, `diagnostics` (see [kits_check()]), `info` (see [kits_info()]).
#' @examples
#' betas <- list(a1 = c(-1, 0), a2 = c(0, 1), a3 = 0.5, a4 = -0.5, a5 = c(-0.5, 0.5), a6 = 1.5)
#' old <- kits_sim(seq(-3, 3, length.out = 400), betas, seed = 1)
#' item_old <- kits_item(old)
#'
#' drifted <- modifyList(betas, list(a3 = 2))   # a3 became harder in the new data
#' new <- kits_sim(seq(-2, 4, length.out = 400), drifted, prefix = "n", seed = 2)
#' drift <- kits_drift(new, item_old)
#' drift$items
#'
#' # Calibrate the new data, fixing the items suggested for fixing
#' item_new <- kits_item(new, fixed_params = drift$fix)
#' coef(item_new)
#' @export
kits_drift <- function(data, ref, cut = 2, min_n = 30, ...) {
  no_extra_args(..., fun = "kits_drift")
  call <- call_text(sys.call())
  data_expr <- substitute(data)
  resp <- resp_long(data, gaps = FALSE)
  ref_scale <- beta_scale_of(ref)
  ref_source <- source_of(ref)
  if (is.data.frame(ref) && !"thr" %in% names(ref)) attr(ref, "rp") <- NULL   # beta 和 RP 無關
  ref <- with_thr(as_steps(ref, "ref"))
  need_steps(ref, "ref")
  ids <- sort(intersect(unique(resp$item_id), ref$item_id))
  abort_if(length(ids) == 0, "No item of the new data is in `ref`.", class = "kits_error_items")
  checks <- drift_item_checks(resp, ref)
  items <- left_join(data.frame(item_id = ids), checks, by = "item_id")   # 兩邊都有的題

  # 1 新資料自由校準：kits_item()，只用能校準的題
  steps_free <- drift_free(resp, checks)
  # 2 能判斷的題的共同步驟：新的門檻（thr_new）對 ref 的門檻（thr_ref）；上下各至少 min_n 人的步驟才用來判斷
  common <- steps_free |>
    filter(item_id %in% items$item_id[items$reason == ""]) |>
    inner_join(ref, by = c("item_id", "item_score"), suffix = c("_new", "_ref")) |>
    arrange(sort_rank(item_id), item_score) |>
    left_join(step_reach(resp, steps_free), by = c("item_id", "item_score")) |>
    mutate(judged = pmin(n_below, n_above) >= min_n)
  items <- mutate(items, reason = if_else(reason == "" & !item_id %in% common$item_id[common$judged],
                                          "too little information", reason))
  use <- filter(common, judged)
  abort_if(nrow(use) < 3,
           "Only ", nrow(use), " item step(s) have at least `min_n` = ", min_n, " persons below and above, ",
           "so drift cannot be judged. At least 3 are needed.", class = "kits_error_fit")
  # 3 判斷的步驟的直線，新資料 -> ref：kits_linking()（mean/sigma，篩一次）
  link <- kits_linking(select(use, item_id, item_score, thr = thr_new),
                       select(use, item_id, item_score, thr = thr_ref), cut = cut)
  # 4 每個步驟到最後那條線的殘差；判斷的步驟被篩掉就是偏移（沒判斷的步驟是 NA）
  points <- common |>
    transmute(item_id, item_score, n_below, n_above, judged, thr_new, thr_ref,
              resid = line_resid(link, thr_new, thr_ref),
              dropped = link$points$dropped[match(paste0(item_id, ":", item_score), link$points$id)])
  items <- drift_items(items, points)
  structure(list(fix = drift_fix(ref, items, ref_scale, ref_source), items = kits_table(items),
                 intercept = link$intercept, slope = link$slope, points = kits_table(points),
                 steps_free = steps_free, scale = ref_scale,
                 diagnostics = drift_diagnostics(items, nrow(use), cut),
                 info = make_info(call, list(cut = cut, min_n = min_n), data_fingerprint(resp),
                                  reference = list(scale = ref_scale, data = ref_source),
                                  input = data_input(data, data_expr))),
            class = "kits_drift")
}

# 新資料每題：能不能自由校準（有 0 分、至少兩種分數、分數沒有缺），兩邊都有的題不能判斷的原因（"" 是可以判斷）
drift_item_checks <- function(resp, ref) {
  resp |>
    summarise(max_new = max(item_score), can_fit = can_calibrate(item_score), gap = has_gap(item_score),
              .by = item_id) |>
    left_join(summarise(ref, max_ref = max(item_score), .by = item_id), by = "item_id") |>
    mutate(reason = case_when(max_new > max_ref ~ "score range differs",
                              !can_fit ~ "no 0 score or one score only",
                              gap ~ "score gap",
                              .default = ""))
}

# 新資料自由校準（kits_item()，只用能校準的題）；門檻表不帶量尺的標籤（自由校準的量尺沒有名字）
drift_free <- function(resp, checks) {
  fit_items <- filter(checks, can_fit, !gap)$item_id
  abort_if(length(fit_items) < 2,
           "Fewer than 2 items of the new data have a 0 score and at least two scores, so the new data cannot be calibrated.",
           class = "kits_error_fit")
  free <- kits_item(filter(resp, item_id %in% fit_items), scale = "new data")
  data.frame(free$thresholds)
}

# 每題一列：最大的殘差、有沒有偏移（不能判斷的題是 NA）、建議固定、原因
drift_items <- function(items, points) {
  by_item <- points |>
    summarise(drift = any(dropped, na.rm = TRUE), max_abs_resid = max(abs(resid)), .by = item_id)
  items |>
    left_join(by_item, by = "item_id") |>
    mutate(judged = reason == "",
           drift = if_else(judged, drift, NA),
           suggest_fix = (judged & !drift) |
             reason %in% c("too little information", "no 0 score or one score only", "score range differs")) |>
    select(item_id, max_abs_resid, drift, suggest_fix, reason)
}

# 建議固定的題在 ref 的參數（給 kits_item(fixed_params = )）；ref 沒有 beta 時 NULL
drift_fix <- function(ref, items, ref_scale, ref_source) {
  if (!"beta" %in% names(ref)) {
    return(NULL)
  }
  fix <- ref |>
    filter(item_id %in% items$item_id[items$suggest_fix]) |>
    select(item_id, item_score, beta) |>
    kits_table()
  structure(fix, scale = ref_scale, source = ref_source)
}

# n_steps：用來判斷的步驟數；cut 個殘差 SD 以外的步驟，就算沒有偏移也約有 2 pnorm(-cut)（cut = 2 時 4.6%）
drift_diagnostics <- function(items, n_steps, cut) {
  drifted <- items$item_id[which(items$drift)]
  not_judged <- sum(is.na(items$drift))
  by_chance <- 2 * stats::pnorm(-cut)
  chance <- sprintf("; with no drift at all, about %.1f of the %d judged steps are flagged by chance (%.1f%%)",
                    n_steps * by_chance, n_steps, 100 * by_chance)
  any_drift <- length(drifted) > 0
  drift_status <- if (any_drift) "note" else "ok"
  drift_items <- if (any_drift) paste0("not suggested for fixing: ", paste(drifted, collapse = ", ")) else "none"
  judged_message <- if (not_judged > 0) "see `reason` in $items" else "every common item was judged"
  diag_table(list(
    diag_row("drifted items", length(drifted), drift_status, paste0(drift_items, chance)),
    diag_row("items not judged", not_judged, "ok", judged_message)))
}


#' @rdname kits_drift
#' @param x A `kits_drift`.
#' @param digits Decimals to print (only the printout is rounded; the object keeps full precision).
#' @param ... Not used.
#' @export
print.kits_drift <- function(x, digits = 2, ...) {
  need_digits(digits)
  it <- x$items
  cat(sprintf("<kits_drift> new data -> ref (scale %s): intercept %s, slope %s\n",
              scale_label(x$scale), fmt(x$intercept, digits), fmt(x$slope, digits)))
  cat(sprintf("  items judged: %d of %d; drifted: %s\n", sum(!is.na(it$drift)), nrow(it),
              if (any(it$drift, na.rm = TRUE)) paste(it$item_id[which(it$drift)], collapse = ", ") else "none"))
  input <- x$info$input$name
  cat(sprintf("  suggested for fixing: %d item(s), in $fix (use kits_item(%s, fixed_params = x$fix))\n", sum(it$suggest_fix),
              if (is.null(input)) "<new data>" else input))
  print(format_table(it, digits), row.names = FALSE, right = TRUE)
  checks_line(x$diagnostics)
  invisible(x)
}

# 只有步驟參數（beta）的表：門檻從 beta 算（和 copied_thresholds() 相同）；每題要有步驟 1, 2, ..., m
with_thr <- function(ref) {
  if ("thr" %in% names(ref) || !all(c("item_id", "item_score", "beta") %in% names(ref))) {
    return(ref)
  }
  ref <- arrange(ref, item_id, item_score)
  bad <- filter(ref, any(item_score != seq_len(n())), .by = item_id)$item_id
  abort_if(length(bad) > 0, "`ref` must give the beta of steps 1, 2, ... of each item; item(s) ",
           "with other steps: ", some_ids(unique(bad), 10), ".")
  mutate(ref, thr = step_thresholds(beta)$thr, .by = item_id)
}
