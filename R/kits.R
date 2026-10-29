#' Item calibration: the PCM by CML in RTMB, plus the threshold of every step
#'
#' Calibrates the items with the partial credit model by conditional maximum likelihood,
#' then computes the threshold `thr` of every step: the ability at which P(X >= k) = .5 (RP50).
#' The origin is the sum of all `beta` = 0, or the fixed items when `fixed_params` is given.
#' The calibration defines a scale, labelled by `scale`; the label travels with the results, so that
#' a wrong link stops.
#'
#' Details: `kits_algorithm("cml")` and tutorial_algorithms/01_calibration.html;
#' the models tutorial (tutorial_models/kits_models.html), sections 2, 3 and 5.
#'
#' @param data Responses, in one of three forms:
#'   * a matrix with one row per person-wave and one column per item: row names are person ids,
#'   column names are item ids, NA means not administered;
#'   * a wide data.frame of the same shape: the person ids are a column named `visit_id` (or else
#'   the row names), and every other column is an item;
#'   * a long data.frame with `visit_id`, `item_id`, `item_score` (recognized by the `item_id` and
#'   `item_score` columns);
#'   * a [kits_resp()] table (the responses with the columns of the child).
#'
#'   Scores are whole numbers from 0, with no gaps within an item (0 and 2 without 1: recode
#'   first); every free item needs a 0 and at least two scores. Include only the items of one
#'   domain. Rows without ids are numbered, empty rows and items are dropped (with a message).
#' @param ... Only `fixed_params`: items to fix at earlier estimates, for example
#'   `fixed_params = dplyr::filter(coef(old), item_id %in% fix)` (a whole `kits_item` is also
#'   accepted). Steps are fixed one by one: a fixed item with a higher score in `data` than in
#'   `fixed_params` keeps its fixed steps (their `beta`), and the new steps are estimated (a new
#'   wave where older children reach a score no younger child had); such an item needs no 0 score,
#'   but its scores above the fixed steps must be consecutive. The thresholds `thr` of the item
#'   change, as they depend on all its steps.
#' @param scale A label for the scale this calibration defines, for example `"lan M03"`. Always give
#'   it; by default the domain and cohort of a [kits_resp()], the name of `data`, or `"scale 1"`.
#'   With `fixed_params` that record a scale, that scale is used (another `scale` stops).
#' @return A `kits_item`: a list that prints briefly. Its elements:
#'   * `thresholds`: one row per step: `item_id`, `item_score` (the step), `beta`, `beta_se`, `thr`,
#'   `thr_se`. Also returned by `coef()`.
#'   * `fit`: the CML fit (class `kits_fit`, see [kits_cml()]).
#'   * `stats`: the sufficient statistics of the CML (class `kits_stats`, see [kits_stats()]).
#'   * `n_persons`: the number of person-waves in the calibration.
#'   * `left_out`: fixed items left out of the calibration (no 0 score, or one score only, and
#'   no score above their fixed steps).
#'   * `design`: the item patterns (`pattern_id`, `item_id`), including left-out items.
#'   * `scale`: the scale label (also `attr(coef(x), "scale")`).
#'   * `itemfit`: the conditional item fit on the calibration data, as [kits_itemfit()] (NULL when
#'   it cannot be computed).
#'   * `diagnostics`: the checks of the calibration (see [kits_check()]).
#'   * `info`: the call, settings and data fingerprint (see [kits_info()]).
#' @examples
#' betas <- list(a1 = c(-1, 0), a2 = c(0, 1), a3 = 0.5, a4 = -0.5)
#' X <- kits_sim(seq(-3, 3, length.out = 400), betas, seed = 1)
#' item <- kits_item(X, scale = "lan M03")
#' item
#' coef(item)
#'
#' # A new wave: fix the items of the first calibration (the new scale is that of the fixed items)
#' betas_new <- c(betas, list(a5 = c(0.5, 1.5)))
#' X_new <- kits_sim(seq(-2, 4, length.out = 400), betas_new, prefix = "n", seed = 2)
#' item_new <- kits_item(X_new, fixed_params = coef(item))
#' print(item_new, n = 8)
#'
#' # One wave at a time: old betas never move and earlier scores are never estimated again.
#' # In wave 2, a2 drifts (it gets a new version, a2_w2), a3 reaches a score of 2 (step 1 stays
#' # fixed, step 2 is estimated; the thresholds of a3 change) and a5 is new.
#' betas_1 <- c(betas, list(a6 = c(-1, 1), a7 = 1))
#' X_1 <- kits_sim(seq(-3, 3, length.out = 400), betas_1, seed = 1)
#' bank <- coef(kits_item(X_1, scale = "lan M03"))[c("item_id", "item_score", "beta", "beta_se")]
#' betas_2 <- c(betas_1[c("a1", "a4", "a6")],
#'              list(a2 = c(1, 2), a7 = 1, a3 = c(0.5, 2), a5 = c(1, 2)))
#' X_2 <- kits_sim(seq(-1, 5, length.out = 600), betas_2, prefix = "m", seed = 3)
#' drift <- kits_drift(X_2, bank)   # stop here when too many common items drift
#' moved <- drift$items$item_id[drift$items$drift %in% TRUE]
#' colnames(X_2)[colnames(X_2) %in% moved] <- paste0(moved, "_w2")
#' item_2 <- kits_item(X_2, fixed_params = drift$fix)
#' added <- dplyr::anti_join(coef(item_2), bank, by = c("item_id", "item_score"))
#' bank <- rbind(bank, added[names(bank)])
#' bank
#'
#' # Thresholds at RP80: the theta at which P(X >= k) = .8
#' coef(item, rp = 0.8)
#' @export
kits_item <- function(data, ..., scale = NULL) {
  check_dots(match.call(expand.dots = FALSE)$..., "kits_item")
  data_expr <- substitute(data)
  resp <- resp_long(data, gaps = is.null(list(...)$fixed_params))   # fixed items may skip a score
  # 預設的標籤只在沒有 scale、固定參數也沒有記錄量尺時才算（default 是 promise）
  calibrate(resp, list(...), scale, default_scale(data_expr, data), call_text(sys.call()),
            data_input(data, data_expr))
}

# kits_item() 不收 predicate（沒有名字的參數）
check_dots <- function(dots, fun) {
  nms <- names(dots)
  abort_if(length(dots) > 0 && (is.null(nms) || any(nms %in% c("", "predicate"))),
           fun, "() takes no predicate; give other arguments by name (e.g. fixed_params = ...). ",
           "To use part of the data, select rows and columns of the matrix first, e.g. X[keep, items].")
}

# 校準（給 kits_item()、kits()）。resp：長格式；args：fixed_params、method；default：預設的量尺標籤
calibrate <- function(resp, args, scale, default, call, input = NULL) {
  check_calibrate_args(args)
  if (!is.null(scale)) need_scale_arg(scale)
  ref <- fixed_reference(args$fixed_params, scale)
  scale <- if (is.null(ref$scale)) default else ref$scale
  data <- data_fingerprint(resp)
  info <- make_info(call, list(fixed_params = fixed_text(ref$fixed), method = "CML", scale = scale), data,
                    reference = ref$reference, input = input)

  stats <- calibration_stats(resp, ref$fixed, info)   # 1 充分統計量 n_ig、n_br（item$stats）
  fit <- cml_estimate(stats)                          # 2 CML（item$fit）
  thr <- item_thresholds(fit)                         # 3 門檻和 SE（item$thresholds）
  item_result(resp, stats, fit, thr, scale, info)
}

# 充分統計量（給 kits_item()、kits_stats()）。CML 要每一題都有 0 分、至少兩種分數：
# 固定的題不符合時，不放進校準（有訊息），參數照抄固定的值；不是固定的題停下來。
calibration_stats <- function(resp, fixed, info = NULL) {
  bad <- left_out_items(resp, fixed)
  calib <- filter(resp, !item_id %in% bad)
  check_consecutive(calib, fixed)
  stats <- cml_stats(calib, fixed, left_out = bad)
  stats$info <- info
  stats
}

# kits_item 的物件：門檻表記下量尺和資料，加上題本的設計、題目適配、診斷
item_result <- function(resp, stats, fit, thr, scale, info) {
  bad <- if (is.null(stats$left_out)) character() else unique(stats$left_out$item_id)
  calib <- filter(resp, !item_id %in% bad)
  books <- calibration_booklets(resp)
  thr <- kits_table(thr)
  attr(thr, "scale") <- scale
  attr(thr, "source") <- info$data
  item <- structure(list(thresholds = thr, fit = fit, stats = stats, n_persons = stats$n_persons, left_out = bad,
                         design = booklet_design(books), scale = scale, itemfit = NULL, diagnostics = NULL,
                         info = info),
                    class = "kits_item")
  itemfit <- item_itemfit(resp, item)
  item["itemfit"] <- list(itemfit$table)   # NULL 時也留著這個元素
  item$diagnostics <- item_diagnostics(fit, calib, bad, stats$n_persons, books$n, itemfit_row(itemfit))
  item
}

# 題目適配（item$itemfit）：和 kits_itemfit(data, item) 同一個表（同樣的程式），只有 info 的呼叫是這次校準的。
# 只有一題一列的彙總，沒有個人的值，所以 item 還是可以分享。
# 算不出來時（例如少於 2 題）table 是 NULL，reason 說原因；不會讓校準停下來
item_itemfit <- function(resp, item) {
  tryCatch(list(table = itemfit_table(resp, item, item$info$call, item$info$data)),
           error = function(e) list(table = NULL, reason = conditionMessage(e)))
}

# kits_item() 的 ... 只收 fixed_params 和 method = "CML"
check_calibrate_args <- function(args) {
  unknown <- setdiff(names(args), c("fixed_params", "method"))
  abort_if(length(unknown) > 0,
           "kits_item() in kits accepts only `fixed_params` (and `method = \"CML\"`); unknown argument(s): ",
           paste(unknown, collapse = ", "), ".")
  abort_if(!is.null(args$method) && !identical(args$method, "CML"),
           "Only CML calibrations are supported (the default method = \"CML\" of kits_item()).")
}

# 固定參數：步驟參數的表（fixed）、來源（reference），和量尺（固定參數記錄的量尺定新的校準的量尺；沒有時是 scale）
fixed_reference <- function(fixed_params, scale) {
  if (is.null(fixed_params)) {
    return(list(fixed = NULL, reference = NULL, scale = scale))
  }
  fixed <- param_frame(fixed_params)
  fixed_scale <- beta_scale_of(fixed_params)
  abort_if(!is.na(fixed_scale) && !is.null(scale) && scale != fixed_scale,
           sprintf("`fixed_params` are on scale %s, so the new calibration is on that scale too, not on %s. ", fixed_scale,
                   scale), "Leave `scale` out, or fix parameters of the scale you mean.", class = "kits_error_scale")
  list(fixed = fixed,
       reference = list(scale = fixed_scale, data = source_of(fixed_params)),
       scale = if (is.na(fixed_scale)) scale else fixed_scale)
}

# 沒有 0 分或只有一種分數的題（照題號排）：是固定的題時不放進校準（有訊息），不是時停下來。
# 固定題這次有比固定的步驟更高的分數時留在校準裡：低的步驟固定，新的步驟沒有 0 分也估得出來
# （估不出來時，例如沒有人在新的步驟以下，check_cml_existence() 停下來）
left_out_items <- function(resp, fixed) {
  scores <- resp |>
    summarise(low = min(item_score), top = max(item_score), ok = can_calibrate(item_score),
              .by = item_id) |>
    mutate(ok = ok | top > coalesce(fixed_top(fixed, item_id), Inf))
  bad <- sort(filter(scores, !ok)$item_id)
  not_fixed <- setdiff(bad, fixed$item_id)
  if (length(not_fixed) > 0) {
    from_one <- filter(scores, item_id %in% not_fixed, low > 0)$item_id
    one_score <- setdiff(not_fixed, from_one)
    kits_abort(score_problem_text(from_one, one_score, all(scores$low > 0), !is.null(fixed)), class = "kits_error_items")
  }
  if (length(bad) > 0) {
    kits_inform("Fixed item(s) with no 0 score or only one score in these data are left out of the calibration ",
                "and keep their fixed parameters: ", paste(bad, collapse = ", "), ".", class = "kits_message_fit")
  }
  abort_if(all(unique(resp$item_id) %in% bad),
           "Every item is fixed and has no 0 score or only one score, so there is nothing to calibrate. ",
           "Score with the fixed parameters directly: kits_subj(X, fixed).", class = "kits_error_fit")
  bad
}

# 沒有 0 分（多半是分數從 1 開始）和只有一種分數的題，各說該怎麼做
score_problem_text <- function(from_one, one_score, all_from_one, has_fixed) {
  fix_it <- if (has_fixed) " or fix it" else ""
  fix_them <- if (has_fixed) " or fix them" else ""
  parts <- character()
  if (all_from_one) {
    parts <- c(parts, paste0("No item has a score of 0: the scores must start at 0. If they are coded 1, 2, 3, ..., ",
                             "subtract 1 so that scores start at 0 (for example X - 1)."))
  } else if (length(from_one) > 0) {
    parts <- c(parts, paste0("Item(s) with no 0 score: ", paste(from_one, collapse = ", "), ". Scores must start at 0: ",
                             "if an item is coded 1, 2, 3, ..., subtract 1 so that its scores start at 0; ",
                             "if nobody scored 0 on it, CML cannot estimate it: leave it out", fix_it, "."))
  }
  if (length(one_score) > 0) {
    parts <- c(parts, paste0("Item(s) with only one score: ", paste(one_score, collapse = ", "),
                             ". CML cannot estimate them; leave them out of the data", fix_them, "."))
  }
  paste(parts, collapse = " ")
}

# 分數要連續（0, 1, 2, ...）：自由估的題從 0 分（都有 0 分）、固定題從固定的最高步驟 low 往上看
# （新的步驟；固定的步驟之間的缺沒關係）。第一個有缺的題（照題號排）停下來
check_consecutive <- function(calib, fixed = NULL) {
  gaps <- calib |>
    distinct(item_id, item_score) |>
    mutate(low = coalesce(fixed_top(fixed, item_id), 0L)) |>
    arrange(sort_rank(item_id), item_score) |>
    summarise(low = low[1], gap = has_gap(c(low, item_score[item_score > low])),
              scores = list(item_score[item_score > low]), .by = item_id) |>
    filter(gap)
  if (nrow(gaps) == 0) {
    return(invisible())
  }
  scores <- gaps$scores[[1]]
  low <- gaps$low[1]
  missing <- paste(setdiff(seq(low + 1, max(scores)), scores), collapse = ", ")
  what <- if (low == 0) "has scores 0, " else sprintf("is fixed up to step %d and has scores ", low)
  how <- if (low == 0) "" else ", or fix the missing step(s) too"
  kits_abort("Item ", gaps$item_id[1], " ", what, paste(scores, collapse = ", "),
             " (missing ", missing, "). Scores must be consecutive (0, 1, 2, ...). ",
             "Recode the item first", how, ".", class = "kits_error_data")
}

# 固定題的最高步驟，照 items 的順序（不是固定的題是 NA）
fixed_top <- function(fixed, items) {
  if (is.null(fixed)) {
    return(rep(NA_integer_, length(items)))
  }
  top <- summarise(fixed, top = as.integer(max(item_score)), .by = item_id)
  top$top[match(items, top$item_id)]
}

# info 的設定裡固定參數的說明
fixed_text <- function(fixed) {
  if (is.null(fixed)) {
    return("none")
  }
  sprintf("%d step(s) of %d item(s)", nrow(fixed), n_distinct(fixed$item_id))
}

# 作答組合（booklet）：每個人答的題目組合（包括沒放進校準的固定題），一個作答組合一列：
# set（題號的鍵）、n（人次）、first_person（第一個人）、pattern_id。
# 編號和 kits_subj() 相同（bk1, bk2, ...，booklet_names()）：人次照 visit_id 排序（C locale，arrange() 的預設），
# 作答組合照出現的順序編號。所以同樣的資料，kits_subj() 和 kits_subj_table(item) 的 pattern_id 相同。
calibration_booklets <- function(resp) {
  number_booklets(person_sets(resp))
}

# 題本的設計表：pattern_id、item_id，一列一題。
# 列名和 kitsPreview 相同：題本第一個人的 visit_id，題本有好幾題時再加上 1, 2, ...（unlist() 的命名）。
# 不帶任何人次編號（題目參數物件要能分享）：列名是 1, 2, ...
booklet_design <- function(books) {
  items <- map(books$set, set_items)
  data.frame(pattern_id = rep(books$pattern_id, lengths(items)), item_id = unlist(items))
}

# kits_item() 的診斷；book_n：每個作答組合（booklet）的人次；fit_row：題目適配那一列
item_diagnostics <- function(fit, calib, left_out, n_persons, book_n, fit_row = NULL) {
  n_books <- length(book_n)
  o <- fit$optim
  rows <- list()
  if (is.null(o) || is.na(o$convergence)) {
    rows$conv <- diag_row("optimizer convergence", NA, "ok", "every step is fixed, so nothing was estimated")
  } else {
    reached <- o$max_gradient <= 1e-4 && isTRUE(o$hessian_pd)
    status <- if (o$convergence == 0 || reached) "ok" else "warn"
    message <- if (o$convergence == 0) {
      sprintf("converged (nlminb: %s)", o$message)
    } else if (reached) {
      sprintf(paste0("converged: nlminb stopped early (code %d) at its very strict tolerance; ",
                     "the Newton steps reached the maximum (max |gradient| %.1e < 1e-4, Hessian ",
                     "positive definite)"),
              o$convergence, o$max_gradient)
    } else {
      sprintf("nlminb stopped with code %d (%s); the estimates may not be at the maximum", o$convergence, o$message)
    }
    rows$conv <- diag_row("optimizer convergence", o$convergence, status, message)
    steep <- o$max_gradient > 1e-4
    grad_status <- if (steep) "warn" else "ok"
    grad_message <- if (steep) "the gradient at the solution is above 1e-4: not at the maximum" else "below 1e-4"
    rows$grad <- diag_row("max |gradient|", o$max_gradient, grad_status, grad_message)
    pd <- isTRUE(o$hessian_pd)
    hess_status <- if (pd) "ok" else "warn"
    hess_message <- if (pd) "yes" else "no: some parameters are not identified; the SEs cannot be trusted"
    rows$hess <- diag_row("Hessian positive definite", as.numeric(o$hessian_pd), hess_status, hess_message)
    rows$extreme <- extreme_row(fit$coef)
  }
  # 每題每個分數的人數（0 到這題觀察到的最高分；沒有人的分數是 0 人）
  cells <- calib |>
    reframe(score = 0:max(item_score), n = tabulate(item_score + 1, nbins = max(item_score) + 1), .by = item_id) |>
    arrange(sort_rank(item_id), score) |>
    rename(item = item_id)
  low <- filter(cells, n < 5)
  rows$cat <- category_row(cells, low)
  any_left <- length(left_out) > 0
  left_status <- if (any_left) "note" else "ok"
  left_message <- if (any_left) {
    paste0("no 0 score or one score; kept at their fixed values: ", paste(left_out, collapse = ", "))
  } else {
    "none"
  }
  rows$left <- diag_row("fixed items left out", length(left_out), left_status, left_message)
  rows$persons <- diag_row("persons", n_persons, "ok", sprintf("%d person-waves in the calibration", n_persons))
  books_text <- sprintf("%d item pattern(s) (the combinations of items answered)", n_books)
  rows$books <- diag_row("item patterns", n_books, "ok", books_text)
  rows$book_n <- booklet_size_row(book_n)
  rows$itemfit <- fit_row
  diag_table(rows)
}

# 每題的 flag（itemfit_flag()，plot() 的範圍）：從 infit、outfit 算，所以沒有 flag 欄的舊物件也可以
itemfit_flags <- function(itemfit) itemfit_flag(itemfit$infit, itemfit$outfit)

# 題目適配在可接受範圍（plot() 的範圍，itemfit_range()）外的題：underfit 的題在前，再來 overfit；
# 算不出來的（NA）不算
itemfit_outside <- function(itemfit) {
  flag <- itemfit_flags(itemfit)
  c(itemfit$item_id[flag %in% "underfit"], itemfit$item_id[flag %in% "overfit"])
}

# 範圍外的題數分成兩種，例如 "2 underfit, 7 overfit"；0 題的那種不寫
itemfit_split <- function(itemfit) {
  flag <- itemfit_flags(itemfit)
  n <- c(underfit = sum(flag %in% "underfit"), overfit = sum(flag %in% "overfit"))
  paste(n[n > 0], names(n)[n > 0], collapse = ", ")
}

# 題目適配的種類，印在 Item fit 後面，例如 "Item fit (conditional)"；沒有記錄（舊的物件）時只有 "Item fit"
itemfit_label <- function(itemfit) {
  type <- attr(itemfit, "type")
  if (is.null(type)) "Item fit" else sprintf("Item fit (%s)", type)
}

# 有 infit 和 outfit 的題數
itemfit_n <- function(itemfit) sum(is.finite(itemfit$infit) & is.finite(itemfit$outfit))

# 範圍的文字，例如 0.6-1.4
itemfit_range_text <- function() paste(format(itemfit_range()), collapse = "-")

# 診斷：題目適配。範圍外的題是 note（KIT 的資料本來就有幾題不適配，所以不是 warn）；算不出來時也是 note
itemfit_row <- function(itemfit) {
  fit <- itemfit$table
  if (is.null(fit)) {
    return(diag_row("item fit", NA, "note", paste0("not computed (x$itemfit is NULL): ", itemfit$reason)))
  }
  bad <- itemfit_outside(fit)
  n <- itemfit_n(fit)
  rng <- itemfit_range_text()
  type <- attr(fit, "type")   # 例如 conditional；舊的物件沒有記錄
  if (length(bad) == 0) {
    what <- paste(c(type, "infit and outfit"), collapse = " ")
    return(diag_row("item fit", 0, "ok", sprintf("%s within %s for all %d items", what, rng, n)))
  }
  flag <- itemfit_flags(fit)
  kinds <- c("underfit", "overfit")
  groups <- map_chr(kinds, \(k) {   # 例如 2 underfit (a3, a5)；0 題的那種不寫
    ids <- fit$item_id[flag %in% k]
    if (length(ids) == 0) return(NA_character_)
    sprintf("%d %s (%s)", length(ids), k, some_ids(ids, 10))
  })
  msg <- sprintf("%d of %d items outside %s (%s): %s; see x$itemfit, %s",
                 length(bad), n, rng, paste(c(type, "infit or outfit"), collapse = " "),
                 paste(stats::na.omit(groups), collapse = ", "), "plot(x, type = \"fit\")")
  diag_row("item fit", length(bad), "note", msg)
}

# 印出的題目適配一行（太長時兩行）；where 是物件的名字（kits() 裡是 x$item）。
# 沒有 itemfit 這個元素（舊的物件）時不印
itemfit_lines <- function(x, where = "x") {
  if (!"itemfit" %in% names(x)) {
    return(character())
  }
  fit <- x$itemfit
  if (is.null(fit)) {
    return(sprintf("Item fit: not computed; see summary(%s).", where))
  }
  bad <- itemfit_outside(fit)
  rng <- itemfit_range_text()
  label <- itemfit_label(fit)
  if (length(bad) == 0) {
    return(sprintf("%s: infit and outfit within %s for all %d items.", label, rng, itemfit_n(fit)))
  }
  head <- sprintf("%s: %d of %d items outside %s (%s: %s);", label, length(bad), itemfit_n(fit), rng,
                  itemfit_split(fit), some_ids(bad, 5))
  see <- sprintf("see %s$itemfit, plot(%s, type = \"fit\").", where, where)
  lines <- strwrap(head, width = 101, exdent = 2)   # 每行最多 100 字；see 不拆開
  last <- paste(lines[length(lines)], see)
  if (nchar(last) <= 100) c(utils::head(lines, -1), last) else c(lines, paste0("  ", see))
}

# 每題每個分數的人數：少於 5 人的分數註記
category_row <- function(cells, low) {
  any_low <- nrow(low) > 0
  status <- if (any_low) "note" else "ok"
  cells_text <- paste0(low$item, ":", low$score, " (", low$n, ")")
  message <- if (any_low) {
    paste0("fewer than 5 persons in: ", some_ids(cells_text, 10))
  } else {
    "every score of every item has at least 5 persons"
  }
  diag_row("smallest category count", min(cells$n), status, message)
}

# 很大的估計或 SE（自由估的步驟）：估計可能不存在或幾乎不存在（例如一題幾乎只有滿分的孩子通過）
extreme_row <- function(cf, max_se = 10, max_beta = 15) {
  free <- cf$SE_beta > 0
  if (!any(free)) {
    return(NULL)
  }
  big <- free & (cf$SE_beta > max_se | abs(cf$beta) > max_beta)
  steps <- paste0(cf$item_id, ":", cf$item_score)[big]
  largest_se <- max(cf$SE_beta[free])
  status <- if (any(big)) "warn" else "ok"
  message <- if (any(big)) {
    sprintf(paste("largest SE %.3g; step(s) with SE above %g or |beta| above %g: %s. The estimate may not exist",
                  "(a step reached or missed almost only by persons with a zero or perfect score); check the data"),
            largest_se, max_se, max_beta, paste(utils::head(steps, 10), collapse = ", "))
  } else {
    sprintf("largest SE %.3g, largest |beta| %.3g", largest_se, max(abs(cf$beta[free])))
  }
  diag_row("extreme estimates", largest_se, status, message)
}

# 題本的人數（Kolen & Brennan p. 304：Rasch 模型每個題本約 400 人）。只看主要的題本（至少 5% 的人次）；
# 其他的題本多半是少答幾題的組合，只在訊息裡說有幾個。
booklet_size_row <- function(book_n, min_n = 400, main = 0.05) {
  is_main <- book_n >= main * sum(book_n)
  if (!any(is_main)) is_main <- book_n == max(book_n)   # 沒有一個題本有 5%（很多題本的設計）：看最大的
  small <- is_main & book_n < min_n
  other <- sum(!is_main)
  others <- if (other > 0) sprintf("; %d other item pattern(s) hold %.0f%%", other,
                                   100 * sum(book_n[!is_main]) / sum(book_n)) else ""
  msg <- sprintf(paste0("%d main item pattern(s) (5%% or more of the person-waves, or the ",
                        "largest): %s person-waves%s"), sum(is_main),
                 paste(sort(book_n[is_main], decreasing = TRUE), collapse = ", "), others)
  if (any(small)) {
    msg <- paste0(msg, sprintf("; %d below %d, the rule of thumb for the Rasch model (Kolen & Brennan p. 304)", sum(small),
                               min_n))
  }
  status <- if (any(small)) "note" else "ok"
  diag_row("smallest main item pattern", min(book_n[is_main]), status, msg)
}


#' @rdname kits_item
#' @param x,object A `kits_item`.
#' @param digits Decimals to print (only the printout is rounded; the object keeps full precision).
#' @param n Steps to print (the first 6 by default; `Inf` for all). Fixed steps are marked in a
#'   `fixed` column.
#' @export
print.kits_item <- function(x, digits = 2, n = 6, ...) {
  need_digits(digits)
  need_rows(n)
  if (!print_item_head(x, digits, n)) {
    return(invisible(x))
  }
  cat(itemfit_lines(x), sep = "\n")
  cat("Next: kits_subj(data, x), kits_linking(x, ...), plot(x, type = \"fit\"), plot(x), ",
      "confint(x).\n", sep = "")
  checks_line(x$diagnostics)
  invisible(x)
}

# print.kits_item() 和 print.kits() 的前半：標題、量尺、門檻表、沒放進校準的題。物件不完整時 FALSE
print_item_head <- function(x, digits, n) {
  thr <- x$thresholds
  if (!is.data.frame(thr) || is.null(x$fit)) {
    cat("<kits_item> (incomplete object)\n")
    return(FALSE)
  }
  cat(sprintf(paste0("<kits_item> PCM (kits, CML): %d person-waves, %d items, %d steps, %d item ",
                     "pattern(s)\n"),
              x$n_persons, length(unique(thr$item_id)), nrow(thr),
              length(unique(x$fit$inputs$design$pattern_id))))
  beta_scale <- beta_scale_of(x)
  if (!identical(beta_scale, scale_of(x))) {
    cat(sprintf("Scale: %s (thr); beta stays on %s\n", scale_label(scale_of(x)), scale_label(beta_scale)))
  } else {
    cat(sprintf("Scale: %s\n", scale_label(scale_of(x))))
  }
  shown <- format_table(utils::head(thr, n), digits)
  fixed <- fixed_steps(x)
  if (length(fixed) > 0) {
    shown$fixed <- ifelse(paste(utils::head(thr$item_id, n), utils::head(thr$item_score, n)) %in% fixed, "yes", "")
    cat(sprintf("Fixed (fixed_params, SE 0): %d step(s) of %d item(s)\n", length(fixed),
                length(unique(sub(" .*$", "", fixed)))))
  }
  print(shown, row.names = FALSE, right = TRUE)
  if (nrow(thr) > n) cat(sprintf("... %d more steps (more: print(x, n = %d); all steps: coef(x))\n", nrow(thr) - n, min(nrow(thr),
      50)))
  if (length(x$left_out) > 0) cat("Fixed, left out of the calibration (no 0 score or one score):", paste(x$left_out,
      collapse = ", "), "\n")
  TRUE
}

# 固定的步驟（"item_id item_score"）：info 的設定記下固定參數；門檻表裡 SE 是 0 的自由步驟不會有
fixed_steps <- function(x) {
  thr <- x$thresholds
  if (identical(x$info$settings$fixed_params, "none") || is.null(x$info$settings$fixed_params)) {
    return(character())
  }
  paste(thr$item_id, thr$item_score)[thr$beta_se == 0]
}

#' @rdname kits_item
#' @param rp The response probability of the thresholds: `thr` of step k is the theta at which
#'   P(X >= k) = `rp`. Default 0.5 (RP50, as calibrated); links use RP50 only.
#' @export
coef.kits_item <- function(object, rp = 0.5, ...) {
  need_rp(rp)
  if (rp == 0.5) {
    return(object$thresholds)
  }
  abort_if(!identical(beta_scale_of(object), scale_of(object)),
           "This kits_item was converted with kits_convert(): compute coef(item, rp = ) on the calibration, ",
           "then convert the table.", class = "kits_error_scale")
  item_thresholds_rp(object, rp)
}


#' Calibrate and score in one step
#'
#' [kits_item()] and then [kits_subj()] on the same data.
#' Use it when no link is needed in between (for example a single cohort).
#' With linking, call the two functions separately.
#' `plot()` draws a Wright map: the abilities and the step thresholds on the same theta axis.
#'
#' Details: `kits_algorithm("cml")`, `kits_algorithm("wle")`, tutorial_algorithms/01_calibration.html
#' and 02_scoring.html.
#'
#' @inheritParams kits_item
#' @inheritParams kits_subj
#' @return A `kits`: a list with `item` (the [kits_item()] result), `subj` (the [kits_subj()] result),
#'   `scale`, `diagnostics` (those of `item` and `subj`) and `info`.
#' @examples
#' betas <- list(a1 = c(-1, 0), a2 = c(0, 1), a3 = 0.5, a4 = -0.5)
#' X <- kits_sim(seq(-3, 3, length.out = 300), betas, seed = 1)
#' k <- kits(X)
#' k
#' head(k$subj)
#' p <- plot(k)
#' if (interactive()) p
#' @export
kits <- function(data, ..., scale = NULL, persons = NULL) {
  check_dots(match.call(expand.dots = FALSE)$..., "kits")
  call <- call_text(sys.call())
  data_expr <- substitute(data)
  resp <- resp_mark(resp_long(data, gaps = is.null(list(...)$fixed_params)))   # 已經檢查過、轉成長格式
  input <- data_input(data, data_expr)
  item <- calibrate(resp, list(...), scale, default_scale(data_expr, data), call, input)
  if (is.null(persons)) persons <- resp_persons(data)   # 作答表：孩子的欄跟著走
  subj <- kits_subj(resp, item, persons = persons)
  info <- attr(subj, "info")
  info$call <- call
  info$input <- input
  attr(subj, "info") <- info
  pre <- function(d, what) {
    d$check <- paste0(what, ": ", d$check)
    d
  }
  diagnostics <- kits_table(rbind(pre(item$diagnostics, "item"), pre(attr(subj, "diagnostics"), "subj")), "kits_diagnostics")
  structure(list(item = item, subj = subj, scale = item$scale, diagnostics = diagnostics, info = item$info), class = "kits")
}


#' @rdname kits
#' @param x A `kits`.
#' @param digits Decimals to print (only the printout is rounded; the object keeps full precision).
#' @export
print.kits <- function(x, digits = 2, ...) {
  need_digits(digits)
  # 題目參數的前半和題目適配；最後只印一行 Checks（item 和 subj 合起來），下一步是計分之後的
  if (print_item_head(x$item, digits, n = 6)) cat(itemfit_lines(x$item, "x$item"), sep = "\n")
  th <- x$subj$theta
  if (length(th) == 0) {
    return(invisible(x))
  }
  cat(sprintf("Abilities (WLE, $subj): %d person-waves, theta mean %s, SD %s, range %s to %s\n",
              length(th), fmt(mean(th), digits), fmt(stats::sd(th), digits), fmt(min(th), digits),
              fmt(max(th), digits)))
  cat("Next: plot(x) (Wright map), x$subj |> kits_subj_table(), kits_algorithm(x),\n",
      "      kits_linking(x$item, ...) to link this scale to another.\n", sep = "")
  checks_line(x$diagnostics)
  invisible(x)
}


#' @rdname kits
#' @param ... Not used.
#' @export
plot.kits <- function(x, ...) {
  wright_plot(x$item$thresholds, x$subj$theta)
}

# Wright map（內部用；plot(kits) 和 plot(item, type = "wright", scores = )）：左邊孩子的能力（分組的人數），
# 右邊每題每個步驟的門檻，題目照 milestone step 由易到難排
wright_plot <- function(thresholds, theta) {
  levels <- thresholds |>
    filter(item_score == 1) |>
    arrange(thr) |>
    pull(item_id)
  thr <- mutate(as.data.frame(thresholds), x = match(item_id, levels))
  symbols <- c("circle", "square", "diamond", "triangle-up", "cross", "x")

  # 先分組算次數再畫，圖裡只有每組的人數，沒有個人的值
  bins <- graphics::hist(theta, breaks = 30, plot = FALSE)
  persons <- plotly::plot_ly(x = bins$counts, y = bins$mids, type = "bar", orientation = "h", marker = list(color = PAL$gray),
                             name = LAB$persons, hovertext = sprintf("%.2f to %.2f: %d", utils::head(bins$breaks, -1),
                                                                         bins$breaks[-1], bins$counts), hoverinfo = "text")
  items <- plotly::plot_ly()
  for (k in sort(unique(thr$item_score))) {
    d <- filter(thr, item_score == k)
    items <- plotly::add_markers(items, x = d$x, y = d$thr, name = paste(LAB$step, k),
                                 marker = list(color = PAL$blue, size = MARKER_SIZE,
                                               symbol = symbols[(k - 1) %% length(symbols) + 1]),
                                 hovertext = sprintf("%s, %s %d: %.2f", d$item_id, LAB$step, k, d$thr), hoverinfo = "text")
  }
  # 兩格之間只留一點空隙：右邊那格和左邊共用 theta 軸，沒有自己的刻度字
  plotly::subplot(persons, items, shareY = TRUE, widths = c(0.25, 0.75), titleX = TRUE, margin = c(0.004, 0.004, 0, 0)) |>
    plotly::layout(xaxis = list(title = LAB$persons, autorange = "reversed"),
                   xaxis2 = list(title = "", tickvals = seq_along(levels), ticktext = levels, tickangle = -90),
                   yaxis = list(title = LAB$theta)) |>
    kits_style(LAB$title_wright)
}
