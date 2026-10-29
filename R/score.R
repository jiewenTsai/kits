#' Raw-score to ability table
#'
#' `scores |> kits_subj_table()` gives the table of each wave from the scores of [kits_subj()]:
#' the WLE and the count of every raw score of the standard form (the item pattern answered by most
#' children of the wave), on the scale of the scores. `kits_subj_table(item)` gives the WLE of
#' every raw score of every item pattern of a calibration, without data.
#'
#' Details: tutorial_algorithms/02_scoring.html.
#'
#' @param x The scores of [kits_subj()] (converted or not, filtered or stacked), or the item
#'   parameters:
#'   a [kits_item()] result, a [kits()], its `fit` ([kits_cml()]), or a data.frame with `item_id`,
#'   `item_score`, `beta`. In `print()`: a `kits_subj_table`.
#' @param ... Not used: an argument the function does not have stops with an error (class
#'   `kits_error_input`); the responses are not an argument (score them first with [kits_subj()]).
#' @return From scores: a `kits_subj_table` (a data.frame); one row is one raw score of the standard
#'   form in one wave: `wave_m`, `n_items`, `raw`, `theta`, `se` (and `se_link`, `se_total` when the
#'   conversions add them), `n` (children with this raw score on the standard form), `pct`,
#'   `cum_pct`, `wave_n` (all children in the wave), `form_pct` (percentage of the wave's children
#'   on the standard form).
#'   `attr(, "forms")` lists every item pattern of each wave: `wave_m`, `form` (1 = the standard
#'   form), `type` (`"missing responses"` or `"other items"`), `n`, `items`.
#'   It prints with a note on the children who are not in the table, and keeps a diagnostics table.
#'
#'   From item parameters: a table (a data.frame) with one row per item pattern and per raw score (0
#'   to the maximum):
#'   `pattern_id`, `raw`, `theta`, `se`.
#' @examples
#' betas <- list(a1 = c(-1, 0), a2 = c(0, 1), a3 = 0.5, a4 = -0.5)
#' X <- kits_sim(seq(-3, 3, length.out = 300), betas, seed = 1)
#' X[1:10, "a4"] <- NA
#' item <- kits_item(X)
#'
#' # The standard form of each wave, from the scores
#' persons <- data.frame(visit_id = rownames(X), wave_m = rep(c(12, 24), each = 150))
#' scores <- kits_subj(X, item, persons = persons)
#' scores |>
#'   kits_subj_table()
#'
#' # Every raw score of every item pattern, without data
#' kits_subj_table(item)
#'
#' # A table of item parameters: all its items are one item pattern
#' parms <- data.frame(item_id = c("a1", "a1", "a2"), item_score = c(1, 2, 1), beta = c(-1, 0, 0.5))
#' kits_subj_table(parms)
#' @export
kits_subj_table <- function(x, ...) {
  no_extra_args(..., fun = "kits_subj_table")
  call <- call_text(sys.call())
  if (inherits(x, "kits_subj")) {
    return(subj_table_scores(x, call))
  }
  wle_table(x, call = call)
}

# 分數表的原始分換算表：每一波的標準題本（人數最多的題目組合，在自己的校準裡比）、它的每個原始分的 WLE
# （和原本 kits_subj_table(item, data, wave_m) 同一個計算），再照分數表記下的換算換尺
subj_table_scores <- function(scores, call) {
  sources <- attr(scores, "raw_theta")
  abort_if(is.null(sources), "These scores keep no raw-score tables ",
           "(made by an older kits, or not by kits_subj()); score them again with kits_subj().")
  need_cols(scores, c("visit_id", "pattern_id", "raw"), "scores")
  rows <- score_sets(as.data.frame(scores), sources)
  waves <- wave_values(rows$wave_m)
  per_wave <- map(waves, \(w) scores_wave_table(filter(rows, wave_m %in% w), w, sources))
  out <- per_wave |>
    map("table") |>
    list_rbind() |>
    kits_table("kits_subj_table")
  sets <- per_wave |>
    map("forms") |>
    list_rbind()
  attr(out, "forms") <- sets
  missing <- sum(filter(sets, type == "missing responses")$n)
  note <- sprintf(paste0("%d person-wave(s) with missing responses are not in the table; ",
                         "their scores are in the scores"), missing)
  diag <- c(list(diag_row("missing responses", missing, "ok", note)),
            compact(map(per_wave, "other_items")))
  attr(out, "scale") <- scale_of(scores)
  attr(out, "diagnostics") <- diag_table(diag)
  ids <- values_fingerprint(scores$visit_id)
  attr(out, "info") <- make_info(call, list(scale = scale_of(scores)), ids)
  out
}

# 每一列：哪一個校準（source）、答的題目組合（set）、原始分、波次
score_sets <- function(d, sources) {
  owner <- sources |>
    imap(\(src, i) data.frame(visit_id = src$persons, source = i)) |>
    list_rbind()
  twice <- unique(owner$visit_id[duplicated(owner$visit_id)])
  abort_if(any(d$visit_id %in% twice), "Scores of two calibrations share person ids, ",
           "so their item patterns cannot be told apart; ",
           "make the raw-score table of each before stacking.", class = "kits_error_data")
  d$source <- owner$source[match(d$visit_id, owner$visit_id)]
  abort_if(anyNA(d$source),
           paste0("Some rows of the scores do not come from kits_subj(); their item patterns ",
                  "are unknown."), class = "kits_error_data")
  books <- sources |>
    imap(\(src, i) mutate(src$booklets, source = i)) |>
    list_rbind()
  d |>
    left_join(books, by = c("source", "pattern_id")) |>
    transmute(visit_id, source, set, key = paste0(source, "|", set), raw,
              wave_m = if ("wave_m" %in% names(d)) wave_m else NA_real_)
}

# 一個波次：標準題本的原始分換算表、每個題目組合、標準題本以外的題的警告
scores_wave_table <- function(rows, w, sources) {
  forms <- rows |>
    mutate(n_wave = n()) |>
    summarise(n = n(), raw = list(raw), n_wave = n_wave[1], source = source[1], set = set[1],
              .by = key) |>
    arrange(desc(n), sort_rank(key)) |>
    mutate(items = map(set, set_items))
  standard <- forms$items[[1]]
  src <- sources[[forms$source[1]]]
  form <- data.frame(pattern_id = "form", item_id = standard)
  wle <- wle_table(filter(src$parms, item_id %in% standard), design = form)
  table <- data.frame(wave_m = w, n_items = length(standard), raw = wle$raw,
                      theta = wle$theta, se = wle$se,
                      n = tabulate(forms$raw[[1]] + 1, nbins = nrow(wle))) |>
    mutate(pct = 100 * n / sum(n),
           cum_pct = cumsum(pct),
           wave_n = forms$n_wave[1],
           form_pct = 100 * sum(n) / wave_n)
  forms <- forms |>
    mutate(form = row_number(),
           inside = map_lgl(items, \(its) all(its %in% standard)),
           type = case_when(form == 1 ~ "standard", inside ~ "missing responses",
                            .default = "other items"))
  list(table = replay_conversions(table, src),
       forms = transmute(forms, wave_m = w, form, type, n,
                         items = map_chr(items, \(its) paste(its, collapse = ", "))),
       other_items = other_items_row(filter(forms, type == "other items"), standard, w))
}

# 分數表換過的尺，換算表也照同樣的連結、同樣的順序換（訊息和警告在換分數時已經說過）
replay_conversions <- function(table, src) {
  table <- kits_table(table)
  attr(table, "scale") <- src$scale
  for (step in src$steps) {
    table <- suppressMessages(suppressWarnings(convert_links(table, step$links, step$from),
                                               classes = "kits_warning_scale"),
                              classes = "kits_message_scale")
  }
  as.data.frame(table)
}

# 幾張表的原始分 -> theta 的來源合在一起（rbind()、bind_rows()）
combine_raw_theta <- function(parts) {
  sources <- unlist(map(parts, \(p) if (is.data.frame(p)) attr(p, "raw_theta")), recursive = FALSE)
  if (length(sources) == 0) NULL else unique(sources)
}

# 有標準題本以外的題的人次：資料可能有問題（題本或波次錯了），所以是警告
other_items_row <- function(other, standard, w) {
  if (nrow(other) == 0) {
    return(NULL)
  }
  extra <- sort(unique(setdiff(unlist(other$items), standard)))
  diag_row(paste("other items, wave", wave_label(w)), sum(other$n), "warn",
           sprintf(paste0("wave %s: %d person-wave(s) answered items outside the standard form ",
                          "(%s); check their wave and item pattern"),
                   wave_label(w), sum(other$n), paste(extra, collapse = ", ")))
}

# 波次，由小到大；一組（NA）時只有 NA
wave_values <- function(w) sort(unique(w), na.last = TRUE)

# 波次的標籤：一組（NA）時是 "all"
wave_label <- function(w) map_chr(unname(w), \(v) if (is.na(v)) "all" else format(v))


#' @rdname kits_subj_table
#' @export
print.kits_pattern_table <- function(x, digits = 2, n = NULL, ...) {
  need_digits(digits)
  need_rows(n)
  cat("<kits_subj_table> every raw score of every item pattern, from the item parameters; ",
      sprintf("scale: %s\n", scale_label(scale_of(x))), sep = "")
  print_rows(x, digits, n, ...)
  invisible(x)
}

#' @rdname kits_subj_table
#' @param digits Decimals to print (only the printout is rounded; the object keeps full precision).
#' @param n Rows to print: by default all rows of a table of up to 30 rows, else the first 10.
#' @param ... In `kits_subj_table()`: not used (an unknown argument stops with an error).
#'   In `print()`: passed to `print()` for a data.frame.
#' @export
print.kits_subj_table <- function(x, digits = 2, n = NULL, ...) {
  need_digits(digits)
  need_rows(n)
  cat("<kits_subj_table> the standard form of each wave of the data, with counts; ",
      sprintf("scale: %s\n", scale_label(scale_of(x))), sep = "")
  print_rows(x, digits, n, ...)
  sets <- attr(x, "forms")
  if (!is.null(sets) && any(sets$type != "standard")) {
    cat("Note: each table converts complete responses to the standard form only.\n")
    for (w in unique(sets$wave_m)) {
      s <- filter(sets, wave_m %in% w)
      n_wave <- sum(s$n)
      missing <- sum(s$n[s$type == "missing responses"])
      other <- filter(s, type == "other items")
      if (missing > 0) {
        cat(sprintf(paste0("  Wave %s: %d of %d person-waves with missing responses are not in the table; score them with ",
                           "kits_subj().\n"),
                    wave_label(w), missing, n_wave))
      }
      if (nrow(other) > 0) {
        standard <- strsplit(s$items[s$type == "standard"], ", ", fixed = TRUE)[[1]]
        extra <- sort(unique(setdiff(unlist(strsplit(other$items, ", ", fixed = TRUE)), standard)))
        cat(sprintf(paste0("  Wave %s: %d of %d person-waves answered items outside the standard form (%s); check their wave ",
                           "and item pattern.\n"),
                    wave_label(w), sum(other$n), n_wave, paste(extra, collapse = ", ")))
      }
    }
  }
  checks_line(attr(x, "diagnostics"))
  invisible(x)
}


# 一個波次的作答組合，照人數由多到少（同樣多時照題號排）；第一個是標準題本（內部用）。
# 每個題本一列：items（題號）、n（人次）、raw（答這個題本的人的原始分）、n_wave（這個波次的總人數）。
wave_forms <- function(resp) {
  resp |>
    person_sets() |>
    mutate(n_wave = n()) |>
    summarise(n = n(), raw = list(raw), n_wave = n_wave[1], .by = set) |>
    arrange(desc(n), sort_rank(set)) |>
    mutate(items = map(set, set_items))
}
