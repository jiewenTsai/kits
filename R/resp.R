# 作答表（kits_resp）：一列一個人次，前面是孩子的欄，後面是題目的欄；kits 記得哪些欄是題目（attr(, "items")）。
# 孩子的資料跟著作答走：kits_subj()、kits_item_table()、kits_dif() 直接拿它的欄，不用自己對齊。

# 孩子的欄（有的話放在前面，照這個順序）
RESP_PERSON_COLS <- c("visit_id", "subj_id", "cohort", "sex", "wave_m", "age_m")

#' Responses with the information of each child
#'
#' `kits_resp()` puts the responses of one domain and the information of each child in one table,
#' one row per person-wave: first the columns of the child (`visit_id`, `subj_id`, `cohort`, `sex`,
#' `wave_m`, `age_m` and any other column that does not change within a person-wave), then one
#' column per item. kits remembers which columns are items (`attr(, "items")`).
#' The later functions ([kits_item()], [kits_subj()], [kits_item_table()], [kits_dif()]) take the
#' information of the child from this table; it is optional, and `filter()` and `mutate()` keep it a
#' `kits_resp`.
#' `print()` gives the counts and the checks, `summary()` the standard form of each wave, and
#' `plot()` the design (which items each wave answered).
#' Details: `vignette("kits-flow", package = "kits")`, section 0.
#'
#' @param data Responses: a long data.frame (`visit_id`, `item_id`, `item_score` and the columns of
#'   the child; one domain at a time), or a matrix or a wide data.frame (see [kits_item()]).
#' @param persons Optional information of each person-wave, a data.frame with `visit_id` (one row
#'   per person-wave; extra rows are fine). Its columns are added to the table.
#' @return A `kits_resp`: a data.frame with one row per person-wave, the columns of the child first,
#'   then the items; `attr(, "items")` (the item columns), `attr(, "domain")` (when the data had
#'   one), diagnostics and [kits_info()].
#' @seealso [kits_item()], [kits_subj()], [kits_item_table()], [kits_dif()].
#' @examples
#' library(dplyr)
#' resp <- kitm03$responses |>
#'   filter(domain == "lan", wave_m <= 12) |>
#'   kits_resp()
#' resp
#' summary(resp)
#' p <- plot(resp)
#' if (interactive()) p
#'
#' # The information of the child goes along
#' item <- kits_item(resp, scale = "lan M03")
#' scores <- kits_subj(resp, item)
#' head(scores)
#'
#' # A matrix with the information given separately
#' X <- kits_sim(rnorm(200), list(a1 = 0, a2 = c(-1, 1), a3 = 0.5), seed = 1)
#' persons <- data.frame(visit_id = rownames(X), wave_m = rep(c(6, 12), each = 100))
#' kits_resp(X, persons = persons)
#' @export
kits_resp <- function(data, persons = NULL) {
  call <- call_text(sys.call())
  if (inherits(data, "kits_resp")) {
    return(data)
  }
  data <- old_ids(data)
  persons <- old_ids(persons, "persons")
  out <- if (is_long(data)) resp_table_long(data) else resp_table_wide(data)
  out <- add_resp_persons(out$table, out$items, persons, out$domain)
  long <- resp_long(out, gaps = FALSE)   # 和其他函數同樣的檢查
  attr(out, "diagnostics") <- resp_diagnostics(out, long)
  attr(out, "info") <- make_info(call, list(domain = attr(out, "domain")), data_fingerprint(long))
  out
}

# 長格式：孩子的欄（同一個人次不變的欄）在前，題目一題一欄
resp_table_long <- function(data) {
  need_cols(data, c("visit_id", "item_id", "item_score"), "data")
  d <- as.data.frame(data)
  domain <- one_domain(d)
  d$domain <- NULL
  d$visit_id <- as.character(d$visit_id)
  d$item_id <- as.character(d$item_id)
  abort_if(vctrs::vec_duplicate_any(d[c("visit_id", "item_id")]),
           "`data` has a person with the same item more than once. ",
           "Each person-wave needs its own visit_id, e.g. paste(child_id, wave_m).", class = "kits_error_data")
  people <- distinct(d[setdiff(names(d), c("item_id", "item_score"))])
  abort_if(anyDuplicated(people$visit_id), "Column(s) ", paste(varying_columns(people), collapse = ", "),
           " of `data` change within a person-wave; ",
           "the columns of the child must be the same in every row of the person-wave.", class = "kits_error_data")
  items <- sort(unique(d$item_id))
  clash <- intersect(items, names(people))
  abort_if(length(clash) > 0, "Item id(s) with the name of a column of the child: ", paste(clash, collapse = ", "), ".",
           class = "kits_error_data")
  wide <- d |>
    select(visit_id, item_id, item_score) |>
    pivot_wider(names_from = item_id, values_from = item_score)
  table <- left_join(people, wide, by = "visit_id")
  list(table = table[c(names(people), items)], items = items, domain = domain)
}

# 一次一個領域：domain 欄有好幾個值時停下來，說怎麼篩
one_domain <- function(d) {
  if (!"domain" %in% names(d)) {
    return(NULL)
  }
  domains <- unique(stats::na.omit(as.character(d$domain)))
  abort_if(length(domains) > 1, "`data` holds several domains (", paste(domains, collapse = ", "),
           "); a kits_resp is one ", "domain. Keep one first, e.g. filter(data, domain == \"", domains[1], "\").")
  if (length(domains) == 0) NULL else domains
}

# 同一個人次有不同值的欄
varying_columns <- function(people) {
  cols <- setdiff(names(people), "visit_id")
  changes <- map_lgl(cols, \(v) anyDuplicated(distinct(people[c("visit_id", v)])$visit_id) > 0)
  cols[changes]
}

# 矩陣或寬的 data.frame：題目的欄；寬的 data.frame 裡孩子的欄拿出來，當成 persons
resp_table_wide <- function(data) {
  own <- NULL
  if (is.data.frame(data)) {
    cols <- intersect(names(data), setdiff(RESP_PERSON_COLS, "visit_id"))
    if (length(cols) > 0) {
      ids <- row_ids(data)
      own <- data.frame(visit_id = ids, as.data.frame(data)[cols])
      data <- as.data.frame(data)[setdiff(names(data), cols)]
    }
  }
  X <- wide_matrix(data, "data")
  if (is.null(rownames(X)) && !is.null(colnames(X))) {   # 沒有編號：照列的順序編 1 到 n
    numbered_message(nrow(X), "data")
    rownames(X) <- as.character(seq_len(nrow(X)))
  }
  table <- data.frame(visit_id = rownames(X), X, check.names = FALSE)
  if (!is.null(own)) table <- left_join(table, own, by = "visit_id")
  list(table = table, items = colnames(X), domain = NULL)
}

# persons 的欄接上；孩子的欄放在前面，題目在後面
add_resp_persons <- function(table, items, persons, domain) {
  if (!is.null(persons)) {
    abort_if(!is.data.frame(persons), "`persons` must be a data.frame with visit_id.")
    need_cols(persons, "visit_id", "persons")
    persons <- as.data.frame(persons)
    persons$visit_id <- as.character(persons$visit_id)
    abort_if(anyDuplicated(persons$visit_id),
             "`persons$visit_id` has repeated ids; give one row per person-wave.", class = "kits_error_data")
    miss <- setdiff(table$visit_id, persons$visit_id)
    abort_if(length(miss) > 0, length(miss), " person-wave(s) of the data have no row in `persons`",
             ids_text(miss, 3), ".", class = "kits_error_data")
    extra <- setdiff(names(persons), c(names(table), items))
    table <- left_join(table, persons[c("visit_id", extra)], by = "visit_id")
  }
  info <- setdiff(names(table), items)
  first <- intersect(RESP_PERSON_COLS, info)
  table <- table[c(first, setdiff(info, first), items)]
  rownames(table) <- NULL
  out <- kits_table(table, "kits_resp")
  attr(out, "items") <- items
  attr(out, "domain") <- domain
  cohorts <- if ("cohort" %in% names(out)) unique(stats::na.omit(as.character(out$cohort)))
  attr(out, "cohort") <- if (length(cohorts) == 1) cohorts
  out
}

# 作答表只有一個領域和一個世代時，預設的量尺標籤是 "領域 世代"（例如 "lan M03"）；否則 NULL
resp_scale <- function(data) {
  known <- inherits(data, "kits_resp") && !is.null(attr(data, "domain")) &&
    !is.null(attr(data, "cohort"))
  if (!known) {
    return(NULL)
  }
  paste(attr(data, "domain"), attr(data, "cohort"))
}

# 作答表的題目欄（照 attr；dplyr 的 select() 拿掉的題不算）
resp_items <- function(x) intersect(attr(x, "items"), names(x))

# 作答表裡孩子的欄
resp_person_cols <- function(x) setdiff(names(x), resp_items(x))

# 作答表的長格式，加上波次（沒有 wave_m 時所有列一組）
resp_by_wave <- function(x) {
  long <- resp_long(x, gaps = FALSE)
  waves <- if ("wave_m" %in% names(x)) x$wave_m else rep(NA_real_, nrow(x))
  long$wave_m <- waves[match(long$visit_id, x$visit_id)]
  long
}

# 每一波的作答組合：第一個是標準題本，其他的是少答（標準題本的一部分）或有其他題
resp_forms <- function(long) {
  waves <- wave_values(long$wave_m)
  map(waves, \(w) typed_forms(wave_forms(filter(long, wave_m %in% w)), w)) |>
    list_rbind()
}

# 作答組合的種類：standard、missing responses、other items（和 kits_subj_table() 相同）
typed_forms <- function(forms, w) {
  standard <- forms$items[[1]]
  forms |>
    mutate(form = row_number(),
           inside = map_lgl(items, \(its) all(its %in% standard)),
           type = case_when(form == 1 ~ "standard", inside ~ "missing responses",
                            .default = "other items")) |>
    transmute(wave_m = w, form, type, n, n_items = lengths(items),
              items = map_chr(items, \(its) paste(its, collapse = ", ")))
}

# 作答表的檢查：少答的人次（訊息）、標準題本以外的題（警告：題本或波次可能錯了）
resp_diagnostics <- function(x, long) {
  forms <- resp_forms(resp_by_wave(x))
  missing <- sum(forms$n[forms$type == "missing responses"])
  share <- 100 * missing / nrow(x)
  rows <- list(diag_row("missing responses", share, "ok",
                        sprintf(paste0("%.1f%% of person-waves answered part of the standard form ",
                                       "of their wave"), share)))
  other <- filter(forms, type == "other items")
  if (nrow(other) > 0) {
    waves <- paste(wave_label(unique(other$wave_m)), collapse = ", ")
    rows <- c(rows, list(diag_row("other items", sum(other$n), "warn",
      sprintf("%d person-wave(s) answered items outside the standard form of their wave (%s); %s",
              sum(other$n), paste("waves", waves), "check their wave and item pattern"))))
  }
  diag_table(rows)
}

#' @rdname kits_resp
#' @param x,object A `kits_resp`.
#' @param ... Not used.
#' @export
print.kits_resp <- function(x, ...) {
  items <- resp_items(x)
  answered <- !is.na(as.matrix(as.data.frame(x)[items]))
  sets <- nrow(unique(answered))
  domain <- attr(x, "domain")
  cat(sprintf("<kits_resp> %d person-waves, %d items, %d item pattern(s)%s\n", nrow(x),
              length(items),
              sets,
              if (is.null(domain)) "" else sprintf("; domain %s", domain)))
  if ("wave_m" %in% names(x)) {
    per_wave <- count(as.data.frame(x), wave_m)
    cat(sprintf("Person-waves per wave: %s\n",
                paste0(per_wave$wave_m, ": ", per_wave$n, collapse = ", ")))
  }
  cat(sprintf("Columns of the child: %s\n", paste(resp_person_cols(x), collapse = ", ")))
  cat("Next: kits_item(x, scale = ), kits_subj(x, item), kits_item_table(thr, norm, data = x), ",
      "kits_dif(x, by = \"sex\"), summary(x), plot(x).\n", sep = "")
  checks_line(attr(x, "diagnostics"))
  invisible(x)
}

#' @rdname kits_resp
#' @export
summary.kits_resp <- function(object, ...) {
  forms <- resp_forms(resp_by_wave(object))
  forms <- mutate(forms, pct = 100 * n / sum(n), .by = wave_m, .after = n)
  structure(list(forms = kits_table(forms), diagnostics = attr(object, "diagnostics")),
            class = "kits_resp_summary")
}

# 印出來的題號：標準題本印前幾題；少答的寫少了哪幾題，有其他題的寫多了哪幾題（物件裡留完整的題號）
short_items <- function(forms, k = 4) {
  ids <- strsplit(forms$items, ", ", fixed = TRUE)
  standard <- ids[match(paste(forms$wave_m, 1), paste(forms$wave_m, forms$form))]
  first_k <- \(its) some_ids(its, k)
  map_chr(seq_along(ids), \(i) {
    if (forms$form[i] == 1) return(first_k(ids[[i]]))
    less <- setdiff(standard[[i]], ids[[i]])
    more <- setdiff(ids[[i]], standard[[i]])
    paste(c(if (length(less) > 0) paste("without", first_k(less)),
            if (length(more) > 0) paste("plus", first_k(more))), collapse = "; ")
  })
}

#' @export
print.kits_resp_summary <- function(x, ...) {
  cat(paste0("Item patterns of each wave (form 1: the standard form, answered by most children of ",
             "the wave):\n"))
  shown <- as.data.frame(x$forms)
  shown$items <- short_items(shown)
  shown <- format_table(shown, 1)
  for (col in c("type", "items")) {   # 文字的欄靠左，欄名也是
    shown[[col]] <- format(shown[[col]])
    names(shown)[names(shown) == col] <- format(col, width = max(nchar(shown[[col]])))
  }
  print(shown, row.names = FALSE, right = TRUE)
  cat("All item ids: summary(x)$forms$items.\n")
  if (!is.null(x$diagnostics)) {
    cat("Diagnostics:\n")
    print(x$diagnostics)
  }
  invisible(x)
}

#' @rdname kits_resp
#' @export
plot.kits_resp <- function(x, ...) {
  long <- resp_by_wave(x) |>
    mutate(wave = wave_label(wave_m))
  cells <- long |>
    summarise(n = n(), pass = mean(item_score >= 1), .by = c(wave, wave_m, item_id)) |>
    arrange(wave_m)
  waves <- unique(cells$wave)
  items <- sort(unique(cells$item_id))
  z <- matrix(NA_real_, length(waves), length(items), dimnames = list(waves, items))
  z[cbind(match(cells$wave, waves), match(cells$item_id, items))] <- cells$pass
  hover <- matrix("", length(waves), length(items))
  hover[cbind(match(cells$wave, waves), match(cells$item_id, items))] <-
    sprintf("wave %s, %s: %.0f%% pass (n = %d)", cells$wave, cells$item_id, 100 * cells$pass,
            cells$n)
  scale <- list(c(0, "#F7FBFF"), c(1, PAL$blue))
  plotly::plot_ly(x = items, y = waves, z = z, type = "heatmap", colorscale = scale, zmin = 0,
                  zmax = 1,
                  text = hover, hoverinfo = "text", colorbar = list(title = "Pass")) |>
    plotly::layout(xaxis = list(title = LAB$item, type = "category"),
                   yaxis = list(title = LAB$wave_axis, type = "category")) |>
    kits_style(LAB$title_resp)
}
