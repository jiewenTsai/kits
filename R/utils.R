# 內部用的小工具（不匯出）。整理資料用 dplyr、purrr（tidyverse 的寫法）；檢查的規則和 kitsPreview 的 utils.R 相同。

# 檢查欄位。缺的時候停下來，說清楚缺哪幾欄。
need_cols <- function(data, cols, arg) {
  abort_if(!is.data.frame(data),
           sprintf("`%s` must be a data.frame with column(s): %s.", arg, paste(cols, collapse = ", ")))
  miss <- setdiff(cols, names(data))
  abort_if(length(miss) > 0,
           sprintf("`%s` is missing column(s): %s. It needs: %s.",
                   arg, paste(miss, collapse = ", "), paste(cols, collapse = ", ")))
}

# 用固定的亂數種子跑，亂數的種類也固定（R 的預設：Mersenne-Twister、Inversion、Rejection），
# 所以同一個 seed 在任何工作階段給同樣的結果。跑完還原使用者原本的亂數種類和狀態。
with_seed <- function(seed, code) {
  had_seed <- exists(".Random.seed", envir = globalenv())
  old <- if (had_seed) get(".Random.seed", envir = globalenv())
  old_kind <- RNGkind()
  on.exit({
    suppressWarnings(RNGkind(old_kind[1], old_kind[2], old_kind[3]))
    if (had_seed) assign(".Random.seed", old, envir = globalenv()) else rm(".Random.seed", envir = globalenv())
  })
  set.seed(seed, kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
  code
}

# 常模表依 segment 分段，依月齡由小到大排。沒有 segment 欄時就是一段。
norm_segments <- function(norm) {
  norm <- as.data.frame(norm)
  segments <- if ("segment" %in% names(norm)) split(norm, norm$segment) else list(all = norm)
  segments <- map(segments, \(s) arrange(s, age_m))
  segments[order(map_dbl(segments, \(s) min(s$age_m)))]
}

# 照 base 的 sort() 排（目前的 locale；和 tapply()、split()、merge() 的順序相同），給 arrange() 用。
# arrange() 自己排文字時用 C locale；用這個，表格的列和 0.14.0 的順序相同（例如連結的點，重抽的結果跟著點的順序）。
sort_rank <- function(x) match(x, sort(unique(x)))

# 作答資料（data）轉成長格式（visit_id、item_id、item_score），並檢查一次（內部用，每個收 data 的函數都用）。
# 收三種：
# * 矩陣：一列一個人次、一欄一題；列名是 visit_id，欄名是題號；NA 是沒施測。
# * 寬的 data.frame：同樣的形狀；visit_id 是叫 visit_id 的欄，沒有這一欄時是列名；其他欄都是題。
# * 長的 data.frame：visit_id、item_id、item_score（有 item_id 或 item_score 欄就當成長格式，少一欄時說清楚）。
# 檢查（錯誤的類別 kits_error_data）：列名（人次）、題號不重複；分數是 0 以上的整數；
# gaps = TRUE 時，每題觀察到的分數中間沒有缺（例如 0 和 2，沒有 1）。
# 沒有作答的題和人次拿掉，有訊息（kits_message_data）。
resp_long <- function(data, arg = "data", gaps = TRUE) {
  if (is.data.frame(data) && isTRUE(attr(data, "kits_long"))) {
    return(data)
  }
  data <- old_ids(data, arg)
  if (is_long(data)) {
    resp <- long_resp(data, arg)
  } else {
    resp <- matrix_resp(wide_matrix(data, arg), arg)
  }
  if (gaps) check_gaps(resp, arg)
  resp
}

# 0.22.0 以前的欄名：person_id 是 visit_id、pid 是 subj_id。這一版照常認得，說一聲；下一版拿掉。
# quiet = TRUE：同一份資料第二次讀時不再說
old_ids <- function(x, arg = "data", quiet = FALSE) {
  if (!is.data.frame(x)) {
    return(x)
  }
  old <- c(person_id = "visit_id", pid = "subj_id")
  found <- names(old)[names(old) %in% names(x) & !old %in% names(x)]
  if (length(found) == 0) {
    return(x)
  }
  if (!quiet) {
    kits_inform(sprintf("`%s` has the old column name(s) %s; kits 0.23.0 calls them %s. ", arg,
                        paste(found, collapse = ", "), paste(old[found], collapse = ", ")),
                "They are renamed for now; rename them in your data.", class = "kits_message_data")
  }
  names(x)[match(found, names(x))] <- old[found]
  x
}

# 有 item_id 和 item_score 欄的 data.frame 是長格式（只有其中一欄時也當成長格式，缺的欄在檢查時說清楚）
is_long <- function(x) is.data.frame(x) && any(c("item_id", "item_score") %in% names(x))

# 寬的 data.frame 轉成矩陣：visit_id 欄（或列名）是列名，其他欄是題
wide_matrix <- function(x, arg) {
  if (is.matrix(x)) {
    return(x)
  }
  abort_if(is.atomic(x) && is.null(dim(x)),
           sprintf("`%s` is a vector, not a matrix. To keep one row as a matrix, use X[1, , drop = FALSE].", arg))
  abort_if(!is.data.frame(x),
           sprintf("`%s` must be a matrix or a wide data.frame (rows = person-waves, columns = items), ", arg),
           "or a long data.frame with visit_id, item_id, item_score.")
  if (inherits(x, "kits_resp")) x <- x[c("visit_id", resp_items(x))]   # 作答表：只拿題目的欄
  x <- as.data.frame(x)
  no_ids <- !"visit_id" %in% names(x) && .row_names_info(x) < 0   # 沒有編號：列名是自動的 1, 2, ...
  if (no_ids) numbered_message(nrow(x), arg)
  ids <- if ("visit_id" %in% names(x)) as.character(x$visit_id) else rownames(x)
  items <- x[setdiff(names(x), "visit_id")]
  person_cols <- intersect(names(items), c("wave_m", "age_m", "sex", "cohort"))
  abort_if(length(person_cols) > 0,
           sprintf(paste0("`%s` has person column(s) %s: a wide data.frame holds only ",
                          "visit_id and the items. "), arg,
                   paste(person_cols, collapse = ", ")),
           "Give the person information separately (for example `wave_m = `, or `persons = ` of kits_subj()).")
  ok <- map_lgl(items, \(v) is.numeric(v) || is.logical(v))
  abort_if(!all(ok),
           sprintf("`%s` has item column(s) that are not numbers: %s. ", arg, paste(utils::head(names(items)[!ok], 10),
                                                                                    collapse = ", ")),
           paste0("Item scores must be 0, 1, 2, ... (NA = not administered); the person ids ",
                  "go in a column named visit_id."))
  X <- as.matrix(items)
  rownames(X) <- ids
  X
}

# 矩陣（或寬的 data.frame）的檢查和長格式
matrix_resp <- function(X, arg) {
  abort_if(nrow(X) == 0 || ncol(X) == 0, sprintf("`%s` has no rows or no columns.", arg), class = "kits_error_data")
  if (is.null(rownames(X)) && !is.null(colnames(X))) {   # 沒有編號：像 TAM 一樣自動編 1, 2, 3, ...
    numbered_message(nrow(X), arg)
    rownames(X) <- as.character(seq_len(nrow(X)))
  }
  if (is.null(rownames(X)) || is.null(colnames(X)) || anyNA(rownames(X)) || anyNA(colnames(X)) ||
      any(rownames(X) == "") || any(colnames(X) == "")) {
    kits_abort(sprintf(paste0("`%s` needs person ids (row names, or a visit_id column) and item ",
                              "ids (column names), none empty."), arg),
               class = "kits_error_data")
  }
  if (anyDuplicated(rownames(X))) {
    dup <- unique(rownames(X)[duplicated(rownames(X))])
    kits_abort(sprintf("`%s` has %d repeated person id(s)%s. ", arg, length(dup), ids_text(dup)),
               "Each row is one person-wave, so use an id per person-wave, e.g. paste(child_id, wave_m). ",
               "Data stacked from several kits_sim() calls: give each call its own `prefix`.", class = "kits_error_data")
  }
  abort_if(anyDuplicated(colnames(X)),
           sprintf("`%s` has repeated item ids (column names).", arg), class = "kits_error_data")
  check_scores(X, arg)
  idx <- which(!is.na(X), arr.ind = TRUE)
  abort_if(nrow(idx) == 0, sprintf("`%s` has no responses (all NA).", arg), class = "kits_error_data")
  idx <- idx[order(idx[, 1], idx[, 2]), , drop = FALSE]
  resp <- data.frame(visit_id = rownames(X)[idx[, 1]], item_id = colnames(X)[idx[, 2]],
                     item_score = as.integer(X[idx]))
  drop_message(setdiff(colnames(X), resp$item_id), setdiff(rownames(X), resp$visit_id))
  resp
}

# 沒有人次的編號時說一聲：照列的順序編 1 到 n
numbered_message <- function(n, arg) {
  kits_inform(sprintf("`%s` has no person ids (row names or a visit_id column): ", arg),
              sprintf("the rows are numbered 1 to %d, in their order.", n),
              class = "kits_message_data")
}

# 長格式的檢查
long_resp <- function(x, arg) {
  need_cols(x, c("visit_id", "item_id", "item_score"), arg)
  x <- as.data.frame(x)
  all_items <- unique(as.character(x$item_id))
  all_persons <- unique(as.character(x$visit_id))
  resp <- x |>
    filter(!is.na(item_score)) |>
    select(visit_id, item_id, item_score)
  abort_if(nrow(resp) == 0, sprintf("`%s` has no responses.", arg), class = "kits_error_data")
  check_scores(resp$item_score, arg)
  resp <- resp |>
    mutate(visit_id = as.character(visit_id), item_id = as.character(item_id),
           item_score = as.integer(item_score))
  empty <- function(v) any(is.na(v) | v == "")
  abort_if(empty(resp$visit_id) || empty(resp$item_id),
           sprintf("`%s` has empty visit_id or item_id.", arg), class = "kits_error_data")
  abort_if(anyDuplicated(resp[c("visit_id", "item_id")]),
           sprintf("`%s` has a person with the same item more than once. ", arg),
           "Each person-wave needs its own visit_id, e.g. paste(child_id, wave_m).", class = "kits_error_data")
  drop_message(setdiff(all_items, resp$item_id), setdiff(all_persons, resp$visit_id))
  rownames(resp) <- NULL
  resp
}

# 拿掉的題和人次（沒有任何作答）。人次只給人數（隱私），見 ids_text()
drop_message <- function(items, persons) {
  if (length(items) > 0) {
    kits_inform(length(items), " item(s) without responses are dropped: ", some_ids(items, 10), ".",
                class = "kits_message_data")
  }
  if (length(persons) > 0) {
    kits_inform(length(persons), " person-wave(s) without responses are dropped", ids_text(persons), ".",
                class = "kits_message_data")
  }
}

# 訊息和錯誤裡的 visit_id（隱私）：預設不列出，只說怎麼看；options(kits.show_ids = TRUE) 時列出前 n 個
ids_text <- function(ids, n = 5) {
  if (!isTRUE(getOption("kits.show_ids", FALSE))) {
    return(" (ids not shown; options(kits.show_ids = TRUE) shows them)")
  }
  paste0(" (for example ", some_ids(ids, n), ")")
}

# 前 k 個（題號、編號等），再多時加 ...；sep 是分隔（預設 ", "）
some_ids <- function(ids, k, sep = ", ") {
  paste0(paste(utils::head(ids, k), collapse = sep), if (length(ids) > k) paste0(sep, "...") else "")
}

# 每題觀察到的分數中間不能有缺（例如 0 和 2，沒有 1）
check_gaps <- function(resp, arg) {
  gaps <- resp |>
    distinct(item_id, item_score) |>
    summarise(gap = has_gap(item_score), scores = paste(sort(item_score), collapse = ", "), .by = item_id) |>
    filter(gap) |>
    arrange(sort_rank(item_id))
  abort_if(nrow(gaps) > 0, "Item(s) with a gap in their observed scores: ",
           paste0(gaps$item_id, " (", gaps$scores, ")", collapse = "; "),
           ". Scores must be consecutive; recode first.", class = "kits_error_data")
}

# 一題觀察到的分數中間有缺（例如 0、1、3）：summarise() 裡每題用（沒有列時 dplyr 用空的向量試算：沒有缺）
has_gap <- function(score) length(score) > 0 && !all(seq(min(score), max(score)) %in% score)

# 一題能自由校準：有 0 分、至少兩種分數（CML 的條件）；summarise() 裡每題用
can_calibrate <- function(score) any(score == 0) && n_distinct(score) > 1

# 沒有參數的題停下來（kits_error_items）；what 是訊息的開頭，後面列出前 10 題
need_params <- function(items, param_items, what, end = ".") {
  miss <- setdiff(items, param_items)
  abort_if(length(miss) > 0, what, some_ids(miss, 10), end, class = "kits_error_items")
}

# 作答分數：0, 1, 2, ... 的整數（NA 是沒施測）
check_scores <- function(x, arg) {
  if (!is.numeric(x) && !is.logical(x)) {
    what <- if (is.factor(x)) "a factor" else typeof(x)
    kits_abort(sprintf("`%s` must hold numbers (item scores 0, 1, 2, ...), not %s.", arg, what),
               if (is.factor(x)) " Convert it first, e.g. item_score = as.integer(as.character(item_score))." else "",
               class = "kits_error_data")
  }
  x <- x[!is.na(x)]
  abort_if(any(!is.finite(x) | x < 0 | x != round(x)),
           sprintf("`%s` has scores that are negative, infinite or not whole numbers.", arg), class = "kits_error_data")
}

# 已經轉成長格式的資料，標記起來，下一個函數不用再轉
resp_mark <- function(resp) {
  attr(resp, "kits_long") <- TRUE
  resp
}

# data 每一列的 visit_id（矩陣和寬的 data.frame：一列一個人次；長格式：一列一個作答）
row_ids <- function(data) {
  if (is.matrix(data)) {
    return(rownames(data))
  }
  data <- old_ids(data, quiet = TRUE)
  if (is_long(data) || "visit_id" %in% names(data)) as.character(data$visit_id) else rownames(data)
}

# 每個人次一個值（例如 wave_m、sex），轉成長格式 resp 每一列的值。
# 收兩種：有名字時，名字是 visit_id；沒有名字時，data 的每一列一個值（同樣的順序；長格式時同一個人要相同）。
per_person <- function(x, data, resp, arg) {
  if (!is.null(names(x))) {
    if (anyDuplicated(names(x))) {
      dup <- unique(names(x)[duplicated(names(x))])
      kits_abort(sprintf("`%s` has %d repeated name(s)%s.", arg, length(dup), ids_text(dup)), class = "kits_error_input")
    }
    no_value <- setdiff(resp$visit_id, names(x))
    abort_if(length(no_value) > 0,
             sprintf("`%s` has no value for %d person-wave(s) of the data%s.", arg, length(no_value), ids_text(no_value)))
    by_person <- data.frame(visit_id = names(x), value = unname(x))
  } else {
    ids <- row_ids(data)
    abort_if(length(x) != length(ids),
             sprintf(paste0("`%s` has %d value(s), but the data have %d row(s). Give one ",
                            "value per row, or name it by visit_id."),
                     arg, length(x), length(ids)))
    by_row <- data.frame(visit_id = ids, value = x)
    if (is_long(data)) {
      mixed <- by_row |>
        summarise(n_values = n_distinct(value), .by = visit_id) |>
        filter(n_values > 1) |>
        arrange(sort_rank(visit_id))
      abort_if(nrow(mixed) > 0,
               sprintf("`%s` has more than one value for %d person-wave(s)%s.", arg,
                       nrow(mixed), ids_text(mixed$visit_id)))
    }
    by_person <- distinct(by_row, visit_id, .keep_all = TRUE)
  }
  resp |>
    select(visit_id) |>
    left_join(by_person, by = "visit_id") |>
    pull(value)
}

# 一個數字的參數，不合時停下來（kits_error_input）："`arg` must be what."。
# lower、upper 是含的界限，above、below 是不含的界限；whole = TRUE：整數；finite = FALSE：Inf 也可以（NA 都不行）；
# not：不可以的值；null = TRUE：NULL 也可以
need_number <- function(x, arg, what, lower = -Inf, upper = Inf, above = NULL, below = NULL, whole = FALSE,
                        finite = TRUE, not = NULL, null = FALSE) {
  if (null && is.null(x)) {
    return(invisible())
  }
  ok <- is.numeric(x) && length(x) == 1 && !is.na(x) && (!finite || is.finite(x)) && x >= lower && x <= upper &&
    (is.null(above) || x > above) && (is.null(below) || x < below) && (!whole || x == round(x)) && !x %in% not
  if (!ok) kits_abort(sprintf("`%s` must be %s.", arg, what), class = "kits_error_input")
}

# theta：一個或多個有限的數字（沒給也停下來）；where 接在訊息後面
need_theta <- function(theta, where = "") {
  if (missing(theta) || !is.numeric(theta) || length(theta) == 0 || any(!is.finite(theta))) {
    kits_abort("`theta` must be one or more finite numbers", where, ".", class = "kits_error_input")
  }
}

# 報表的小數位數：只在印出來的時候限縮，物件裡的數值不變
need_digits <- function(digits) {
  need_number(digits, "digits", "one whole number from 0 to 15", lower = 0, upper = 15, whole = TRUE)
}

# 一個數字印成小數點下 digits 位
fmt <- function(x, digits) {
  x <- round(x, digits)
  x[!is.na(x) & x == 0] <- 0   # 不印 -0.00
  sprintf("%.*f", as.integer(digits), x)
}

# 表格裡不是整數的數值欄，印成小數點下 digits 位
format_table <- function(d, digits) {
  need_digits(digits)
  d <- as.data.frame(d)
  class(d) <- "data.frame"
  for (j in seq_along(d)) {
    v <- d[[j]]
    if (is.double(v) && any(!is.na(v) & is.finite(v) & v != round(v))) d[[j]] <- fmt(v, digits)
  }
  d
}

# kits 回傳的表格：還是 data.frame，只是印出來時小數點下 2 位
kits_table <- function(d, cls = NULL) {
  d <- as.data.frame(d)   # tibble 也變回 data.frame：kits 的表都是 data.frame
  attr(d, "kits_long") <- NULL   # dplyr 的動詞會帶著輸入的屬性；長格式的標記不屬於結果
  class(d) <- c(cls, "kits_table", "data.frame")
  d
}

# dexter 給題本的名字：bk 加上編號，位數是 ceiling(log10(題本數 + 1))（長格式資料的規則）
booklet_names <- function(n) {
  sprintf(paste0("bk%0", ceiling(log10(n + 1)), "i"), seq_len(n))
}

# 每個人次答的題目組合，用一個字串當鍵（題號裡不會有的分隔字元）
set_key <- function(items) {
  paste(sort(items), collapse = "\u001f")
}
set_items <- function(key) {
  strsplit(key, "\u001f", fixed = TRUE)[[1]]
}

# 每個人次（照 resp 裡出現的順序）：答的題目組合（set）、原始分；給了題目參數時再加上滿分（raw_max）
person_sets <- function(resp, parms = NULL) {
  if (!is.null(parms)) {
    item_max <- summarise(parms, item_max = max(item_score), .by = item_id)
    resp <- left_join(resp, item_max, by = "item_id")
  }
  summarise(resp, set = set_key(item_id), raw = sum(item_score),
            raw_max = if (!is.null(parms)) sum(item_max), .by = visit_id)   # NULL：沒有這一欄
}

# 作答組合（booklet）的編號和 dexter 相同：人次照 visit_id 排序（C locale，arrange() 的預設），作答組合照出現的
# 順序編號（bk1, bk2, ...，booklet_names()）。一個作答組合一列：set、n（人次）、first_person（第一個人）、pattern_id
number_booklets <- function(persons) {
  persons |>
    arrange(visit_id) |>
    summarise(n = n(), first_person = first(visit_id), .by = set) |>
    mutate(pattern_id = booklet_names(n()))
}
