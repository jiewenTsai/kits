# 來源紀錄（provenance）：呼叫、kits 版本、設定、資料的指紋。沒有時間戳記：同樣的輸入給完全相同的物件。

# 呼叫寫成一行文字（太長時截斷），不帶資料（隱私）。照使用者寫的樣子寫，只換掉資料：
# * 參數是名字（symbol）時照寫；
# * 參數是值（例如 do.call()、map() 給進來的矩陣或向量）時寫 <data>；
# * 參數是程式（call）時照寫，但給資料的參數裡的常數（數字、字串）換成 <data>。
# 給資料的參數：沒有預設值的參數、人的資料（DATA_ARGS）和 ...；其他有預設值的參數是設定，單一的常數照寫。
call_text <- function(call, name = NULL) {
  call <- tidy_call(call)
  if (!is.null(name)) call[[1]] <- as.symbol(name)   # S3 的方法（例如 kits_norms.default）寫成泛型的名字
  txt <- paste(trimws(deparse(call, width.cutoff = 500L, nlines = 5L)), collapse = " ")
  txt <- gsub("`<data>`", "<data>", txt, fixed = TRUE)
  if (nchar(txt) > 300) txt <- paste0(substr(txt, 1, 297), "...")
  txt
}

# 有預設值（NULL），但給的是人的資料的參數
DATA_ARGS <- c("persons", "wave_m", "age_m", "cohort", "sex", "d", "D", "se_from", "se_to", "theta_dist",
               "persons_from", "persons_to")

# 呼叫裡的資料換成 <data>（內部用，給 call_text()）
tidy_call <- function(call) {
  if (!is.call(call)) {
    return(as.symbol("<data>"))
  }
  fun <- call_function(call[[1]])
  call[[1]] <- fun$name
  formal <- arg_formals(call, fun$fun)
  settings <- if (is.null(fun$fun)) character() else setting_args(fun$fun)
  for (i in seq_along(call)[-1]) {
    if (is.null(call[[i]])) next
    call[[i]] <- scrub_arg(call[[i]], data = !formal[i - 1] %in% settings)
  }
  call
}

# 每個參數對到哪一個正式參數（沒有名字的照位置）：先把每個參數換成記號 .kits_arg_i，再用 match.call() 對
arg_formals <- function(call, fun) {
  n <- length(call) - 1
  out <- rep("", n)
  if (is.null(fun) || n == 0) {
    return(out)
  }
  marked <- call
  for (i in seq_len(n)) marked[[i + 1]] <- as.symbol(paste0(".kits_arg_", i))
  matched <- tryCatch(as.list(match.call(fun, marked, expand.dots = FALSE)), error = function(e) NULL)
  if (is.null(matched)) {
    return(out)
  }
  for (nm in setdiff(names(matched), c("", "..."))) {
    i <- match(as.character(matched[[nm]]), paste0(".kits_arg_", seq_len(n)))
    if (!is.na(i)) out[i] <- nm
  }
  out
}

# 呼叫的函數：名字和函數本身（do.call() 給的是函數本身時，找回 kits 裡的名字）
call_function <- function(f) {
  ns <- asNamespace("kits")
  if (is.symbol(f)) {
    name <- as.character(f)
    return(list(name = f, fun = if (exists(name, envir = ns, inherits = FALSE)) get(name, envir = ns) else NULL))
  }
  if (is.function(f)) {
    hit <- Filter(\(n) identical(get(n, envir = ns), f), ls(ns))
    if (length(hit) > 0) {
      return(list(name = as.symbol(hit[1]), fun = f))
    }
    return(list(name = as.symbol("<function>"), fun = f))
  }
  list(name = f, fun = NULL)
}

# 設定的參數：有預設值、不是人的資料（例如 n_boot = 0、scale = NULL、cut = 2）
setting_args <- function(fun) {
  f <- as.list(formals(fun))
  has_default <- map_lgl(f, \(v) !(is.symbol(v) && identical(as.character(v), "")))
  ages <- intersect("ages", names(f))   # 常模的月齡是設定（沒有預設值），不是孩子的資料
  setdiff(c(names(f)[has_default], ages), c("...", DATA_ARGS))
}

# 一個參數：名字照寫；程式照寫（給資料的參數裡的常數換成 <data>）；值換成 <data>。設定的單一常數和短的向量照寫
scrub_arg <- function(x, data) {
  if (is.symbol(x)) {
    return(x)
  }
  if (is.call(x)) {
    for (i in seq_along(x)[-1]) {
      if (!is.null(x[[i]])) x[[i]] <- scrub_arg(x[[i]], data)
    }
    return(x)
  }
  if (is.null(x) || (is.logical(x) && length(x) == 1 && is.null(attributes(x)))) {
    return(x)
  }
  if (!data && is.atomic(x) && is.null(attributes(x)) && length(x) <= 20) {
    return(x)
  }
  as.symbol("<data>")
}

# 一個 R 物件的 md5（base R）：serialize() 的位元組去掉檔頭（檔頭有 R 的版本），寫到暫存檔再算。
# 字串先轉成 UTF-8（enc2utf8()），所以 md5 只和內容有關，和字串的編碼標記無關。
hash_of <- function(x) {
  f <- tempfile()
  on.exit(unlink(f))
  writeBin(serialize(utf8_deep(x), NULL, version = 2L)[-(1:14)], f)
  unname(tools::md5sum(f))
}

# 每個字串（含 list 裡的、名字）轉成 UTF-8
utf8_deep <- function(x) {
  if (is.list(x)) {
    x[] <- lapply(x, utf8_deep)
  } else if (is.character(x)) {
    x[] <- enc2utf8(x)
  }
  if (!is.null(names(x))) names(x) <- enc2utf8(names(x))
  x
}

# 作答資料的指紋：人次數、題號、排序後 visit_id 的 md5、作答（visit_id、item_id、分數）的 md5
data_fingerprint <- function(resp) {
  sorted <- resp |>
    transmute(person = as.character(visit_id), item = as.character(item_id),
              score = as.integer(item_score)) |>
    arrange(person, item)   # C locale（arrange() 的預設），和 method = "radix" 相同，所以和 locale 無關
  persons <- unique(sorted$person)
  list(n_rows = length(persons), items = sort(unique(sorted$item), method = "radix"),
       persons_md5 = hash_of(persons), scores_md5 = hash_of(list(sorted$person, sorted$item, sorted$score)))
}

# 數值輸入（連結的點、常模的分數和月齡）的指紋
values_fingerprint <- function(...) {
  v <- list(...)
  list(n_rows = length(v[[1]]), values_md5 = hash_of(v))
}

# 校準資料的指紋：kits_item 的 info$data；門檻表、drift$fix 的 attr(, "source")
source_of <- function(x) {
  if (inherits(x, "kits")) x <- x$item
  if (inherits(x, "kits_item")) {
    return(x$info$data)
  }
  attr(x, "source")
}

make_info <- function(call, settings = list(), data = NULL, ...) {
  structure(c(list(call = call, version = as.character(utils::packageVersion("kits")), settings = settings, data = data),
              Filter(Negate(is.null), list(...))), class = "kits_info")
}

# 人的資料（persons 的欄、wave_m、sex）的指紋：照 visit_id 排序後的 md5
persons_fingerprint <- function(p) {
  columns <- p |>
    as.data.frame() |>
    arrange(as.character(visit_id)) |>   # C locale，和 method = "radix" 相同
    map(\(v) if (is.factor(v)) as.character(v) else v)
  list(n_rows = length(columns$visit_id), columns = setdiff(names(columns), "visit_id"),
       values_md5 = hash_of(columns))
}

# 輸入是 data.frame（寬或長的作答）或 kits_subj 時，記下它的名字和形式
# （記在 info 裡：怎麼從它拿到 X，或 D、wave_m、age_m）。名字只記 symbol；其他寫 <data>（不帶資料）
data_input <- function(x, expr) {
  if (!is.data.frame(x) || isTRUE(attr(x, "kits_long"))) {
    return(NULL)
  }
  name <- if (is.symbol(expr)) as.character(expr) else "<data>"
  if (inherits(x, "kits_subj")) {
    return(list(name = name, class = "kits_subj", columns = names(x)))
  }
  if (inherits(x, "kits_resp")) {
    return(list(name = name, class = "kits_resp", columns = names(x)))
  }
  list(name = name, class = if (is_long(x)) "long data.frame" else "wide data.frame", columns = names(x))
}

info_of <- function(x) {
  if (is.data.frame(x) || is.atomic(x)) attr(x, "info") else x$info
}


#' Provenance of a result
#'
#' Every main result of kits records how it was made: the call (object names only; data given by
#' value show as `<data>`), the kits version, the settings, and a fingerprint of the data (counts,
#' item ids, md5 hashes of the sorted ids and responses). The record holds no data and no time
#' stamp, so a result can be shared and the same input gives an identical object.
#'
#' @param x A result of kits, for example of [kits_item()], [kits_subj()], [kits_linking()] or
#'   [kits_norms()].
#' @return A `kits_info` list (printed briefly): `call`, `version`, `settings`, `data`, and, when
#'   they apply:
#'   * `reference`: the fixed parameters; `item_data`: the data the item parameters came from;
#'   * `persons`: a fingerprint of the person information (`persons` of [kits_subj()], or `wave_m`);
#'   * `input`: the name and form of a data.frame or [kits_subj()] given as data;
#'   * `conversions` (links applied by [kits_convert()]) and `links` (their `intercept`, `slope`).
#' @examples
#' betas <- list(a1 = c(-1, 0), a2 = c(0, 1), a3 = 0.5, a4 = -0.5)
#' X <- kits_sim(seq(-3, 3, length.out = 300), betas, seed = 1)
#' item <- kits_item(X, scale = "demo")
#' kits_info(item)
#' identical(kits_item(X, scale = "demo"), item)
#'
#' # Person information given to kits_subj() is fingerprinted too
#' s <- kits_subj(X, item, persons = data.frame(visit_id = rownames(X), wave_m = 12))
#' kits_info(s)
#' @export
kits_info <- function(x) {
  info <- info_of(x)
  abort_if(is.null(info), "`x` has no record of how it was made (no `info`). Give a result of a kits function.")
  info
}

#' @rdname kits_info
#' @param ... Not used.
#' @export
print.kits_info <- function(x, ...) {
  cat("<kits_info>\n")
  cat("  Call:", x$call, "\n")
  cat("  kits version:", x$version, "\n")
  if (length(x$settings) > 0) {
    vals <- map_chr(x$settings, \(v) if (is.null(v)) "NULL" else paste(format(v, trim = TRUE), collapse = ", "))
    cat("  Settings:", paste0(names(x$settings), " = ", vals, collapse = "; "), "\n")
  }
  if (!is.null(x$input)) cat(sprintf("  Input: %s (%s)\n", x$input$name, x$input$class))
  print_fingerprint("Data", x$data)
  if (!is.null(x$persons)) {
    cat(sprintf("  Persons table: %d person-waves, columns %s; md5 %s\n", x$persons$n_rows, paste(x$persons$columns,
        collapse = ", "),
                x$persons$values_md5))
  }
  if (!is.null(x$reference)) {
    cat("  Reference (fixed_params): scale", scale_label(x$reference$scale), "\n")
    print_fingerprint("  Reference data", x$reference$data)
  }
  print_fingerprint("Item parameters from data", x$item_data)
  if (length(x$conversions) > 0) cat("  Converted:", paste(x$conversions, collapse = "; "), "\n")
  invisible(x)
}

print_fingerprint <- function(label, fp) {
  if (is.null(fp)) {
    return(invisible())
  }
  if (!is.null(fp$scores_md5)) {
    items <- some_ids(fp$items, 6)
    cat(sprintf("  %s: %d person-waves, %d items (%s)\n    person ids md5 %s, responses md5 %s\n",
                label, fp$n_rows, length(fp$items), items, fp$persons_md5, fp$scores_md5))
  } else {
    cat(sprintf("  %s: %d values, md5 %s\n", label, fp$n_rows, fp$values_md5))
  }
}
