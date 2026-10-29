# 量尺：每個結果記得自己在哪一條量尺上（內部用的小工具，加上 kits_convert()）。
# * kits_item()：item$scale（門檻 thr 的量尺）；換算過的 kits_item 另外有 beta_scale（beta 沒換，還在校準的量尺）。
# * 表格（門檻表、kits_subj、kits_subj_table、drift$fix）：attr(, "scale")。
# * 連結：link$from_scale、link$to_scale。

# 一個物件的量尺標籤；沒有記錄時 NA
scale_of <- function(x) {
  if (inherits(x, "kits")) x <- x$item
  s <- if (inherits(x, "kits_item")) x$scale else attr(x, "scale")
  if (is.null(s) || length(s) == 0) NA_character_ else s
}

# beta 所在的量尺：換算過的 kits_item（或門檻表）門檻換了、beta 沒換
beta_scale_of <- function(x) {
  if (inherits(x, "kits")) x <- x$item
  b <- if (inherits(x, "kits_item")) x$beta_scale else attr(x, "beta_scale")
  if (is.null(b)) scale_of(x) else b
}

# 兩個量尺標籤都有記錄、又不同時停下來（kits_error_scale）；... 是錯誤訊息
need_same_scale <- function(a, b, ...) {
  abort_if(!is.na(a) && !is.na(b) && a != b, ..., class = "kits_error_scale")
}

# 印出來的量尺：沒有記錄時寫 "(not recorded)"
scale_label <- function(s) {
  if (length(s) == 0 || all(is.na(s))) "(not recorded)" else paste(ifelse(is.na(s), "(not recorded)", s), collapse = ", ")
}

# 預設的量尺標籤：資料物件的名字（symbol）。資料不是一個有名字的物件時（例如 do.call()、map()、X[keep, ]），
# 標籤是 "scale 1"，並說一聲（標籤不帶資料；兩個這樣的校準會有同樣的標籤，連結時會停下來）
default_scale <- function(expr, data = NULL) {
  label <- resp_scale(data)
  if (!is.null(label)) {
    kits_inform("No `scale` given: the scale is labelled \"", label,
                "\" (the domain and cohort of the responses).",
                class = "kits_message_scale")
    return(label)
  }
  if (is.symbol(expr)) {
    return(as.character(expr))
  }
  kits_inform("No `scale` given and the data are not a named object, so the scale is labelled \"scale 1\". ",
              "Give `scale`, for example scale = \"lan M03\".", class = "kits_message_scale")
  "scale 1"
}

need_scale_arg <- function(scale, arg = "scale") {
  abort_if(!is.character(scale) || length(scale) != 1 || is.na(scale) || scale == "",
           "`", arg, "` must be one label, for example \"lan M03\".")
}

# 設定量尺（kits_item 連門檻表一起）
set_scale <- function(x, scale) {
  if (inherits(x, "kits_item")) {
    x$scale <- scale
    attr(x$thresholds, "scale") <- scale
    return(x)
  }
  attr(x, "scale") <- scale
  x
}


#' Convert scores or thresholds with one or more links
#'
#' Applies the links in order, each a straight line: new value = intercept + slope x old value, and
#' new SE = slope x old SE. Before each link the scale of `x` must be one the link starts from
#' (`link$from_scale`), else it stops (class `kits_error_scale`); an unrecorded scale warns.
#' When every link has bootstrap draws (`n_boot > 0`), a table also gets `se_link` (linking error
#' only) and `se_total = sqrt(se^2 + se_link^2)`.
#'
#' Details: `kits_algorithm("wle")` (the conversion) and tutorial_algorithms/02_scoring.html.
#'
#' @param x Scores or thresholds: a [kits_subj()] result or a data.frame with `theta` (and `se`);
#'   a [kits_subj_table()]; a [kits_item()] or its `coef()` (`thr`, `thr_se`; `beta` is not
#'   converted); or a numeric vector of theta (give `from`). Other columns are kept.
#' @param ... One or more links ([kits_linking()], [kits_link_steps]), applied in this order.
#' @param from The scale of `x` when `x` does not record it (for example a numeric vector).
#' @return `x` on the scale where the last link ends (recorded in `attr(, "scale")`), with the same
#'   column names; with draws also `se_link` and `se_total`.
#' @examples
#' betas <- list(a1 = -2, a2 = -1, a3 = c(-0.5, 0.5), a4 = 0, a5 = c(0.5, 1.5), a6 = 2, a7 = 3)
#' X03 <- kits_sim(rnorm(500, 0), betas[1:5], seed = 1)
#' X36 <- kits_sim(rnorm(500, 2), betas[3:7], prefix = "q", seed = 2)
#' item03 <- kits_item(X03, scale = "lan M03")
#' item36 <- kits_item(X36, scale = "lan M36")
#' co <- kits_linking(item36, item03, n_boot = 0)
#' anchors <- data.frame(item_id = c("a1", "a2", "a3", "a4", "a5"), tau = c(40, 47, 51, 55, 60))
#' d_link <- kits_linking(item03, anchors, n_boot = 0)
#'
#' s36 <- kits_subj(X36, item36)
#'
#' # One link per line: M36 -> M03, then M03 -> D
#' D36 <- s36 |>
#'   kits_convert(co) |>
#'   kits_convert(d_link)
#' head(D36)
#'
#' # The same in one call
#' all.equal(D36, kits_convert(s36, co, d_link))
#'
#' # A numeric vector: give its scale
#' c(-1, 0, 1) |>
#'   kits_convert(d_link, from = "lan M03")
#'
#' # Chosen theta values with their linking error (links made with n_boot > 0)
#' co_boot <- kits_linking(item36, item03, n_boot = 100)
#' d_boot <- kits_linking(item03, anchors, n_boot = 100)
#' data.frame(theta = c(-1, 0, 1)) |>
#'   kits_convert(co_boot, d_boot, from = "lan M36")
#'
#' # The wrong order stops
#' try(s36 |> kits_convert(d_link))
#' @export
kits_convert <- function(x, ..., from = NULL) {
  convert_links(x, list(...), from)
}

# kits_convert() 的本體：一串連結（list），照順序套用並檢查量尺
convert_links <- function(x, links, from = NULL) {
  abort_if(length(links) == 0, "Give at least one link after `x`, e.g. kits_convert(s36, co, d_link).")
  is_link <- map_lgl(links, \(l) inherits(l, "kits_link"))
  abort_if(!all(is_link), "Argument ", which(!is_link)[1] + 1, " is not a link: give results of kits_linking() ",
           "or kits_link_line() after `x`.")
  abort_if(inherits(x, "kits"), "`x` is a kits(): convert its scores or its items, e.g. kits_convert(k$subj, link).")
  abort_if(!is.numeric(x) && !is.data.frame(x) && !inherits(x, "kits_item"),
           "`x` must be scores (kits_subj(), a data.frame with theta), thresholds (kits_item(), coef()) ",
           "or a numeric vector.")
  abort_if(is.data.frame(x) && !any(c("theta", "thr") %in% names(x)), "`x` has no `theta` or `thr` column to convert.")
  scale <- scale_of(x)   # 數字向量：前一個 kits_convert() 記下的量尺（管線裡下一條連結照它檢查）
  if (!is.null(from)) {
    need_scale_arg(from, "from")
    need_same_scale(scale, from, sprintf("`x` records scale %s, but `from` is %s.", scale, from))
    scale <- from
  }
  # 門檻換算時，記下 beta 的量尺（beta 不換）。kits_item 的門檻表也記下，所以 coef() 拿出來計分時量尺還是對的
  if (inherits(x, "kits_item")) {
    x$beta_scale <- beta_scale_of(x)
    attr(x$thresholds, "beta_scale") <- x$beta_scale
  } else if (is.data.frame(x) && "thr" %in% names(x)) {
    attr(x, "beta_scale") <- beta_scale_of(x)
  }

  # 連結誤差（se_link）：從換算前的值算；x 已經有 se_link 時，接著算（它是前面的連結的誤差）
  x0 <- NULL
  v0 <- 0
  had_se <- FALSE
  if (is.data.frame(x)) {
    x0 <- if ("theta" %in% names(x)) x$theta else x$thr
    had_se <- "se_link" %in% names(x)
    if (had_se) v0 <- ifelse(is.finite(x$se_link), x$se_link^2, NA_real_)
  }
  # 前面的換算裡，誤差沒有帶下來的連結（沒有重抽，或換算的是 kits_item）：之後都不加 se_link
  earlier <- no_se_record(x)

  chain <- links |>
    map_chr(\(l) sprintf("%s -> %s", scale_label(l$from_scale), scale_label(l$to_scale))) |>
    paste(collapse = ", then ")
  done <- character()
  lines <- list()   # 套用的線：from、to、a、b（記在結果裡）
  for (i in seq_along(links)) {
    link <- links[[i]]
    from_scale <- if (is.null(link$from_scale)) NA_character_ else link$from_scale
    to_scale <- if (is.null(link$to_scale)) NA_character_ else link$to_scale
    if (is.na(scale) || all(is.na(from_scale))) {
      kits_warn(sprintf("Link %d (%s -> %s): scale not recorded, not checked.", i, scale_label(from_scale),
                        scale_label(to_scale)),
                class = "kits_warning_scale")
    } else if (!scale %in% from_scale) {
      kits_abort(sprintf("`x` is on scale %s, but link %d converts %s -> %s. ", scale, i, scale_label(from_scale),
                         scale_label(to_scale)),
                 sprintf("The links given: %s. Give the links in order, starting with a link from %s.", chain, scale),
                 class = "kits_error_scale")
    }
    x <- convert_one(x, link$intercept, link$slope)
    done <- c(done, sprintf("%s -> %s", scale_label(scale), scale_label(to_scale)))
    lines[[i]] <- data.frame(from = scale_label(scale), to = scale_label(to_scale),
                             intercept = link$intercept, slope = link$slope)
    scale <- to_scale
  }
  lines <- list_rbind(lines)
  if (is.numeric(x)) {
    attr(x, "scale") <- scale
    return(x)
  }
  x <- add_link_se(x, x0, v0, links, earlier)
  no_se_link_message(x, links, had_se, earlier)
  x <- set_no_se_record(x, bind_rows(earlier, new_no_se_record(x, links, done)))
  if (inherits(x, "kits_item")) {
    x$info$conversions <- c(x$info$conversions, done)
    x$info$links <- bind_rows(x$info$links, lines)
  } else if (!is.null(attr(x, "info"))) {
    info <- attr(x, "info")
    info$conversions <- c(info$conversions, done)
    info$links <- bind_rows(info$links, lines)
    attr(x, "info") <- info
  }
  x <- record_conversion(x, links, from)
  set_scale(x, scale)
}

# 分數表記下這次的換算（給 kits_subj_table(scores)：原始分 -> theta 的對照照同樣的連結換）
record_conversion <- function(x, links, from) {
  sources <- attr(x, "raw_theta")
  if (is.null(sources)) {
    return(x)
  }
  step <- list(links = links, from = from)
  attr(x, "raw_theta") <- map(sources, \(src) {
    src$steps <- c(src$steps, list(step))
    src
  })
  x
}

# 連結誤差沒有帶下來的連結的紀錄（attr(, "links_without_se")）：一張表，link（"from -> to"）和 reason。
# reason："no draws"（連結沒有重抽），"kits_item"（連結有重抽，但換算的是 kits_item，它不帶 se_link）。
# kits_item 的紀錄放在門檻表上，所以 coef(item) 拿出來再換算時也看得到。
# 有紀錄的表再換算時不加 se_link：否則 se_link 只含後面的連結，低估了連結誤差。
no_se_record <- function(x) {
  if (inherits(x, "kits_item")) x <- x$thresholds
  rec <- attr(x, "links_without_se")
  if (is.null(rec)) empty_no_se_record() else rec
}

empty_no_se_record <- function() {
  data.frame(link = character(), reason = character())
}

# 這次換算新增的紀錄：沒有重抽的連結；換算 kits_item 時，每一條連結
new_no_se_record <- function(x, links, done) {
  has_draws <- map_lgl(links, link_has_draws)
  keep <- !has_draws | inherits(x, "kits_item")
  reason <- ifelse(has_draws, "kits_item", "no draws")
  data.frame(link = done[keep], reason = reason[keep])
}

set_no_se_record <- function(x, rec) {
  rec <- if (nrow(rec) == 0) NULL else unique(rec)
  if (inherits(x, "kits_item")) {
    attr(x$thresholds, "links_without_se") <- rec
    return(x)
  }
  attr(x, "links_without_se") <- rec
  x
}

# 幾張表的紀錄合在一起（rbind()、bind_rows()）
combine_no_se_records <- function(parts) {
  recs <- map(parts, \(p) if (is.data.frame(p)) attr(p, "links_without_se"))
  rec <- unique(bind_rows(compact(recs)))
  if (nrow(rec) == 0) NULL else rec
}

# 連結誤差沒有帶下來，所以 se_link 不加（或拿掉）時，說一聲是哪幾條連結。
# 整串連結（前面的換算加上這次的）都沒有重抽時不說，和一次給所有連結相同。
no_se_link_message <- function(x, links, had_se, earlier) {
  if (!is.data.frame(x) || "se_link" %in% names(x)) {
    return(invisible())
  }
  has_draws <- map_lgl(links, link_has_draws)
  earlier_draws <- any(earlier$reason == "kits_item")
  if (!any(has_draws) && !had_se && !earlier_draws) {
    return(invisible())
  }
  parts <- c(this_call_no_draws_text(has_draws), earlier_no_se_text(earlier))
  kits_inform("No se_link column: ", paste(parts, collapse = "; "),
              ". Make every link with n_boot > 0 (e.g. 200) to get the linking error.", class = "kits_message_scale")
}

this_call_no_draws_text <- function(has_draws) {
  if (all(has_draws)) {
    return(character())
  }
  paste0("link(s) ", paste(which(!has_draws), collapse = ", "), " of this conversion have no bootstrap draws (n_boot = 0)")
}

earlier_no_se_text <- function(earlier) {
  if (nrow(earlier) == 0) {
    return(character())
  }
  item_text <- "applied to a kits_item(), which keeps no se_link; convert coef(item) instead"
  why <- ifelse(earlier$reason == "kits_item", item_text, "no bootstrap draws, n_boot = 0")
  paste0("an earlier conversion applied link(s) whose linking error is not kept: ",
         paste0(earlier$link, " (", why, ")", collapse = ", "))
}

# 用一串連結（list）換算；和 kits_convert(x, link_1, link_2, ..., from = from) 相同
convert_by_links <- function(x, links, from = NULL) {
  if (length(links) == 0) {
    return(x)
  }
  args <- c(list(x), links, list(from = from))
  as.vector(do.call(kits_convert, args))
}

# 一條線：new = intercept + slope old；SE 乘 slope。按類別分派（內部的 S3）：數字向量、kits_item（換它的門檻表）、表
convert_one <- function(x, intercept, slope) UseMethod("convert_one")

#' @exportS3Method
convert_one.default <- function(x, intercept, slope) {
  line_apply(list(intercept = intercept, slope = slope), as.vector(x))
}

#' @exportS3Method
convert_one.kits_item <- function(x, intercept, slope) {
  x$thresholds <- convert_one(x$thresholds, intercept, slope)
  x
}

# 表：theta 和 se（分數），thr 和 thr_se（門檻）
#' @exportS3Method
convert_one.data.frame <- function(x, intercept, slope) {
  line <- list(intercept = intercept, slope = slope)
  if ("theta" %in% names(x)) {
    x$theta <- line_apply(line, x$theta)
    if ("se" %in% names(x)) x$se <- line_stretch(line, x$se)
  }
  if ("thr" %in% names(x)) {
    x$thr <- line_apply(line, x$thr)
    if ("thr_se" %in% names(x)) x$thr_se <- line_stretch(line, x$thr_se)
  }
  x
}
