# 題目表（kits_item_table）：一題一列，D、通過月齡、常模中位數到達的月齡、DIF、題幹；表裡帶著常模，
# 所以 plot() 直接畫成長曲線（type = "growth"）和成長階梯（type = "ladder"）。

#' Item table: each item with its D, ages, DIF and stem, for the report plots
#'
#' One row per item: the D of its milestone step, the age at which a share of the children pass it,
#' the age at which the norm median reaches its D, its sex DIF and its stem. The table keeps the
#' norm (summary numbers, no child's data), so `plot()` draws the growth curves with the items
#' (`type = "growth"`, one panel per domain) or the developmental ladder (`type = "ladder"`).
#'
#' Details: `vignette("kits-advanced", package = "kits")`, sections 7.3 to 7.5 and 10.9.
#'
#' @param thr The thresholds, on the scale of the norm: a [kits_item()] converted with
#'   [kits_convert()], or its thresholds table (`coef()`, `item_id`, `item_score`, `thr`, optionally
#'   `domain` and `column`), or a list of them (stacked in order).
#' @param norms The norm: a [kits_norms()] result for one domain, or by domain (`by = "domain"`);
#'   P10, P50 and P90 are used.
#' @param data Optional responses with ages: a [kits_resp()], or what it takes (a long data.frame,
#'   or a wide data.frame of item scores), with `age_m` and `wave_m`. For each item, a logistic
#'   regression of passing (a score of at least `step`) on age gives the age at which a share `rp`
#'   pass (NA when the pass rate of the waves does not cross `rp`).
#' @param dif Optional: a [kits_dif_combine()] or a [kits_dif()] result.
#' @param stems Optional: the item stems, a data.frame with `item_id` and `stem`, or a character
#'   vector named by item.
#' @param rp The response probability of `thr` and of the ages (default 0.5). A [kits_item()] gives
#'   its thresholds at this RP (`coef(item, rp = rp)`); a table must hold thresholds at this RP (its
#'   `attr(, "rp")`, RP50 when not recorded), otherwise the function stops.
#' @param dif_cut With one [kits_dif()] as `dif`: the size of DIF (logits) beyond which an item
#'   favors a group (default 0.3, as in [kits_dif_combine()]). Not used for a [kits_dif_combine()]
#'   result, whose `DIF` column is already decided.
#' @param step The pass step: a score of at least `step` passes, and `D` is the threshold of that
#'   step. One whole number for every item (default 1, the milestone step), or a vector named by
#'   item (items not named use 1), e.g. `c(cogD01 = 3)` for a reversed 0-4 item whose pass is 3 or
#'   more. Items whose steps differ are not the same kind of milestone; the `step` column says which.
#' @return A `kits_item_table`, a data.frame with one row per item; `attr(, "norm")` keeps the norm,
#'   `attr(, "rp")` the RP; with `data`, the diagnostics list the items with perfect separation by
#'   age (their age is only known between the oldest child who fails and the youngest who passes).
#'   Its columns:
#'   * `D`: the threshold of the pass step (`item_score == step`), on the scale of the norm.
#'   * `step`: the pass step, only when `step` is not 1 for every item.
#'   * `domain`, `column`: the domain of the item and its column of the ladder (default the domain).
#'   * `age50` (`age75`, ... for another `rp`): the age at which a share `rp` pass, from `data`.
#'   * `DIF`: `"none"`, `"favors girls"`, `"favors boys"` or `"inconsistent"`, from `dif`.
#'   * `stem`: the item stem (any language), from `stems`; the item id when not given.
#'   * `age_reach`, `below_range`, `above_range`, `empirical`, `x`: the place on the growth chart
#'   (the empirical age when there is one, else the age at which the norm median reaches `D`).
#' @examples
#' set.seed(1)
#' age <- runif(800, 3, 24)
#' D <- 20 + 15 * log(age) + rnorm(800, 0, 3)
#' norm <- kits_norms(D, age, ages = 3:24, degrees = c(3, 2, 1, 1),
#'                    by = data.frame(domain = rep("lan", 800)))
#' thr <- data.frame(item_id = c("a1", "a2", "a3", "a4"), item_score = 1, thr = c(40, 55, 65, 60),
#'                   domain = "lan")
#' pass <- sapply(c(40, 55, 65, 60), \(b) rbinom(800, 1, plogis((D - b) / 3)))
#' colnames(pass) <- thr$item_id
#' resp <- data.frame(visit_id = 1:800, age_m = age, wave_m = round(age), pass)
#' dif <- data.frame(item_id = c("a2", "a3", "a4"),
#'                   DIF = c("favors girls", "favors boys", "inconsistent"))
#' stems <- c(a1 = "Lifts head", a2 = "Says two words", a3 = "Names colours",
#'            a4 = "Points to pictures")
#' item_tab <- kits_item_table(thr, norm, data = resp, dif = dif, stems = stems)
#' item_tab
#'
#' # A 0-4 item whose pass is a score of 3 or more
#' thr4 <- rbind(thr, data.frame(item_id = "a5", item_score = 1:4, thr = c(30, 45, 58, 70),
#'                               domain = "lan"))
#' kits_item_table(thr4, norm, step = c(a5 = 3))
#' p <- plot(item_tab)
#' p <- plot(item_tab, type = "ladder", ages = c(6, 12, 18))
#' if (interactive()) p
#' @export
kits_item_table <- function(thr, norms, data = NULL, dif = NULL, stems = NULL, rp = 0.5,
                            dif_cut = 0.3, step = 1) {
  need_rp(rp)
  need_pass_step(step)
  norm <- table_norm(norms)
  steps <- item_steps(thr, norm, rp, step)
  ages <- table_ages(data, rp, step)
  items <- steps |>
    left_join(ages, by = "item_id") |>
    left_join(table_dif(dif, dif_cut), by = "item_id") |>
    left_join(table_stems(stems), by = "item_id") |>
    mutate(DIF = coalesce(DIF, "none"), stem = coalesce(stem, item_id))
  attr(items, "rp") <- if (rp != 0.5) rp
  out <- item_positions(items, norm, rp)
  out <- kits_table(as.data.frame(out), "kits_item_table")
  if (rp != 0.5) attr(out, "rp") <- rp
  attr(out, "dif_labels") <- table_dif_labels(dif)   # 兩個方向的標籤（plot() 的顏色照方向）
  attr(out, "norm") <- norm
  attr(out, "scale") <- scale_of(norms)
  attr(out, "diagnostics") <- attr(ages, "diagnostics")
  out
}

# 常模：一個領域一個（by = "domain"，或一個常模）；有其他的組時停下來
table_norm <- function(norms) {
  abort_if(!is.data.frame(norms), "`norms` must be a kits_norms() result.")
  others <- setdiff(attr(norms, "by"), "domain")
  abort_if(length(others) > 0, "The item table uses one norm per domain; this norm is also by ",
           paste(others, collapse = ", "), ". Fit it by domain only, e.g. kits_norms(scores, by = \"domain\").")
  norm <- as.data.frame(norms)
  need_cols(norm, c("age_m", "P10", "P50", "P90"), "norms")
  if (!"domain" %in% names(norm)) norm$domain <- NA_character_
  norm
}

# 每題通過的 step 的門檻（D）、領域、階梯的欄；同一題給了兩次時留第一列
item_steps <- function(thr, norm, rp = 0.5, step = 1) {
  one <- is.data.frame(thr) || inherits(thr, "kits_item") || !is.list(thr)
  parts <- if (one) list(thr) else thr
  steps <- parts |>
    map(\(x) as.data.frame(rp_steps(x, rp))) |>
    list_rbind()
  need_cols(steps, c("item_id", "item_score", "thr"), "thr")
  ids <- unique(as.character(steps$item_id))
  extra <- setdiff(names(step), ids)
  abort_if(length(extra) > 0, "`step` names item(s) that are not in `thr`: ", some_ids(extra, 10), ".")
  pass_step <- step_of(ids, step)
  no_step <- ids[!map_lgl(seq_along(ids), \(j) any(steps$item_id == ids[j] & steps$item_score == pass_step[j]))]
  abort_if(length(no_step) > 0, "Item(s) without their pass step in `thr`: ", some_ids(no_step, 10), ".")
  steps <- filter(steps, item_score == pass_step[match(as.character(item_id), ids)])
  twice <- unique(steps$item_id[duplicated(steps$item_id)])
  if (length(twice) > 0) {
    kits_inform("Item(s) given twice keep their first row: ",
                paste(utils::head(twice, 10), collapse = ", "), ".",
                class = "kits_message_data")
    steps <- filter(steps, !duplicated(item_id))
  }
  domains <- unique(norm$domain)
  if (!"domain" %in% names(steps)) {
    abort_if(length(domains) > 1,
             "The norm has several domains (", paste(domains, collapse = ", "), "): give `thr` a `domain` column.")
    steps$domain <- domains
  }
  if (!"column" %in% names(steps)) steps$column <- steps$domain
  out <- transmute(steps, item_id = as.character(item_id), domain, column, D = thr)
  if (any(steps$item_score != 1)) out$step <- steps$item_score
  out
}

# step：一個整數（每題都用），或用題號當名字的向量（沒列的題用 1）
need_pass_step <- function(step) {
  ok <- is.numeric(step) && length(step) >= 1 && all(is.finite(step)) && all(step >= 1) &&
    all(step == round(step)) && (length(step) == 1 || !is.null(names(step)))
  abort_if(!ok, "`step` must be one whole number >= 1, or whole numbers >= 1 named by item, ",
           "e.g. c(cogD01 = 3).")
}

# 每題通過的 step
step_of <- function(items, step) {
  if (is.null(names(step))) return(rep(step, length(items)))
  out <- rep(1, length(items))
  hit <- items %in% names(step)
  out[hit] <- step[items[hit]]
  out
}

# 一張門檻表在 RP rp：kits_item() 用 coef(x, rp = rp)；表要記著同一個 RP（沒記的是 RP50），不然停下來
rp_steps <- function(x, rp) {
  if (inherits(x, "kits_item")) {
    return(coef(x, rp = rp))
  }
  abort_if(is.data.frame(x) && rp_of(x) != rp,
           sprintf("`thr` holds RP%g thresholds, but `rp` = %g. ", 100 * rp_of(x), rp),
           sprintf("Give the RP%g thresholds (coef(item, rp = %g)) or set `rp` = %g.",
                   100 * rp, rp, rp_of(x)))
  x
}

# 通過月齡：data 的作答（分數 >= step 是通過）算出 rp 的孩子通過的月齡，欄名 age50、age80……；
# 完全分離的題記在 attr "diagnostics"
table_ages <- function(data, rp, step = 1) {
  out <- data.frame(item_id = character(), age = numeric())
  if (!is.null(data)) {
    data <- kits_resp(data)
    need_cols(data, c("age_m", "wave_m"), "data")
    X <- wide_matrix(data, "data")
    ap <- age_pass(X, data$age_m, data$wave_m, step = step_of(colnames(X), step), p = rp)
    out <- structure(data.frame(item_id = ap$item_id, age = ap[[2]]),
                     diagnostics = attr(ap, "diagnostics"))
  }
  names(out)[2] <- age_col(rp)
  out
}

# DIF：kits_dif_combine() 的 DIF 欄，或一個 kits_dif() 照 cut（logit）判斷方向，標籤用它的組
table_dif <- function(dif, cut = 0.3) {
  if (is.null(dif)) {
    return(data.frame(item_id = character(), DIF = character()))
  }
  if ("DIF" %in% names(dif)) {
    return(data.frame(item_id = as.character(dif$item_id), DIF = as.character(dif$DIF)))
  }
  need_cols(dif, c("item_id", "dif"), "dif")
  need_dif_cut(cut, "dif_cut")
  g <- dif_sides(dif)
  side <- coalesce(dif_direction(dif$dif, cut, g$focal, g$ref), "none")
  data.frame(item_id = as.character(dif$item_id), DIF = side)
}

# DIF 兩個方向的標籤（up：對 focal 有利；down：對 ref 有利）；沒有 DIF 時 NULL
table_dif_labels <- function(dif) {
  if (is.null(dif)) return(NULL)
  g <- dif_sides(dif)
  dif_labels(g$focal, g$ref)
}

# 題幹：data.frame（item_id、stem）或用題號當名字的文字向量
table_stems <- function(stems) {
  if (is.null(stems)) {
    return(data.frame(item_id = character(), stem = character()))
  }
  if (is.character(stems) && !is.null(names(stems))) {
    return(data.frame(item_id = names(stems), stem = unname(stems)))
  }
  need_cols(stems, c("item_id", "stem"), "stems")
  data.frame(item_id = as.character(stems$item_id), stem = as.character(stems$stem))
}

#' @rdname kits_item_table
#' @param x A `kits_item_table`.
#' @param type `"growth"` (the default; growth curves with the items) or `"ladder"` (the
#'   developmental ladder).
#' @param items Which items to draw (item ids); by default all of them.
#' @param ages For `type = "ladder"`: the ages of the ticks on the right.
#' @param age_ticks For `type = "growth"`: the age ticks.
#' @param y_range For `type = "growth"`: the range of the D axis, the same in every panel; by
#'   default from the data.
#' @param ... Not used.
#' @export
plot.kits_item_table <- function(x, type = c("growth", "ladder"), items = NULL,
                                 ages = c(3, 6, 9, 12, 18, 24, 36, 48, 60, 72, 84),
                                 age_ticks = AGE_TICKS,
                                 y_range = NULL, ...) {
  type <- match_option(type, c("growth", "ladder"), "type")
  norm <- attr(x, "norm")
  abort_if(is.null(norm), "This table keeps no norm; make it with kits_item_table().")
  shown <- x[x$item_id %in% pick_items(x$item_id, items), , drop = FALSE]
  attr(shown, "rp") <- attr(x, "rp")
  domains <- unique(shown$domain)
  if (type == "growth") {
    return(growth_plot(norm, shown, domains, age_ticks, y_range))
  }
  marks <- ladder_marks(norm, domains, ages)
  ladder_plot(shown, unique(shown$column), marks)
}
