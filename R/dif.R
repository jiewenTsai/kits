# 組間的差：Andersen (1973) 的概似比檢定，兩組時再加每題的 DIF（Bechger & Maris, 2015）。
#
# 把人次分組，每組各自用 CML 校準，和全部一起校準的條件概似比：
#   LR = 2 (sum_g log L_g - log L_all)，df = sum_g (組 g 自由的步驟參數) - (全部自由的步驟參數)
# 條件概似只和步驟參數有關（原始分是能力的充分統計量），所以各組的能力分布不同也沒關係。
#
# 一組裡觀察不到的步驟（例如年幼的波次沒有人拿到最高分）：這組只估它觀察得到的步驟。
# 這組的條件概似的最大值在參數空間的邊上（那個 beta 是正或負無限大），極限正好是把這題截短
# （或往下平移）之後的條件概似，所以 log L_g 用截短的校準算，這些步驟不算這組的自由參數。
# 作答不變，所以 log L_all 照常是全部作答的條件概似。
# 分數中間有缺（例如這組只有 0、1、3 分）時不是 PCM 了：這組不用這題，全部一起的校準也不用這些作答。
# unobserved = "drop"：和 eRm 的 LRtest() 相同，一題在任何一組缺分數類別，就每一組都不用、全部也不用。
#
# 兩組時每題的 DIF 和 dexter 的 DIF() 相同：每組用它所有的作答各自校準（不只有資訊的人次），
# 兩組都有、分數也相同的題才比較。

#' Differences between groups: LR test and item DIF
#'
#' Calibrates each group of person-waves on its own by CML and tests whether the step parameters
#' are the same (Andersen, 1973): LR = 2 (sum of log L_g - log L_all), on the steps each group
#' observes. It returns each group's steps on a common origin and the Wald z of the difference
#' between every two groups; with large N, read these sizes rather than the p value.
#' With two groups it also gives each item's DIF (the item-pair method of Bechger & Maris, 2015, as
#' dexter's `DIF()`): `dif` is the mean over the item's steps, in logits, positive is harder for
#' the `ref` group (with `ref = "boy"`, favors girls). `print(d)` shows the LR test.
#'
#' Details: `kits_algorithm("dif")` and tutorial_algorithms/05_checks.html.
#'
#' @param data Responses: a matrix (one row per person-wave, one column per item, row names are
#'   person ids), a wide or long data.frame, or a [kits_resp()] (see [kits_item()]).
#' @param group The groups:
#'   * a vector or factor with one group per person-wave: one value per row of `data` (for a long
#'     data.frame, one per response row), or named by visit_id; NA is left out, with a message;
#'   * `"raw"`: low and high raw scores, split at the median within each item pattern (the median
#'     goes to the low group);
#'   * `"wave"`: one group per wave (`wave_m` of a [kits_resp()], or of `persons`);
#'   * the name of another column of a [kits_resp()] or of `persons`, for example `"sex"`.
#'
#'   Not given: the column `by` of a [kits_resp()].
#' @param item Optional: the [kits_item()] calibration of the same `data`, used as the joint model
#'   (a calibration of other data stops with an error). By default it is calibrated here.
#' @param ref The reference group of the item DIF (two groups). By default `"boy"` when the groups
#'   are boy and girl, otherwise the first group.
#' @param by With a [kits_resp()] and no `group`: the column with the groups, by default `"sex"`.
#' @param persons Optional data.frame with `visit_id` and the column named by `group` (for example
#'   `wave_m`), when `data` is not a [kits_resp()] with that column.
#' @param unobserved Steps that a group does not observe: `"steps"` (default) calibrates each group
#'   on the steps it observes; `"drop"` leaves an item that lacks a score category in any group out
#'   of every group and of the joint calibration (as eRm's `LRtest()`).
#' @return A data.frame (class `kits_dif`): with two groups one row per compared item, `item_id`
#'   and `dif` (logits; no rows with more groups), and the attributes
#'   * `lrtest`: `LR`, `df`, `p`, `split` (how the groups were formed), `rule` (`unobserved`) and
#'     `groups`, one row per group: `group`, `n`, `n_info` (informative person-waves), `n_items`,
#'     `n_steps` (steps estimated), `n_par` (free step parameters), `logLik`;
#'   * `steps`: one row per group and estimated step: `group`, `item_id`, `item_score`, `beta`, `se`
#'     (on the origin of the joint calibration);
#'   * `full`: the joint calibration: `item_id`, `item_score`, `beta`, `se`, with attributes
#'     `logLik` and `n_par`;
#'   * `diff`: one row per pair of groups and step: `item_id`, `item_score`, `group_1`, `group_2`,
#'     `diff` (beta in `group_2` minus beta in `group_1`), `se`, `z`;
#'   * `unobserved`: the steps of answered items that a group does not estimate: `group`, `item_id`,
#'     `item_score`, `reason`, `left_out` (the item's responses in that group are left out of the
#'     joint calibration too);
#'   * with two groups (otherwise NULL): `overall` (`stat`, `df`, `p`: the test that any item has
#'     DIF), `groups` (`ref`, `focal` and the person-waves `n` of each), `se` (of each `dif`) and
#'     `left_out` (items not compared: in one group only, or with different scores);
#'   * `diagnostics` (see [kits_check()]) and `info` (see [kits_info()]).
#' @seealso [kits_dif_combine()] to combine two cohorts, [kits_drift()], [kits_itemfit()].
#' @references
#' Andersen, E. B. (1973). A goodness of fit test for the Rasch model. *Psychometrika, 38*, 123-140.
#' \doi{10.1007/BF02291180}
#'
#' Bechger, T. M., & Maris, G. (2015). A statistical test for differential item pair functioning.
#' *Psychometrika, 80*(2), 317-340. \doi{10.1007/s11336-014-9408-y}
#' @examples
#' betas <- list(a1 = c(-1, 0), a2 = c(0, 1), a3 = 0.5, a4 = -0.5, a5 = 0)
#' theta <- seq(-3, 3, length.out = 300)
#' boys <- kits_sim(theta, betas, prefix = "b", seed = 1)
#' girls <- kits_sim(theta, modifyList(betas, list(a5 = 1)), prefix = "g", seed = 2)
#' X <- rbind(boys, girls)
#'
#' # a5 is harder for girls: negative, favors boys
#' dif <- X |> kits_dif(rep(c("boy", "girl"), each = 300))
#' dif
#' attr(dif, "lrtest")$p
#'
#' # Andersen's median split; do the steps stay the same across waves?
#' kits_dif(X, "raw")
#' resp <- kitm03$responses |>
#'   dplyr::filter(domain == "lan", wave_m <= 12) |>
#'   kits_resp()
#' kits_dif(resp, "wave")
#' @export
kits_dif <- function(data, group, item = NULL, ref = NULL, by = NULL, persons = NULL,
                     unobserved = c("steps", "drop"), ...) {
  no_extra_args(..., fun = "kits_dif")
  call <- call_text(sys.call())
  data_expr <- substitute(data)
  unobserved <- match_option(unobserved, c("steps", "drop"), "unobserved")
  if (missing(group)) group <- dif_by(data, by, persons)
  all_resp <- resp_long(data)
  item <- lrt_check_item(item, all_resp)

  # 1 分組（沒有組的人次不用）；每一題都要有 0 分、至少兩種分數（和 kits_item() 相同）
  sp <- lrt_split(group, data, all_resp, persons)
  resp <- as.data.frame(all_resp)[c("visit_id", "item_id", "item_score")]
  resp$group <- sp$group
  n_na <- lrt_n_na(resp)
  resp <- resp[!is.na(resp$group), ]
  left_out_items(resp, NULL)
  levels <- sp$levels[sp$levels %in% resp$group]
  lrt_need_groups(levels)
  ref <- dif_ref(ref, levels)

  # 2 每組估得出來的步驟；分數有缺的題（和 "drop" 時缺類別的題）這組不用，全部一起的校準也不用
  use <- lrt_usable(resp, levels, unobserved)
  lrt_need_groups(names(use$groups))
  lrt_tell_unobserved(use$unobserved)

  # 3 每組的 CML，和全部一起的 CML（同一份作答）
  fits <- imap(use$groups, \(d, g) lrt_group_fit(d, n = n_distinct(resp$visit_id[resp$group == g])))
  joint <- use$resp[c("visit_id", "item_id", "item_score")]
  same <- !is.null(item) && nrow(joint) == nrow(all_resp)
  full <- if (same) item$fit else lrt_joint_fit(joint, nrow(joint) < nrow(resp))
  res <- lrt_result(fits, full, names(use$groups))

  # 4 兩組時每題的 DIF
  out <- dif_items(resp, levels, ref)
  attr(out, "lrtest") <- c(res[c("LR", "df", "p")], list(split = sp$label, rule = unobserved,
                                                         groups = res$groups))
  for (a in c("steps", "full", "diff")) attr(out, a) <- res[[a]]
  attr(out, "unobserved") <- kits_table(use$unobserved)
  attr(out, "diagnostics") <- dif_diagnostics(out, n_na)
  attr(out, "info") <- make_info(call, list(group = sp$label, ref = ref, unobserved = unobserved),
                                 data_fingerprint(all_resp),
                                 persons = persons_fingerprint(distinct(resp, visit_id, group)),
                                 input = data_input(data, data_expr))
  out
}

# 沒有給 group：kits_resp 的 by 欄（預設 sex）
dif_by <- function(data, by, persons) {
  abort_if(!inherits(data, "kits_resp") && is.null(persons),
           "Give `group` (one value per row, \"raw\", \"wave\" or a column name), or the ",
           "responses as a kits_resp(): kits_dif(resp, by = \"sex\").")
  if (is.null(by)) "sex" else by
}

# item：NULL（在這裡校準），或同一份資料的 kits_item（kits() 的結果也可以）。
# 回傳可以直接用的 item（沒有固定參數的自由校準），否則 NULL
lrt_check_item <- function(item, all_resp) {
  if (is.null(item)) {
    return(NULL)
  }
  if (inherits(item, "kits")) item <- item$item
  abort_if(!inherits(item, "kits_item"),
           "`item` must be the kits_item() calibration of `data`, or NULL (the joint model is ",
           "then calibrated here).")
  abort_if(!identical(item$info$data, data_fingerprint(all_resp)),
           "`item` was calibrated on other data than `data`; the LR test needs the joint ",
           "calibration of the same responses. Give kits_item(data), or leave `item` out.")
  if (!identical(item$info$settings$fixed_params, "none") || length(item$left_out) > 0) {
    kits_inform("`item` has fixed parameters, so it is not the joint CML of these data; ",
                "the joint model is calibrated here without them.", class = "kits_message_fit")
    return(NULL)
  }
  item
}

# ---- 1 分組 ------------------------------------------------------------------------------------
# 回傳 group（長格式每一列的組，文字）、levels（組的順序）、label（怎麼分的）
lrt_split <- function(group, data, resp, persons) {
  abort_if(is.null(group) || length(group) == 0 || is.list(group) || !is.atomic(group),
           "`group` must be \"raw\", \"wave\", a column name, or one group per person-wave.")
  if (is.character(group) && length(group) == 1 && is.null(names(group))) {
    abort_if(is.na(group), "`group` is NA; give \"raw\", \"wave\" or a column name.")
    if (group == "raw") {
      return(lrt_raw_split(resp))
    }
    col <- if (group == "wave") "wave_m" else group
    value <- lrt_person_column(col, data, resp, persons)
    label <- col
  } else {
    value <- per_person(group, data, resp, "group")
    label <- "the vector `group`"
  }
  lv <- if (is.factor(value)) levels(droplevels(value)) else as.character(sort(unique(value)))
  list(group = as.character(value), levels = lv, label = label)
}

# 依原始分：每個作答組合（答的題）裡，原始分 <= 中位數的是 low，其他是 high（同分歸 low，和 eRm 相同）。
# 不同作答組合的原始分不能比（短的題本 5 分是高分，長的題本是低分），所以中位數在作答組合裡算。
lrt_raw_split <- function(resp) {
  persons <- resp |>
    person_sets() |>
    mutate(group = if_else(raw <= stats::median(raw), "low", "high"), .by = set)
  list(group = persons$group[match(resp$visit_id, persons$visit_id)], levels = c("low", "high"),
       label = "raw score (median within each item pattern; ties low)")
}

# 孩子的一欄（wave_m、sex、...）：從 persons，或作答表（kits_resp）的欄，照 visit_id 對上長格式每一列
lrt_person_column <- function(col, data, resp, persons) {
  if (!is.null(persons)) {
    src <- old_ids(persons, "persons", quiet = TRUE)
    need_cols(src, c("visit_id", col), "persons")
  } else if (inherits(data, "kits_resp")) {
    abort_if(!col %in% resp_person_cols(data),
             sprintf("The responses have no column `%s` (columns of the child: %s).", col,
                     paste(resp_person_cols(data), collapse = ", ")))
    src <- as.data.frame(data)
  } else {
    kits_abort(sprintf("group = \"%s\" needs the column %s: give the responses as a kits_resp(), ",
                       col, col),
               sprintf("or `persons` (a data.frame with visit_id and %s).", col),
               class = "kits_error_input")
  }
  ids <- as.character(src$visit_id)
  abort_if(anyDuplicated(ids),
           "`persons` has repeated visit_id; give one row per person-wave.", class = "kits_error_data")
  miss <- setdiff(unique(resp$visit_id), ids)
  abort_if(length(miss) > 0,
           sprintf("%d person-wave(s) of the data have no row in `persons`%s.", length(miss), ids_text(miss)))
  src[[col]][match(resp$visit_id, ids)]
}

# 沒有組（NA）的人次數；有的話說一聲（這些人次不用）
lrt_n_na <- function(resp) {
  n <- n_distinct(resp$visit_id[is.na(resp$group)])
  if (n > 0) {
    kits_inform(sprintf("%d person-wave(s) without a group (NA in `group`) are left out.", n),
                class = "kits_message_data")
  }
  n
}

lrt_need_groups <- function(levels) {
  abort_if(length(levels) < 2, "The LR test needs at least two groups with items that can be calibrated; there is ",
           length(levels), if (length(levels) == 1) sprintf(" (%s)", levels) else "", ".", class = "kits_error_data")
}

# DIF 的參照組：沒給時，boy、girl 兩組是 boy，否則第一組
dif_ref <- function(ref, levels) {
  if (is.null(ref)) {
    return(if (setequal(levels, c("boy", "girl"))) "boy" else levels[1])
  }
  abort_if(length(ref) != 1 || is.na(ref) || !as.character(ref) %in% levels,
           sprintf("`ref` = %s is not a group. The groups are %s; set the reference group ",
                   deparse(ref)[1], paste(levels, collapse = ", ")),
           sprintf("with `ref`, for example ref = \"%s\".", levels[1]))
  as.character(ref)
}

# ---- 2 每組估得出來的步驟 ------------------------------------------------------------------------
# 回傳 resp（全部一起校準的作答）、groups（每組校準的作答：分數減去 lo，和每題的 lo）、
# unobserved（每組答了、但沒有估的步驟和原因）
lrt_usable <- function(resp, levels, unobserved) {
  top <- tapply(resp$item_score, resp$item_id, max)
  removed <- data.frame(group = character(), item_id = character(), reason = character())
  if (unobserved == "drop") {   # eRm：缺分數類別的題，每一組都不用、全部也不用
    removed <- lrt_erm_drop(resp, levels, top)
    resp <- resp[!resp$item_id %in% removed$item_id, ]
  }
  answered <- distinct(resp, group, item_id)
  groups <- list()
  for (g in levels) {
    found <- lrt_group_steps(resp[resp$group == g, ], top)
    if (nrow(found$gaps) > 0) {   # 分數有缺：這組不用這題，全部一起的校準也不用這些作答
      resp <- resp[!(resp$group == g & resp$item_id %in% found$gaps$item_id), ]
      removed <- bind_rows(removed, found$gaps)
    }
    if (nrow(found$data) > 0) groups[[g]] <- found
  }
  list(resp = resp, groups = groups, unobserved = lrt_unobserved_table(groups, answered, removed, top,
                                                                        levels))
}

# eRm 的規則（unobserved = "drop"）：一題在答它的任何一組缺了分數類別，每一組都不用
lrt_erm_drop <- function(resp, levels, top) {
  cats <- resp |>
    summarise(k = n_distinct(item_score), scores = paste(sort(unique(item_score)), collapse = ", "),
              .by = c(group, item_id))
  bad <- cats[cats$k < top[cats$item_id] + 1, ]   # 沒有時回傳 0 列的表
  first <- bad[!duplicated(bad$item_id), ]
  out <- cats[cats$item_id %in% bad$item_id, c("group", "item_id")]
  j <- match(out$item_id, first$item_id)
  out$reason <- sprintf(paste0("unobserved = \"drop\": left out of every group (group %s has scores ",
                               "%s only)"), first$group[j], first$scores[j])
  out
}

# 一組估得出來的步驟（重複到不再改變）：
# * 只看有資訊的人次：答了至少兩題（這組校準的題），原始分不是 0 也不是滿分；其他人次的條件概似是 1；
# * 每題的分數範圍 lo..hi（有資訊的人次裡）：估 lo + 1 .. hi 步，分數減去 lo；只有一種分數的題不估；
# * 分數有缺的題（例如 0、1、3）不是 PCM：記在 gaps，這組不用它（重新開始）。
# 回傳 data（visit_id、item_id、item_score = 分數 - lo）、lo、hi、gaps、n_info
lrt_group_steps <- function(rg, top) {
  gaps <- data.frame(group = character(), item_id = character(), reason = character())
  repeat {
    d <- rg[!rg$item_id %in% gaps$item_id, c("visit_id", "item_id", "item_score")]
    lo <- stats::setNames(numeric(length(top)), names(top))
    hi <- top
    repeat {
      d <- lrt_informative(d, lo, hi)
      rng <- d |>
        summarise(lo = min(item_score), hi = max(item_score), gap = has_gap(item_score),
                  scores = paste(sort(unique(item_score)), collapse = ", "), .by = item_id)
      gap <- rng[rng$gap, ]
      if (nrow(gap) > 0) break
      new_lo <- stats::setNames(rng$lo, rng$item_id)
      new_hi <- stats::setNames(rng$hi, rng$item_id)
      keep <- rng$item_id[rng$hi > rng$lo]
      changed <- any(new_lo != lo[rng$item_id]) || any(new_hi != hi[rng$item_id]) ||
        !all(unique(d$item_id) %in% keep)
      lo[rng$item_id] <- new_lo
      hi[rng$item_id] <- new_hi
      d <- d[d$item_id %in% keep, ]
      if (!changed) break
    }
    if (nrow(gap) == 0) break
    gaps <- bind_rows(gaps, data.frame(
      group = rg$group[1], item_id = gap$item_id,
      reason = sprintf("scores %s in this group (a gap): left out of this group", gap$scores)))
  }
  items <- unique(d$item_id)
  calib <- d
  calib$item_score <- calib$item_score - lo[calib$item_id]
  list(data = calib, lo = lo[items], hi = hi[items], gaps = gaps, n_info = n_distinct(d$visit_id))
}

# 有資訊的人次：答了至少兩題，0 < 原始分 < 滿分（滿分照目前每題的 lo..hi 算）
lrt_informative <- function(d, lo, hi) {
  p <- d |>
    summarise(n = n(), raw = sum(item_score - lo[item_id]), top = sum(hi[item_id] - lo[item_id]),
              .by = visit_id)
  ok <- p$visit_id[p$n >= 2 & p$raw > 0 & p$raw < p$top]
  d[d$visit_id %in% ok, ]
}

# 每組答了、但沒有估的步驟（照原本的步驟編號）和原因；left_out：這些作答全部一起的校準也不用
lrt_unobserved_table <- function(groups, answered, removed, top, levels) {
  steps_of <- function(its) data.frame(item_id = rep(its, top[its]), item_score = sequence(top[its]))
  rows <- map(levels, \(g) {
    its <- answered$item_id[answered$group == g]
    if (length(its) == 0) {
      return(NULL)
    }
    s <- steps_of(its)
    s$group <- g
    f <- groups[[g]]
    lo <- if (is.null(f)) numeric() else f$lo
    hi <- if (is.null(f)) numeric() else f$hi
    estimated <- s$item_id %in% names(lo) & s$item_score > lo[s$item_id] &
      s$item_score <= hi[s$item_id]
    out <- s[!estimated, ]
    gone <- removed[removed$group == g, ]
    out$left_out <- out$item_id %in% gone$item_id
    out$reason <- case_when(
      out$left_out ~ gone$reason[match(out$item_id, gone$item_id)],
      !out$item_id %in% names(lo) ~ "no information in this group (no informative person-wave answered it)",
      out$item_score <= lo[out$item_id] ~ "reached by every informative person-wave of this group",
      .default = "reached by no informative person-wave of this group")
    out
  })
  out <- list_rbind(c(list(data.frame(group = character(), item_id = character(),
                                      item_score = integer(), reason = character(),
                                      left_out = logical())), rows))
  # eRm 規則拿掉的題：這組不再有它（answered 是拿掉之後的），補上它的步驟
  extra <- removed[!paste(removed$group, removed$item_id) %in% paste(out$group, out$item_id), ]
  if (nrow(extra) > 0) {
    k <- rep(seq_len(nrow(extra)), top[extra$item_id])   # 每一題的每一步一列
    out <- bind_rows(out, data.frame(group = extra$group[k], steps_of(extra$item_id),
                                     reason = extra$reason[k], left_out = TRUE))
  }
  out <- out[order(match(out$group, levels), sort_rank(out$item_id), out$item_score),
             c("group", "item_id", "item_score", "reason", "left_out")]
  rownames(out) <- NULL
  out
}

# 沒有估的步驟（"item:score (group)"）和不用的題（"item (group)"），給訊息、診斷、print
lrt_unobserved_ids <- function(tab) {
  steps <- tab[!tab$left_out, ]
  gone <- distinct(tab[tab$left_out, ], group, item_id)
  list(steps = sprintf("%s:%s (%s)", steps$item_id, steps$item_score, steps$group),
       items = sprintf("%s (%s)", gone$item_id, gone$group))
}

# 訊息（沒有估的步驟）和警告（不用的作答）
lrt_tell_unobserved <- function(tab) {
  ids <- lrt_unobserved_ids(tab)
  if (length(ids$steps) > 0) {
    kits_inform(sprintf("%d step(s) are not observed in a group, so that group does not estimate ",
                        length(ids$steps)),
                "them (each group uses the steps it observes): ", some_ids(ids$steps, 6),
                ". See attr(x, \"unobserved\").", class = "kits_message_data")
  }
  if (length(ids$items) > 0) {
    kits_warn("Item(s) left out of a group and of the joint calibration: ", some_ids(ids$items, 10),
              ". See attr(x, \"unobserved\") for the reasons.", class = "kits_warning_data")
  }
}

# ---- 3 每組的 CML -------------------------------------------------------------------------------
# 全部一起的 CML；拿掉了作答（分數有缺、"drop"）以後估不出來時，說是拿掉作答的關係
lrt_joint_fit <- function(joint, reduced) {
  if (!reduced) {
    return(cml_fit(joint))
  }
  tryCatch(cml_fit(joint), kits_error = function(e) {
    kits_abort("After leaving out the items named in the warning, the joint calibration fails: ",
               conditionMessage(e), " With unobserved = \"drop\", try unobserved = \"steps\".",
               class = "kits_error_data")
  })
}

# 一組的題分成幾片沒有共同題的組時，每一片各自校準（各自的原點），自由的參數是每片的步驟數減 1。
# 步驟照原本的編號（校準的分數減去了 lo，所以步驟 s 是原本的 lo + s）
lrt_group_fit <- function(found, n) {
  d <- found$data
  items <- sort(unique(d$item_id))
  sets <- unique(person_sets(d)$set)
  part <- item_groups(map(sets, \(s) match(set_items(s), items)), length(items))
  fits <- map(split(items, part), \(its) cml_fit(d[d$item_id %in% its, ]))
  coef <- imap(unname(fits), \(f, k) {
    data.frame(f$coef[c("item_id", "item_score", "beta", "SE_beta")], part = k)
  }) |>
    list_rbind()
  coef$item_score <- as.integer(coef$item_score + found$lo[coef$item_id])
  list(coef = coef, nll = sum(map_dbl(fits, \(f) f$est$nll)), n_par = nrow(coef) - length(fits),
       n = n, n_info = found$n_info)
}

# 結果：檢定、組、每組的步驟參數（共同的原點）、全部一起的校準、組間的差和 Wald z
lrt_result <- function(fits, full, levels) {
  fc <- full$coef
  full_tab <- data.frame(item_id = fc$item_id, item_score = fc$item_score, beta = fc$beta,
                         se = fc$SE_beta)
  n_par_all <- nrow(fc) - 1L
  ll_g <- map_dbl(fits, \(f) -f$nll)
  LR <- 2 * (sum(ll_g) + full$est$nll)
  df <- as.integer(sum(map_dbl(fits, "n_par")) - n_par_all)
  abort_if(df < 1, "The groups have no more free step parameters than the joint calibration, so there ",
           "is nothing to test.", class = "kits_error_data")
  groups <- data.frame(group = levels, n = map_int(fits, "n"), n_info = map_int(fits, "n_info"),
                       n_items = map_int(fits, \(f) n_distinct(f$coef$item_id)),
                       n_steps = map_int(fits, \(f) nrow(f$coef)),
                       n_par = as.integer(map_dbl(fits, "n_par")), logLik = unname(ll_g))
  steps <- imap(fits, \(f, g) lrt_align(f$coef, full_tab, g)) |>
    unname() |>
    list_rbind()
  list(LR = LR, df = df, p = stats::pchisq(LR, df, lower.tail = FALSE),
       groups = kits_table(groups), steps = kits_table(steps),
       full = structure(kits_table(full_tab), logLik = -full$est$nll, n_par = n_par_all),
       diff = kits_table(lrt_diff(steps, levels)))
}

# 一組的步驟參數放到共同的原點：每一片平移，使它的 beta 的平均等於全部一起的 beta（同樣的步驟）的平均
lrt_align <- function(cf, full_tab, g) {
  key <- paste(cf$item_id, cf$item_score)
  ref <- full_tab$beta[match(key, paste(full_tab$item_id, full_tab$item_score))]
  shift <- stats::ave(ref, cf$part) - stats::ave(cf$beta, cf$part)
  data.frame(group = g, item_id = cf$item_id, item_score = cf$item_score, beta = cf$beta + shift,
             se = cf$SE_beta)
}

# 每兩組（照組的順序 g < h）都有的步驟：差 = beta_h - beta_g，SE = sqrt(SE_g^2 + SE_h^2)，Wald z
lrt_diff <- function(steps, levels) {
  pairs <- utils::combn(levels, 2, simplify = FALSE)
  map(pairs, \(pr) {
    a <- steps[steps$group == pr[1], ]
    b <- steps[steps$group == pr[2], ]
    d <- inner_join(a, b, by = c("item_id", "item_score"), suffix = c("_1", "_2"))
    data.frame(item_id = d$item_id, item_score = d$item_score, group_1 = d$group_1,
               group_2 = d$group_2, diff = d$beta_2 - d$beta_1, se = sqrt(d$se_1^2 + d$se_2^2)) |>
      mutate(z = diff / se)
  }) |>
    list_rbind() |>
    arrange(match(group_1, levels), match(group_2, levels), sort_rank(item_id), item_score)
}

# ---- 4 兩組時每題的 DIF ----------------------------------------------------------------------------
# 每組各自用它所有的作答校準（和 dexter 的 DIF() 相同）；這組不能校準的題（沒有 0 分、只有一種分數、
# 分數有缺）這組不用。兩組都有、分數也相同的題才比較。多於兩組（或可以比的題少於 2 題、一組不能
# 校準）時沒有列
dif_items <- function(resp, levels, ref) {
  out <- kits_table(data.frame(item_id = character(), dif = numeric()), "kits_dif")
  if (length(levels) != 2) {
    return(out)
  }
  by_group <- split(resp[c("visit_id", "item_id", "item_score")], resp$group)
  groups <- names(by_group)
  focal <- setdiff(groups, ref)
  n <- map_int(by_group, \(g) n_distinct(g$visit_id))
  items <- resp |>
    summarise(ok = can_calibrate(item_score) && !has_gap(item_score),
              scores = paste(sort(unique(item_score)), collapse = ", "), .by = c(group, item_id))
  by_group <- imap(by_group, \(g, k) g[g$item_id %in% items$item_id[items$group == k & items$ok], ])
  compared <- dif_comparable(items[items$ok, ], groups)
  if (is.null(compared)) {
    return(out)
  }
  # 一組不能校準（例如題分成沒有共同題的兩片）：沒有每題的 DIF（lapply：錯誤保留 kits 的類別）
  steps <- tryCatch(lapply(by_group, group_steps, compared$common), kits_error = function(e) {
    kits_inform("No item DIF: ", conditionMessage(e), class = "kits_message_fit")
  })
  if (is.null(steps)) {
    return(out)
  }
  rel <- map(steps, \(st) st$beta - mean(st$beta))   # 每組的原點不同：減去這組所有步驟的平均
  out <- data.frame(item_id = steps[[ref]]$item_id, dif = rel[[ref]] - rel[[focal]]) |>
    summarise(dif = mean(dif), .by = item_id) |>
    arrange(sort_rank(item_id)) |>
    kits_table("kits_dif")
  attr(out, "overall") <- dif_overall(steps)
  attr(out, "se") <- dif_se(steps, out$item_id)
  attr(out, "groups") <- list(ref = ref, focal = focal, n = n)
  attr(out, "left_out") <- c(compared$only_one, compared$differ)
  out
}

# 比較的題：兩組都有、分數也相同（common）；只有一組有（或能校準）的題（only_one，有訊息）、
# 分數不同的題（differ，有警告）。可以比的題少於 2 題時 NULL（有訊息）
dif_comparable <- function(items, groups) {
  scores <- items |>
    mutate(group = match(group, groups)) |>
    select(item_id, group, scores) |>
    pivot_wider(names_from = group, values_from = scores, names_prefix = "scores_") |>
    arrange(sort_rank(item_id)) |>
    mutate(in_both = !is.na(scores_1) & !is.na(scores_2),
           differ = in_both & scores_1 != scores_2)
  only_one <- filter(scores, !in_both)$item_id
  differ <- filter(scores, differ)
  if (length(only_one) > 0) {
    kits_inform("Item(s) answered (with a 0 score, two scores and no gap) in one group only are ",
                "left out of the item DIF: ", paste(only_one, collapse = ", "), ".",
                class = "kits_message_data")
  }
  if (nrow(differ) > 0) {
    kits_warn("Item(s) with different scores in the two groups are left out: ",
              paste0(differ$item_id, " (", groups[1], ": ", differ$scores_1, "; ", groups[2], ": ", differ$scores_2, ")",
                     collapse = "; "),
              ".", class = "kits_warning_data")
  }
  common <- filter(scores, in_both, !differ)$item_id
  if (length(common) < 2) {
    kits_inform("Fewer than 2 items can be compared between the two groups: no item DIF.",
                class = "kits_message_data")
    return(NULL)
  }
  list(common = common, only_one = only_one, differ = differ$item_id)
}

# 一組的 CML：比較的題的步驟參數和它們的共變異
group_steps <- function(g, common) {
  fit <- cml_fit(g)
  keep <- fit$coef$item_id %in% common
  list(item_id = fit$coef$item_id[keep], beta = fit$est$beta[keep], V = fit$est$acov.beta[keep, keep, drop = FALSE])
}

# 整體檢定：比較的步驟先各自置中（減去這組比較的步驟的平均，C = I - 11'/S），再算拿掉第一個步驟的 Mahalanobis 距離。
# 每組各自校準，原點（所有步驟的和是 0）定在這組答的所有題上；有題沒有比較（只有一組有、分數不同）時，
# 兩組比較的步驟差了一個未知的常數。置中把這個常數拿掉，所以檢定只看相對的差（和 dif 相同）。
# 沒有拿掉題時，兩組比較的步驟的和都是 0，C d = d、C S C' = S，和 dexter 的 DIF_overall 相同。
dif_overall <- function(steps) {
  n <- length(steps[[1]]$beta)
  C <- diag(n) - matrix(1 / n, n, n)
  d <- drop(C %*% (steps[[1]]$beta - steps[[2]]$beta))
  S <- C %*% (steps[[1]]$V + steps[[2]]$V) %*% t(C)
  stat <- stats::mahalanobis(d[-1], rep(0, n - 1), S[-1, -1, drop = FALSE])
  list(stat = stat, df = n - 1L, p = stats::pchisq(stat, n - 1L, lower.tail = FALSE))
}

# 每題 dif 的 SE（給 plot() 的信賴區間）：置中的步驟差 C d 的共變異 C (V_1 + V_2) C'，
# 一題的 dif 是它的步驟的平均（權重 w），所以 SE = sqrt(w' C (V_1 + V_2) C' w)
dif_se <- function(steps, items) {
  n <- length(steps[[1]]$beta)
  C <- diag(n) - matrix(1 / n, n, n)
  S <- C %*% (steps[[1]]$V + steps[[2]]$V) %*% t(C)
  ids <- steps[[1]]$item_id
  se <- map_dbl(items, \(it) {
    w <- (ids == it) / sum(ids == it)
    sqrt(drop(t(w) %*% S %*% w))
  })
  stats::setNames(se, items)
}

dif_diagnostics <- function(x, n_na) {
  lr <- attr(x, "lrtest")
  ov <- attr(x, "overall")
  big <- lrt_largest(attr(x, "diff"))
  ids <- lrt_unobserved_ids(attr(x, "unobserved"))
  listed <- function(check, txt, what) {   # 有的話 note，列出前 6 個
    diag_flag(check, length(txt), length(txt) > 0, "note", paste(what, some_ids(txt, 6)), "none")
  }
  test_row <- function(check, what, stat, df, p, detail) {
    diag_row(check, p, if (p < 0.05) "note" else "ok",
             sprintf("%s = %.2f, df = %d, p = %s%s", what, stat, df, format.pval(p, digits = 3),
                     if (p < 0.05) detail else ""))
  }
  item_rows <- if (!is.null(ov)) list(
    test_row("overall DIF test", "chi-square", ov$stat, ov$df, ov$p,
             ": some items differ between the groups"),
    diag_row("largest |dif|", max(abs(x$dif)), "ok",
             sprintf("item %s (logits)", x$item_id[which.max(abs(x$dif))])),
    listed("items not compared", attr(x, "left_out"), "not in both groups with the same scores:"))
  diag_table(c(list(
    test_row("LR test", "LR", lr$LR, lr$df, lr$p,
             ": the steps differ between the groups; look at the size"),
    diag_row("largest |step difference|", big$value, "ok", paste(big$text, "(logits)"))),
    item_rows, list(
      listed("steps not observed in a group", ids$steps, "not estimated in that group:"),
      listed("items left out", ids$items, "of a group and of the joint calibration:"),
      diag_flag("person-waves without a group", n_na, n_na > 0, "note", "NA in `group`: left out",
                "none"),
      diag_row("groups", sum(lr$groups$n), "ok",
               paste0(lr$groups$group, ": ", lr$groups$n, collapse = ", ")))))
}

# 最大的 |步驟差|：值和說明（哪一步、哪兩組、Wald z）
lrt_largest <- function(d) {
  if (is.null(d) || nrow(d) == 0) {
    return(list(value = NA_real_, text = "no step in two groups"))
  }
  i <- which.max(abs(d$diff))
  list(value = abs(d$diff[i]), text = sprintf("%s:%d, %s minus %s, Wald z = %.1f", d$item_id[i],
                                              d$item_score[i], d$group_2[i], d$group_1[i], d$z[i]))
}


#' @rdname kits_dif
#' @param x A `kits_dif` (or its summary).
#' @param digits Decimals to print (only the printout is rounded; the object keeps full precision).
#' @param n Rows to print (see [print.kits_table()]); for `summary()`, the number of step
#'   differences to show (largest |difference| first).
#' @param ... In `kits_dif()`: not used (an argument the function does not have stops with an
#'   error).
#'   In `print()`: passed to `print()` for a data.frame.
#' @export
print.kits_dif <- function(x, digits = 2, n = NULL, ...) {
  need_digits(digits)
  lr <- attr(x, "lrtest")
  if (!is.null(lr)) {
    big <- lrt_largest(attr(x, "diff"))
    cat(sprintf("<kits_dif> groups by %s\n", lr$split))
    cat(strwrap(paste0("Groups: ", paste0(lr$groups$group, " (n = ", lr$groups$n, ")", collapse = ", ")),
                width = 100, exdent = 2), sep = "\n")
    cat(sprintf("LR test (Andersen): LR = %s, df = %d, p = %s\n", fmt(lr$LR, digits), lr$df,
                format.pval(lr$p, digits = 3)))
    cat(sprintf("Largest |step difference|: %s logits (%s)\n", fmt(big$value, digits), big$text))
  }
  g <- attr(x, "groups")
  if (!is.null(g)) {
    ov <- attr(x, "overall")
    cat(sprintf("Item DIF in logits: %s minus %s (positive = harder for %s, favors %s)\n", g$ref,
                g$focal, g$ref, g$focal))
    cat(sprintf("Overall DIF test: chi-square = %s, df = %d, p = %s\n", fmt(ov$stat, digits), ov$df,
                format.pval(ov$p, digits = 3)))
    print_rows(x, digits, n, ...)
    left <- attr(x, "left_out")
    if (length(left) > 0) cat("Left out (not in both groups with the same scores):", left, "\n")
  }
  if (NROW(attr(x, "unobserved")) > 0) {
    ids <- lrt_unobserved_ids(attr(x, "unobserved"))
    cat(sprintf("Not estimated in a group: %d step(s), %d item(s); see attr(x, \"unobserved\")\n",
                length(ids$steps), length(ids$items)))
  }
  checks_line(attr(x, "diagnostics"))
  invisible(x)
}

#' @rdname kits_dif
#' @param object A `kits_dif`.
#' @export
summary.kits_dif <- function(object, digits = 2, n = 10, ...) {
  need_digits(digits)
  need_rows(n)
  d <- as.data.frame(attr(object, "diff"))
  structure(list(object = object, largest = kits_table(utils::head(d[order(-abs(d$diff)), ], n)),
                 n_diff = nrow(d), digits = digits), class = "kits_dif_summary")
}

#' @rdname kits_dif
#' @export
print.kits_dif_summary <- function(x, ...) {
  digits <- x$digits
  print(x$object, digits = digits)
  cat(paste0("\nGroups (n = person-waves; n_info = informative; n_par = free step parameters; ",
             "logLik = conditional):\n"))
  print(format_table(attr(x$object, "lrtest")$groups, digits), row.names = FALSE, right = TRUE)
  cat(sprintf("\nLargest step differences (%d of %d; diff = beta in group_2 minus group_1):\n",
              nrow(x$largest), x$n_diff))
  print(format_table(x$largest, digits), row.names = FALSE, right = TRUE)
  u <- as.data.frame(attr(x$object, "unobserved"))
  if (nrow(u) > 0) {
    u <- summarise(u, n = n(), .by = group)
    cat(sprintf("\nSteps not estimated in a group: %s\n", paste0(u$group, ": ", u$n, collapse = ", ")))
  }
  cat("\nDiagnostics:\n")
  print(attr(x$object, "diagnostics"))
  invisible(x)
}


#' Combined DIF decision over two cohorts
#'
#' Decides the direction of each item in each cohort (beyond `cut` favors the focal group, below
#' `-cut` the reference group, else `"none"`) and combines the two: an item in one cohort only takes
#' that cohort's result, opposite signs with one side beyond `cut` are `"inconsistent"`, otherwise
#' the direction that exceeds `cut`. Both [kits_dif()] results must use the same `ref`.
#'
#' Details: `kits_algorithm("dif")` and tutorial_algorithms/05_checks.html.
#'
#' @param d1,d2 The [kits_dif()] results of the two cohorts (plain tables are read as
#'   `ref = "boy"`).
#' @param cut The cutoff for the decision (logits).
#' @param cohorts The names of the two cohorts, used in the output column names (`dif_<name>`,
#'   `dir_<name>`); by default the names of the two arguments, otherwise `"1"` and `"2"`.
#' @return A data.frame, one row per item: `item_id`, `dif_` and `dir_` for both cohorts, and the
#'   combined `DIF` (`"none"`, `"favors girls"`, `"favors boys"` or `"inconsistent"` with
#'   `ref = "boy"`), which [kits_item_table()] uses to colour the items.
#' @examples
#' d1 <- data.frame(item_id = c("a1", "a2", "a3", "a4"), dif = c(0.1, 0.5, 0.4, -0.5))
#' d2 <- data.frame(item_id = c("a1", "a2", "a3", "a5"), dif = c(0.2, 0.1, -0.2, 0.6))
#' kits_dif_combine(d1, d2)
#' @export
kits_dif_combine <- function(d1, d2, cut = 0.3, cohorts = NULL) {
  if (is.null(cohorts)) cohorts <- arg_names(substitute(d1), substitute(d2))
  rule <- function(r1, r2, dir1, dir2) {
    if (is.na(r2)) return(dir1)
    if (is.na(r1)) return(dir2)
    if (dir1 == "none" && dir2 == "none") return("none")
    if (r1 * r2 < 0) return("inconsistent")
    if (dir1 != "none") dir1 else dir2
  }

  need_cols(d1, c("item_id", "dif"), "d1")
  need_cols(d2, c("item_id", "dif"), "d2")
  abort_if(anyDuplicated(d1$item_id) || anyDuplicated(d2$item_id), "`d1` and `d2` must have one row per item.")
  abort_if(length(cohorts) != 2 || anyNA(cohorts) || cohorts[1] == cohorts[2], "`cohorts` needs two different names.")
  need_dif_cut(cut)
  g <- dif_sides(d1)
  g2 <- dif_sides(d2)
  abort_if(!identical(c(g$ref, g$focal), c(g2$ref, g2$focal)),
           sprintf("`d1` compares %s with %s, `d2` %s with %s. ", g$ref, g$focal, g2$ref,
                   g2$focal), "Use the same `ref` in both kits_dif() calls.")
  direction <- function(r) dif_direction(r, cut, g$focal, g$ref)
  out <- bind_rows(data.frame(item_id = d1$item_id, cohort = 1L, dif = d1$dif),
                   data.frame(item_id = d2$item_id, cohort = 2L, dif = d2$dif)) |>
    pivot_wider(names_from = cohort, values_from = dif, names_prefix = "dif_") |>
    arrange(sort_rank(item_id)) |>
    mutate(dir_1 = direction(dif_1),
           dir_2 = direction(dif_2),
           DIF = pmap_chr(list(dif_1, dif_2, dir_1, dir_2), rule)) |>
    select(item_id, dif_1, dif_2, dir_1, dir_2, DIF)
  names(out) <- c("item_id", paste0("dif_", cohorts), paste0("dir_", cohorts), "DIF")
  kits_table(out)
}

# 兩個參數的名字（kits_dif_combine() 的預設世代名）：都是不同的物件名字時用名字，否則 "1"、"2"
arg_names <- function(e1, e2) {
  if (is.symbol(e1) && is.symbol(e2) && !identical(e1, e2)) c(as.character(e1), as.character(e2)) else c("1", "2")
}

# DIF 的方向（kits_dif_combine()、kits_item_table()、plot() 共用）：大於 cut 對 focal 有利，小於 -cut
# 對 ref 有利，其餘 none；NA 留 NA。標籤用組的名字：boy、girl 寫成 "favors boys"、"favors girls"
dif_direction <- function(estimate, cut = 0.3, focal = "girl", ref = "boy") {
  labels <- dif_labels(focal, ref)
  case_when(is.na(estimate) ~ NA_character_,
            estimate > cut ~ labels[["up"]],
            estimate < -cut ~ labels[["down"]],
            .default = "none")
}

# 兩個方向的標籤：up 對 focal 有利，down 對 ref 有利
dif_labels <- function(focal, ref) {
  group <- c(up = focal, down = ref)
  plural <- group %in% c("boy", "girl")
  group[plural] <- paste0(group[plural], "s")
  stats::setNames(paste("favors", group), c("up", "down"))
}

# 一個 kits_dif() 的兩組（ref、focal）；沒記錄的表（例如自己做的）照 ref = "boy"
dif_sides <- function(dif) {
  g <- attr(dif, "groups")
  if (is.null(g)) list(ref = "boy", focal = "girl") else list(ref = g$ref, focal = g$focal)
}

need_dif_cut <- function(cut, arg = "cut") {
  need_number(cut, arg, "one number of at least 0 (logits)", lower = 0, finite = FALSE)
}
