# PCM 的條件最大概似估計（CML），用 RTMB 從頭寫（內部用，給 kits_item()）。
#
# PCM：P(X_pi = k | theta_p) = exp(k theta_p - sum_{g <= k} beta_ig) / sum_h exp(h theta_p - sum_{g <= h} beta_ig)
# 原始分是能力的充分統計量，所以給定原始分，作答型態的機率只和 beta 有關：
#   P(x | r) = exp(-sum_i sum_{g <= x_i} beta_ig) / gamma_r
#   gamma_r = 原始分 r 的所有作答型態的 exp(-sum 達到的 beta) 加總（elementary symmetric function）
# gamma_r 用多項式卷積算：每一題是一個多項式 (1, exp(-beta_i1), exp(-beta_i1 - beta_i2), ...)，
# 題本的多項式是各題的乘積，第 r 項就是 gamma_r。
# 充分統計量：n_ig = 這題得到 g 分以上的人數；n_br = 題本 b 原始分 r 的人數。
#   log L = -sum_ig beta_ig n_ig - sum_b sum_r n_br log gamma_br
# 原點：所有 beta 的和是 0（和 dexter 相同）；有固定題時，固定題定原點，其他 beta 自由估。
# 共變異矩陣：Hessian 的反矩陣（A V A'）。
#
# 三步，一步一個函數（kits_item() 照這個順序呼叫）：
#   1 cml_stats()：充分統計量（kits_stats()，item$stats）
#   2 cml_estimate()：條件概似的最大值（kits_cml()，item$fit；計算在 algo_cml()，algo.R）
#   3 item_thresholds()：門檻（kits_thresholds()，item$thresholds；在 calibrate.R，計算在 algo_thresholds()）

# resp：長格式作答（只有要校準的題）；fixed：item_id、item_score、beta 的表，或 NULL
cml_fit <- function(resp, fixed = NULL) {
  cml_stats(resp, fixed) |>
    cml_estimate()
}

# ---- 1 充分統計量 ------------------------------------------------------------------------------------------
# n_ig（每個步驟，得到 g 分以上的人數）和 n_br（每個題本，原始分 r 的人數）。
# 題照 sort() 排，每題的步驟照順序；題本：人次照 visit_id 排（C locale），題本照出現的順序編號（和 dexter 相同）。
# left_out：沒放進校準的固定題（這些資料裡沒有 0 分或只有一種分數），門檻照固定的值算（kits_thresholds()）。
cml_stats <- function(resp, fixed = NULL, left_out = character()) {
  resp <- mutate(resp, item_id = as.character(item_id), visit_id = as.character(visit_id))
  items <- sort(unique(resp$item_id))
  fixed_all <- fixed
  if (!is.null(fixed)) fixed <- filter(fixed, item_id %in% items)
  m <- cml_step_counts(resp, items, fixed)
  persons <- person_sets(arrange(resp, visit_id))   # 每個人次照 visit_id 排（C locale）
  booklets <- booklet_table(persons, items, m)
  persons <- left_join(persons, select(booklets, set, pattern_id, max_raw), by = "set")
  steps <- step_table(resp, items, m)
  stats <- structure(list(items = items, steps = steps, patterns = select(booklets, -set),
                          raw = raw_table(persons, booklets), fixed = fixed,
                          left_out = if (length(left_out) > 0) filter(fixed_all, item_id %in% left_out),
                          n_persons = nrow(persons)),
                     class = "kits_stats")
  inp <- cml_inputs(stats)
  check_cml_design(inp$book_items, items, m, inp$step_item, inp$step_k, beta0_fixed = fixed)
  check_cml_existence(resp, persons, steps, fixed)
  stats
}

# 每題的步驟數：觀察到的最高分；固定題至少照固定參數的步驟數。
# 固定題這次出現更高的分數（舊的波次沒有人拿到）時，舊的步驟固定、新的步驟自由估（cml_constraint() 照步驟分）
cml_step_counts <- function(resp, items, fixed) {
  m <- summarise(resp, m = max(item_score), .by = item_id)
  m <- as.integer(pmax(m$m[match(items, m$item_id)], fixed_top(fixed, items), na.rm = TRUE))
  abort_if(sum(m) < 2, "At least 2 item steps are needed for a calibration.", class = "kits_error_fit")
  m
}

# 每個步驟（steps 的每一列：item_id、item_score = k）在 resp 的作答裡：分數低於 k 的人數（n_below）和
# 至少 k 的人數（n_above），照 steps 的順序。沒有人答的題兩個都是 0。
# 誰算（例如 CML 只算有資訊的人次）由呼叫的函數先篩 resp；規則（哪些步驟不行）也留在呼叫的函數：
# step_table()（n_ig）、check_cml_existence()、kits_drift() 的共同步驟（上下各幾人）。
step_reach <- function(resp, steps) {
  counts <- count(resp, item_id, score = item_score)
  steps |>
    select(item_id, item_score) |>
    left_join(counts, by = "item_id", relationship = "many-to-many") |>
    summarise(n_below = sum(n[score < item_score], na.rm = TRUE),
              n_above = sum(n[score >= item_score], na.rm = TRUE), .by = c(item_id, item_score))
}

# n_ig：每個步驟，得到 g 分以上的人數（一列一個步驟）
step_table <- function(resp, items, m) {
  data.frame(item_id = rep(items, m), item_score = sequence(m)) |>
    step_reach(resp = resp) |>
    transmute(item_id, item_score, n = n_above) |>
    as.data.frame()
}

# 題本：一個題目組合一列，照出現的順序編號；items 是題號（照 sort()），max_raw 是滿分
booklet_table <- function(persons, all_items, m) {
  number_booklets(persons) |>
    mutate(items = map(set, set_items),
           max_raw = map_int(items, \(its) sum(m[match(its, all_items)]))) |>
    select(pattern_id, items, n, max_raw, set) |>
    as.data.frame()
}

# n_br：每個題本，原始分 r = 0..滿分 的人數（沒有人的分數是 0）
raw_table <- function(persons, booklets) {
  counts <- count(persons, pattern_id, raw)
  booklets |>
    reframe(raw = 0:max_raw, .by = pattern_id) |>
    left_join(counts, by = c("pattern_id", "raw")) |>
    mutate(n = coalesce(n, 0L)) |>
    as.data.frame()
}

# 數值的部分要的向量：題、步驟、n_ig、每個題本的題（題的位置）和 n_br
cml_inputs <- function(stats) {
  items <- stats$items
  list(items = items,
       m = as.integer(tabulate(match(stats$steps$item_id, items), nbins = length(items))),
       step_item = stats$steps$item_id,
       step_k = stats$steps$item_score,
       n_ge = stats$steps$n,
       book_items = map(stats$patterns$items, \(its) match(its, items)),
       n_raw = unname(split(stats$raw$n,
                            factor(stats$raw$pattern_id, levels = stats$patterns$pattern_id))))
}

# ---- 2 條件最大概似 ----------------------------------------------------------------------------------------
# algo_cml()（algo.R）：cml_inputs() -> cml_constraint() -> cml_nll() -> cml_optimize()；這裡包成 kits_fit：
# 步驟參數的表（algo_cml() 的 se 叫 SE_beta）、共變異、題本的設計
cml_estimate <- function(stats) {
  cml <- algo_cml(stats)
  cf <- cml$coef
  names(cf)[names(cf) == "se"] <- "SE_beta"
  design <- data.frame(pattern_id = rep(stats$patterns$pattern_id, lengths(stats$patterns$items)),
                       item_id = unlist(stats$patterns$items, use.names = FALSE))
  structure(list(coef = cf,
                 est = list(beta = cml$coef$beta, acov.beta = cml$vcov, nll = cml$nll),
                 inputs = list(design = design, method = "CML", n_persons = stats$n_persons, left_out = stats$left_out),
                 optim = cml$optim),
            class = "kits_fit")
}

# beta = beta0 + A b：固定的步驟用固定的值，其他的估；沒有固定題時所有 beta 的和是 0
cml_constraint <- function(inp, fixed) {
  n_step <- length(inp$step_item)
  if (is.null(fixed) || nrow(fixed) == 0) {
    return(list(beta0 = numeric(n_step), A = rbind(diag(n_step - 1), -1)))
  }
  key <- paste(inp$step_item, inp$step_k)
  fixed_key <- paste(fixed$item_id, fixed$item_score)
  is_fixed <- key %in% fixed_key
  list(beta0 = ifelse(is_fixed, fixed$beta[match(key, fixed_key)], 0),
       A = diag(n_step)[, !is_fixed, drop = FALSE])
}

# 條件對數概似（RTMB 自動微分）：log L = -sum_ig beta_ig n_ig - sum_b sum_r n_br log gamma_br
# log gamma_r 用置中的 beta 算：gamma_r(beta) = exp(-r c) gamma_r(beta - c)，c 是這個題本的步驟的平均。
# 這樣長的題本、或固定題離原點很遠時，多項式的係數不會下溢（變成 0）或上溢（數學上完全相同）。
# 所有題本一起算（一列一個題本的矩陣 G，一欄一個原始分）：一題一題乘上它的多項式，只乘答這題的題本。
# 和一個題本一個迴圈（教學程式 cml.R 的寫法）的數學完全相同；一起算時 RTMB 記錄的速度快很多（題本多的時候）。
cml_nll <- function(inp, con) {
  design <- cml_book_design(inp)
  function(par) {
    beta <- con$beta0 + drop(con$A %*% par$b)
    centre <- cml_centre(beta, design)                                 # c_b：每個題本的步驟的平均
    G <- cml_gamma_matrix(beta, centre, inp, design)                    # gamma_br(beta - c_b)
    log_gamma <- log(G[design$cells]) - design$raw * centre[design$book]
    sum(beta * inp$n_ge) + sum(design$n * log_gamma)                    # -log L
  }
}

# 置中用的 c_b：每個題本的步驟參數的平均（CML 的概似和 algo_esf() 共用這一行；RTMB 也可以微分）
cml_centre <- function(beta, design) drop(design$has_step %*% beta) / design$n_steps

# 題本的設計（數值的部分要的索引）：每個題本有哪些題、哪些步驟、滿分，和 n_br 的每一格（題本、原始分）
cml_book_design <- function(inp) {
  n_book <- length(inp$book_items)
  has_item <- vapply(seq_along(inp$items), \(i) map_lgl(inp$book_items, \(its) i %in% its), logical(n_book))
  has_item <- matrix(has_item, nrow = n_book)
  has_step <- has_item[, rep(seq_along(inp$items), inp$m), drop = FALSE] * 1
  max_raw <- map_int(inp$book_items, \(its) sum(inp$m[its]))
  book <- rep(seq_len(n_book), max_raw + 1)
  raw <- sequence(max_raw + 1) - 1
  list(has_item = has_item, has_step = has_step, n_steps = rowSums(has_step), width = max(max_raw) + 1,
       steps_of = split(seq_along(inp$step_item), factor(inp$step_item, levels = inp$items)),
       cells = cbind(book, raw + 1), book = book, raw = raw, n = unlist(inp$n_raw, use.names = FALSE))
}

# gamma 的多項式，所有題本一起：G 一列一個題本，第 r + 1 欄是 z^r 的係數。
# 一開始每個題本是 1；第 i 題把答它的題本乘上 (1, exp(-(beta_i1 - c)), exp(-(beta_i1 - c) - (beta_i2 - c)), ...)：
# 乘上 z^k 的項就是整列往右移 k 欄。beta 可以是數值（algo_esf()）或 RTMB 的 advector（cml_nll()）。
cml_gamma_matrix <- function(beta, centre, inp, design) {
  "[<-" <- RTMB::ADoverload("[<-")   # 套件裡要明說：把 advector 放進數值矩陣
  G <- matrix(0, nrow(design$has_item), design$width)
  G[, 1] <- 1
  for (i in seq_along(inp$items)) {
    rows <- which(design$has_item[, i])
    G[rows, ] <- times_item_poly(G[rows, , drop = FALSE], beta[design$steps_of[[i]]], centre[rows])
  }
  G
}

# 一題的多項式乘上去（這些題本）：g(z) (1 + e_1 z + e_2 z^2 + ...)，e_k = exp(-sum_{h <= k} (beta_h - c))；
# 相乘用 poly_times()（pcm_core.R），寬度不變
times_item_poly <- function(g, beta_i, centre) {
  e <- list(1)
  cum <- 0
  for (k in seq_along(beta_i)) {
    cum <- cum + (beta_i[k] - centre)
    e[[k + 1]] <- exp(-cum)
  }
  poly_times(g, e, ncol(g))
}

# 最佳化：交給 fit_ad()（RTMB 給梯度和 Hessian）；共變異用最後那一點的 Hessian。
# 收斂的紀錄（給 kits_item() 的診斷）：nlminb 的代碼和訊息、最後的最大梯度、Hessian 是否正定
cml_optimize <- function(nll, con) {
  n_step <- length(con$beta0)
  if (ncol(con$A) == 0) {   # 每個步驟都固定：沒有東西要估
    return(list(beta = con$beta0, V = matrix(0, n_step, n_step), value = nll(list(b = numeric(0))),
                optim = list(convergence = NA_integer_, message = "", max_gradient = NA_real_, hessian_pd = NA)))
  }
  obj <- RTMB::MakeADFun(nll, list(b = numeric(ncol(con$A))), silent = TRUE)
  fit <- fit_ad(obj)
  H <- fit$hessian
  H_inv <- tryCatch(if (all(is.finite(H)) && rcond(H) > 1e-14) solve(H) else stop("reciprocal condition number below 1e-14"),
                    error = function(e) {
    kits_abort("The Hessian of the calibration is singular, so the SEs cannot be computed: ", conditionMessage(e),
               ". Check items with very few persons in a score category.", class = "kits_error_fit")
  })
  list(beta = con$beta0 + drop(con$A %*% fit$par), V = con$A %*% H_inv %*% t(con$A), value = fit$value,
       optim = fit[c("convergence", "message", "max_gradient", "hessian_pd")])
}

# 一個 RTMB 模型的最佳化（內部用；CML、同齡線的潛在常態共用）。
# obj 是 RTMB::MakeADFun() 的結果：模型只寫 -log L，梯度 obj$gr 和 Hessian obj$he 都由 RTMB 自動微分。
# 1. nlminb 只用梯度（準牛頓法；參數多時每一步也很快）。
# 2. 再用精確的 Hessian 走牛頓步 par - H^-1 g：nlminb 在很嚴的容許誤差下常提早停（false convergence），
#    牛頓步才走到最大值。-log L 沒有變大（超過捨入誤差 1e-9 + 1e-12 |-log L|）才走，最多 max_newton 步，
#    步子小於 1e-10 就停（最後的 Hessian 是走這一步前算的，差別小到可以忽略）。Hessian 最花時間，所以只算幾次。
# 走太遠溢位時，目標函數給很大的值、梯度的 NaN 給 0，讓 nlminb 退回來。
# 回傳 par、value（-log L）、hessian（最後一點），和收斂的紀錄 convergence、message、max_gradient、hessian_pd
fit_ad <- function(obj, max_newton = 5) {
  fn <- function(p) {
    v <- obj$fn(p)
    if (is.finite(v)) v else 1e300
  }
  gr <- function(p) {
    g <- obj$gr(p)
    g[!is.finite(g)] <- 0
    g
  }
  opt <- stats::nlminb(obj$par, fn, gr, control = list(rel.tol = 1e-14, eval.max = 2000, iter.max = 2000))
  par <- opt$par
  value <- obj$fn(par)
  H <- obj$he(par)
  for (step in seq_len(max_newton)) {
    if (any(!is.finite(H)) || rcond(H) < 1e-14) break   # 奇異的 Hessian：不走牛頓步（估不出來的設計在前面就停下來了）
    delta <- tryCatch(solve(H, drop(obj$gr(par))), error = function(e) NULL)
    if (is.null(delta)) break
    par_new <- par - delta
    v_new <- obj$fn(par_new)
    if (!is.finite(v_new) || v_new > value + 1e-9 + 1e-12 * abs(value)) break
    par <- par_new
    value <- v_new
    if (max(abs(delta)) < 1e-10) break
    H <- obj$he(par)
  }
  list(par = par, value = value, hessian = H, convergence = opt$convergence, message = opt$message,
       max_gradient = max(abs(obj$gr(par))), hessian_pd = is_pd(H))
}


#' @export
coef.kits_fit <- function(object, ...) {
  object$coef
}

#' @export
print.kits_fit <- function(x, ...) {
  cat(sprintf(paste0("<kits_fit> PCM by CML in RTMB: %d person-waves, %d items, %d steps, %d item ",
                     "pattern(s), -log L = %s\n"),
              x$inputs$n_persons, length(unique(x$coef$item_id)), nrow(x$coef),
              length(unique(x$inputs$design$pattern_id)), format(x$est$nll, digits = 8)))
  invisible(x)
}


# 校準的設計估得出來嗎（內部用，給 cml_fit()）。
# CML 以原始分為條件，所以作答組合（每個人答的題目）之間要有共同題，才能放在同一條量尺上：
# 題本經過共同題連成一片（沒有固定題時只能有一片；有固定題時，每一片都要有固定題來定原點）。
# 每個人只答一題時，原始分就是那一題的分數，條件概似沒有資訊。
check_cml_design <- function(book_items, items, m, step_item, step_k, beta0_fixed = NULL) {
  abort_if(all(lengths(book_items) == 1),
           "Every person-wave answered only one item, so CML has no information: given the raw score, ",
           "the response is known. CML needs at least two items per person-wave; calibrate the items together with ",
           "data where children answered several items.", class = "kits_error_data")
  group <- item_groups(book_items, length(items))
  fixed_items <- if (is.null(beta0_fixed)) character() else unique(as.character(beta0_fixed$item_id))
  parts <- split(items, group)
  free_parts <- parts[!map_lgl(parts, \(its) any(its %in% fixed_items))]
  too_many <- if (length(fixed_items) == 0) length(parts) > 1 else length(free_parts) > 0
  if (too_many) {
    show <- if (length(fixed_items) == 0) parts else free_parts
    txt <- map_chr(show, \(its) paste0("{", some_ids(its, 6), "}"))
    kits_abort("The item patterns (combinations of items answered) form ", length(parts),
               " separate groups that share no item: ",
               some_ids(txt, 5, sep = "; "), ". ",
               if (length(fixed_items) == 0) "CML cannot place separate groups on one scale. "
               else "A group without a fixed item has no origin. ",
               "Calibrate each group on its own and link them with kits_linking(), or add common items",
               if (length(fixed_items) > 0) " or fix an item of each group" else "", ".", class = "kits_error_data")
  }
}

# 題目分組：同一個題本的題連在一起（union-find，不用 igraph）；回傳每一題的組號
item_groups <- function(book_items, n_items) {
  parent <- seq_len(n_items)
  find <- function(i) {
    while (parent[i] != i) i <- parent[i]
    i
  }
  for (its in book_items) {
    if (length(its) < 2) next
    root <- find(its[1])
    for (i in its[-1]) {
      r <- find(i)
      if (r != root) parent[r] <- root
    }
  }
  vapply(seq_len(n_items), find, integer(1))
}

# 每個自由的步驟要有有限的 CML 估計：在原始分不是 0、也不是題本滿分的人次裡（只有他們帶有條件的資訊），
# 答這題的人要有人達到這個步驟，也要有人沒達到（Fischer, 1981 的必要條件）。
# 例如一題只有滿分的孩子通過：它的 beta 是無限大，CML 找不到最大值（以前會給 SE 上萬的數字）。
check_cml_existence <- function(resp, persons, steps, fixed) {
  inner <- filter(persons, raw > 0, raw < max_raw)$visit_id
  fixed_key <- if (is.null(fixed)) character() else paste(fixed$item_id, fixed$item_score)
  bad <- steps |>
    filter(!paste(item_id, item_score) %in% fixed_key) |>   # 固定的步驟不用估
    step_reach(resp = filter(resp, visit_id %in% inner)) |>
    rename(reach = n_above, miss = n_below) |>
    filter(reach == 0 | miss == 0) |>
    mutate(text = paste0(item_id, ":", item_score, " (", if_else(reach == 0, "reached by none", "missed by none"), ")"))
  abort_if(nrow(bad) > 0, "Item step(s) with no finite CML estimate: ", some_ids(bad$text, 10), ". ",
           paste0("Only person-waves whose raw score is neither 0 nor the top score of their ",
                  "item pattern carry information in CML, "),
           "and among them these steps are reached by none (or missed by none), for example an item passed only by children ",
           paste0("with a perfect score. Leave the item out, merge its top score with the score below, or fix it with ",
                  "`fixed_params`."), class = "kits_error_items")
}

# 對稱矩陣是否正定
is_pd <- function(H) {
  ev <- eigen((H + t(H)) / 2, symmetric = TRUE, only.values = TRUE)$values
  all(is.finite(ev)) && min(ev) > 0
}
