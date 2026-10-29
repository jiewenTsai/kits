# dplyr::bind_rows() 疊 kits 的表時，也檢查量尺（和 rbind.kits_table() 相同）。
# bind_rows() 用 vctrs 疊表：先用 vec_ptype2() 求共同的型別，再用 vec_cast() 把每張表轉過去。
# vctrs 只看類別的第一個（不看繼承），所以每一對類別都註冊一個方法（.onLoad() 裡註冊）。
# * 兩張 kits 的表、量尺不同：停下來（kits_error_scale）。
# * 同一種 kits 的表：結果留著類別和共同的量尺（診斷和來源拿掉，和 dplyr 的動詞相同）。
# * kits 的表和一般的 data.frame（或 tibble）：警告（kits_warning_scale），結果是 kits 的表。
# * 不同種的 kits 的表（例如 kits_subj 和門檻表）：結果是一般的 data.frame。

# kits 的表的類別（類別的第一個）
KITS_TABLE_CLASSES <- c("kits_table", "kits_subj", "kits_subj_table", "kits_pattern_table",
                        "kits_norms", "kits_itemfit", "kits_dif",
                        "kits_vertical", "kits_precision", "kits_diagnostics", "kits_check")

# 一般的表：沒有量尺
PLAIN_TABLE_CLASSES <- c("data.frame", "tbl_df")

# 兩張 kits 的表的量尺要相同（沒有記錄的不比）
check_bind_scales <- function(x, y) {
  sx <- scale_of(x)
  sy <- scale_of(y)
  need_same_scale(sx, sy, "The tables are on different scales: ", sx, " and ", sy,
                  ". Convert them to one scale with kits_convert() first.")
  abort_if(rp_of(x) != rp_of(y),
           "The tables hold thresholds at different response probabilities (RP): ", 100 * rp_of(x), " and ", 100 * rp_of(y),
           ".")
}

# kits 的表的屬性（量尺、beta 的量尺、RP）留在共同的型別上；量尺從有記錄的那一張拿
kits_attrs <- function(out, x, y = NULL) {
  for (a in c("scale", "beta_scale", "rp")) {
    v <- attr(x, a)
    if (is.null(v) && !is.null(y)) v <- attr(y, a)
    attr(out, a) <- v
  }
  attr(out, "links_without_se") <- combine_no_se_records(list(x, y))   # 連結誤差沒有帶下來的連結：兩張表的合起來
  attr(out, "raw_theta") <- combine_raw_theta(list(x, y))
  out
}

# 共同的型別：同一種 kits 的表
ptype2_same <- function(x, y, ...) {
  check_bind_scales(x, y)
  out <- vctrs::df_ptype2(x, y, ...)
  class(out) <- class(x)
  kits_attrs(out, x, y)
}

# 共同的型別：不同種的 kits 的表，給一般的 data.frame
ptype2_mixed <- function(x, y, ...) {
  check_bind_scales(x, y)
  vctrs::df_ptype2(x, y, ...)
}

# 共同的型別：kits 的表和一般的表（量尺沒有檢查，警告）
ptype2_plain <- function(x, y, ..., kits_first = TRUE) {
  k <- if (kits_first) x else y
  kits_warn("A table without a recorded scale is bound to a table on scale ", scale_label(scale_of(k)),
            ": its scale was not checked.", class = "kits_warning_scale")
  out <- vctrs::df_ptype2(x, y, ...)
  class(out) <- class(k)
  kits_attrs(out, k)
}

# 轉到 kits 的表：欄照共同的型別，屬性照 to
cast_to_kits <- function(x, to, ...) {
  out <- vctrs::df_cast(x, to, ...)
  class(out) <- class(to)
  kits_attrs(out, to)
}

# 轉到一般的 data.frame（或 tibble）
cast_to_plain <- function(x, to, ...) {
  if (inherits(to, "tbl_df")) vctrs::tib_cast(x, to, ...) else vctrs::df_cast(x, to, ...)
}

# vctrs 的方法（.onLoad() 呼叫）：每一對類別各一個
register_vctrs_methods <- function() {
  ns <- asNamespace("vctrs")
  reg <- function(generic, x_class, y_class, fun) {
    registerS3method(generic, paste(x_class, y_class, sep = "."), fun, envir = ns)
  }
  for (a in KITS_TABLE_CLASSES) {
    for (b in KITS_TABLE_CLASSES) {
      reg("vec_ptype2", a, b, if (a == b) ptype2_same else ptype2_mixed)
      reg("vec_cast", a, b, if (a == b) cast_to_kits else cast_to_plain)
    }
    for (p in PLAIN_TABLE_CLASSES) {
      reg("vec_ptype2", a, p, function(x, y, ...) ptype2_plain(x, y, ..., kits_first = TRUE))
      reg("vec_ptype2", p, a, function(x, y, ...) ptype2_plain(x, y, ..., kits_first = FALSE))
      reg("vec_cast", a, p, cast_to_kits)   # 一般的表轉成 kits 的表（vec_cast.<to>.<from>）
      reg("vec_cast", p, a, cast_to_plain)
    }
  }
  invisible()
}

.onLoad <- function(libname, pkgname) {
  register_vctrs_methods()
}
