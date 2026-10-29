# 有類別的錯誤、警告、訊息（內部用）。使用者可以用 tryCatch(..., kits_error_scale = function(e) ...) 接。
# 訊息的各段像 stop() 一樣接起來：kits_abort("a", b, "c", class = "kits_error_data")。
# 錯誤：kits_error_input、kits_error_data、kits_error_items、kits_error_scale、kits_error_fit
# 警告：kits_warning_convergence、kits_warning_scale、kits_warning_norms、kits_warning_data、kits_warning_anchors
# 訊息：kits_message_data、kits_message_scale、kits_message_fit

kits_condition <- function(msg, class, type) {
  structure(class = c(class, paste0("kits_", type), type, "condition"),
            list(message = paste0(msg, collapse = ""), call = NULL))
}

kits_abort <- function(..., class = "kits_error_input") {
  stop(kits_condition(paste0(...), class, "error"))
}

# 檢查：cond 是 TRUE 時停下來，訊息和類別同 kits_abort()（輸入檢查的 if (...) { kits_abort(...) } 寫成一行）
abort_if <- function(cond, ..., class = "kits_error_input") {
  if (cond) kits_abort(..., class = class)
}

kits_warn <- function(..., class) {
  warning(kits_condition(paste0(...), class, "warning"))
}

kits_inform <- function(..., class) {
  message(kits_condition(paste0(paste0(...), "\n"), class, "message"))
}

# 函數沒有的參數（... 只是為了放舊的名字）：說清楚是哪幾個
no_extra_args <- function(..., fun) {
  if (...length() == 0) {
    return(invisible())
  }
  given <- ...names()
  if (is.null(given)) given <- rep("", ...length())
  named <- given[given != ""]
  unnamed <- sum(given == "")
  what <- c(if (length(named) > 0) paste0("`", named, "`", collapse = ", "),
            if (unnamed > 0) sprintf("%d unnamed value(s)", unnamed))
  kits_abort(sprintf("%s() got argument(s) it does not have: %s. See ?%s.", fun, paste(what, collapse = " and "), fun),
             class = "kits_error_input")
}

