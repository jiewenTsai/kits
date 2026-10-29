# D 連結的錨題（KIT 題目對到 GSED 題目，和那個 GSED 題目在 dscore key 裡的 tau）。
# 表存在 R/sysdata.rda（kits_anchor_table，37 列：gsed2406 29 題、gsed2510 8 題），只給一個 key。

#' The anchor items of the link to D (GSED), for one dscore key
#'
#' The anchor items of KIT for the link to the D score of the Global Scales for Early Development (GSED):
#' each row is a KIT item, the GSED item it is matched to, and the difficulty `tau` of that GSED
#' item in a key of the dscore package (van Buuren et al., 2025). Give the table to
#' [kits_linking()] as `to`; use one key per link (the link records it).
#' Details: `vignette("kits-flow", package = "kits")`, section 3.
#'
#' @param key The dscore key: `"gsed2406"` (default, 29 anchors) or `"gsed2510"` (8 anchors,
#'   another D metric; gives a warning of class `kits_warning_anchors`).
#' @return A data.frame with one row per anchor: `item_id` (the KIT item), `gsed_item` (the GSED
#'   item), `key`, `tau` (on the D scale), `note` (the matching decision) and `domain` (`"lan"`,
#'   `"motor"` or `"cog"`).
#' @references
#' van Buuren, S., Eekhout, I., McCray, G., Lancaster, G. A., Waldman, M. R., McCoy, D. C., Gladstone, M., et al. (2025).
#' Enhancing comparability in early child development assessment with the D-score.
#' *International Journal of Behavioral Development, 49*(4), 348-364. \doi{10.1177/01650254241294033}
#'
#' van Buuren, S., Eekhout, I., Huizing, A., & Seiden, J. (2025). *dscore: D-score for child development*.
#' R package. <https://CRAN.R-project.org/package=dscore>
#' @seealso [kits_linking()]
#' @examples
#' head(kits_anchors())
#'
#' # A link to D from a simulated calibration whose items carry KIT item ids
#' anchors <- kits_anchors()
#' ids <- anchors$item_id[1:8]
#' betas <- as.list((anchors$tau[1:8] - 55) / 10)
#' names(betas) <- ids
#' X <- kits_sim(rnorm(600, 0, 1.5), betas, seed = 1)
#' item <- kits_item(X, scale = "sim M03")
#' kits_linking(item, anchors)
#' @export
kits_anchors <- function(key = c("gsed2406", "gsed2510")) {
  chosen <- match_option(key, c("gsed2406", "gsed2510"), "key")
  if (chosen == "gsed2510") {
    kits_warn("Key gsed2510 (dscore: \"In preparation\") covers only the GSED SF and LF items and has another D metric ",
              "than gsed2406: only 8 of the 29 KIT anchors have a tau in it (2 language, 6 motor, no cognitive), ",
              "and its D values are not interchangeable with those of gsed2406. See ?kits_anchors.",
              class = "kits_warning_anchors")
  }
  out <- filter(kits_anchor_table, key == chosen) |>
    mutate(domain = anchor_domain(item_id))
  rownames(out) <- NULL
  out
}

# 錨題的領域，從 KIT 的題號看：lan 語言、body 動作（粗動作和精細動作）、cog 認知
anchor_domain <- function(item_id) {
  case_when(startsWith(item_id, "lan") ~ "lan",
            startsWith(item_id, "body") ~ "motor",
            startsWith(item_id, "cog") ~ "cog",
            .default = NA_character_)
}

# 錨題表的 key（kits_anchors() 的表有 key 欄，只有一個值時記在連結裡）；沒有時 NULL
anchor_key <- function(anchors) {
  if (!"key" %in% names(anchors)) {
    return(NULL)
  }
  keys <- unique(as.character(anchors$key))
  if (length(keys) == 1 && !is.na(keys)) keys else NULL
}
