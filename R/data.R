#' Simulated KIT data for teaching
#'
#' Simulated versions of the two Kids in Taiwan (KIT) cohorts, language and motor only, for
#' examples and tutorials. No row comes from a real child.
#' The response columns are those of the real KIT response file, renamed as everywhere in kits:
#' `visit_id` (KIT: `person_id`), `subj_id` (KIT: `pid`) and `cohort` (KIT: `sample`). On the real data,
#' run `rename(visit_id = person_id, subj_id = pid, cohort = sample)` once and the same code runs.
#'
#' Each cohort is a list of two tables:
#'
#' * `responses`: the responses, in long format (see below).
#' * `anchors`: the anchor key for the items of this cohort, from [kits_anchors()] (key `gsed2406`):
#'   KIT items matched to GSED items, with their GSED difficulty `tau`. Use it to link to D with
#'   `kits_linking(items, anchors)`. `kitm03` has 27 anchors, `kitm36` 12.
#'
#' The data were simulated from aggregates of the real KIT data (`data-raw/kit_demo.R`); the waves
#' at 60 and 72 months of `kitm03` go beyond the real M03 follow-up (48 months) and are not observed.
#' Details: tutorial_algorithms/01_calibration.html.
#'
#' @format Each is a list with two elements.
#'
#' `responses` is a tibble with one row per child x wave x item answered:
#' \describe{
#'   \item{visit_id}{Visit (person-wave) id: one child at one wave, for example `"M03_c55df9b_M12"`.}
#'   \item{subj_id}{Child id: 7 random hex characters, like a short git hash, for example `"c55df9b"`.}
#'   \item{cohort}{Cohort: `"M03"` (3 to 72 months) or `"M36"` (36 to 96 months).}
#'   \item{sex}{`"boy"` or `"girl"`.}
#'   \item{wave_m}{Wave, in months.}
#'   \item{age_m}{Age at the visit, in months.}
#'   \item{item_id}{Item id, as in KIT.}
#'   \item{item_score}{Score, 0 to the maximum of the item.}
#'   \item{domain}{`"lan"` (language) or `"motor"`.}
#' }
#'
#' `anchors` is a data.frame with the columns of [kits_anchors()]: `item_id`, `gsed_item`, `key`,
#' `tau`, `note`.
#' @examples
#' library(dplyr)
#' kitm03$responses |>
#'   distinct(visit_id, wave_m) |>
#'   count(wave_m)
#' kitm03$anchors |>
#'   count(key)
#' @name kitm03
"kitm03"

#' @rdname kitm03
"kitm36"
