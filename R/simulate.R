#' Simulate PCM responses
#'
#' Simulates the scores (starting at 0) of persons with abilities `theta` on PCM items with step
#' parameters `betas`, for small examples.
#'
#' @param theta Ability of each person.
#' @param betas A list with one element per item: the names are item ids, each element holds the
#'   step parameters of that item (one beta for a dichotomous item).
#' @param prefix Start of the person ids (`"y0001"` for `prefix = "y"`); give each cohort its own.
#'   `NULL` (default): `"p"`, 4 hex digits made from the simulated data, and `"_"` (`"p3fa2_0001"`).
#' @param seed Random seed; when given, your current random state is left unchanged. `NULL`
#'   (default) uses the current random state, so you can call `set.seed()` first.
#' @return A matrix: one row per person and one column per item; row names are visit_id, column
#'   names are item ids. It can be passed directly to the functions in this package.
#' @examples
#' betas <- list(a1 = 0, a2 = c(-1, 1))
#' kits_sim(theta = c(-2, 0, 2), betas = betas)
#'
#' # Two cohorts: give each its own prefix, so the ids differ
#' young <- kits_sim(rnorm(5), betas, prefix = "y", seed = 1)
#' old <- kits_sim(rnorm(5, 1), betas, prefix = "o", seed = 2)
#' rownames(rbind(young, old))
#' @export
kits_sim <- function(theta, betas, prefix = NULL, seed = NULL) {
  abort_if(is.null(names(betas)) || any(names(betas) == "") || anyDuplicated(names(betas)),
           "`betas` must be a named list: one element per item, named by different item ids.")
  abort_if(length(theta) == 0 || !is.numeric(theta) || any(!is.finite(theta)),
           "`theta` must be one or more finite numbers.")
  abort_if(any(!map_lgl(betas, \(b) is.numeric(b) && length(b) > 0 && all(is.finite(b)))),
           "Each element of `betas` must be one or more finite numbers.")
  abort_if(!is.null(prefix) && (!is.character(prefix) || length(prefix) != 1 || is.na(prefix)),
           "`prefix` must be one string, for example \"y\", or NULL.")

  # 每個人一列的類別權重 exp(eta - 最大值)，來自 pcm_kernel()（pcm_core.R）；sample() 自己標準化
  one_item <- function(item) {
    ex <- pcm_kernel(theta, betas[[item]])$ex
    scores <- seq_len(ncol(ex)) - 1
    vapply(seq_along(theta), function(r) sample(scores, 1, prob = ex[r, ]), numeric(1))
  }

  simulate <- function() vapply(names(betas), one_item, numeric(length(theta)))
  out <- if (is.null(seed)) simulate() else with_seed(seed, simulate())
  out <- matrix(as.integer(out), nrow = length(theta), dimnames = list(NULL, names(betas)))
  if (is.null(prefix)) prefix <- paste0("p", substr(hash_of(list(theta, betas, out)), 1, 4), "_")   # 從資料來的前綴：不同的呼叫不撞名
  rownames(out) <- sprintf("%s%04d", prefix, seq_along(theta))
  out
}
