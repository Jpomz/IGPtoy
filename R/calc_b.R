#' calculate parameter b which determines asymptotic level of prey recruitment
#'
#' @param r_max maximum per capita recruitment rate for basal species B
#' @param k carrying capacity
#'
#' @return parameter b
#' @export
#'
#' @examples
calc_b <- function(r_max, k){
  (r_max - 1) / k
}
