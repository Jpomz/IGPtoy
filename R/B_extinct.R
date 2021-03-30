#' Set abundances to 0 when B goes extinct
#'
#' @param N Community abundance matrix with `nrow = n_sp = 3` and `ncol = n_patch`
#'
#' @details This is an internal function which sets higher trophic levels (species C and P), to 0 in a patch when the abundance of species B is 0.
#' @return Abundance matrix, N.
#' @export
#'
#' @examples
#' n_sp = 3
#' n_patch = 5
#' N <- matrix(rpois(n = n_sp*n_patch, lambda = 100), nrow = n_sp, ncol = n_patch)
#' N
#' # make a few Row 1 == 0
#' N[1, 1:2] <- 0
#' N
#' B_extinct(N)
#'
B_extinct <- function (N){
  N[,N[1,] == 0] <- 0
  N
}
