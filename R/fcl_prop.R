#' Title
#'
#' @param fcl Vector of FCL states, output from `get_fcl()`
#'
#' @details Designed to accept the `get_fcl()` output as an input vector. `fcl` should only be composed of values = `c(0, 1, 2, 2.5, 3)`
#'
#' @return vector of length 5, with each element being the proportion of patches in that FCL state.
#' @export
#'
#' @examples
#' # simulate abundance matrix
#' # nrow = n_sp = 3
#' # ncol = n_patch = 10
#' N <- matrix(rpois(3 * 10, lambda = 1), nrow = 3, ncol = 10)
#'
#' # make sure abundance of C and P = 0 when abundance of B = 0
#' N[,N[1,] == 0] <- 0
#' # calculate food chain length state
#' fcl_state <- get_fcl(N)
#' # calculate proportion of patches in a given state
#' fcl_prop(fcl = fcl_state)
#'
#'
fcl_prop <- function(fcl){
  prop_0 = length(fcl[fcl == 0]) / length(fcl)
  prop_1 = length(fcl[fcl == 1]) / length(fcl)
  prop_2 = length(fcl[fcl == 2]) / length(fcl)
  prop_25 = length(fcl[fcl == 2.5]) / length(fcl)
  prop_3 = length(fcl[fcl == 3]) / length(fcl)
  fcl_prop = (c("p0" = prop_0,
                "p1" = prop_1,
                "p2" = prop_2,
                "p2.5" = prop_25,
                "p3" = prop_3))
  fcl_prop
}
