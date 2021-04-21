# helper function for development, not to be included in final package

sim_N <- function(n_patch = 5, fcl_5_states = FALSE){
  n_sp = 3
  if(fcl_5_states == FALSE){
  N = matrix(rpois(n_patch * n_sp,
                   lambda = c(100, 50, 10)),
             nrow = n_sp,
             ncol = n_patch)
  }
  if(fcl_5_states == TRUE){
    if(n_patch < 5) stop("n_patch must be >= 5 for all FCL states to be present")
    N = matrix(rpois(n_patch * n_sp,
                     lambda = c(100, 50, 10)),
               nrow = n_sp,
               ncol = n_patch)
    N[,1] <- 0
    N[2:3,2] <- 0
    N[2,3] <- 0
    N[3,4] <- 0
  }
  N
}
