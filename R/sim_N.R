# helper function for development, not to be included in final package

sim_N <- function(n_patch){
  n_sp = 3
  N = matrix(rpois(n_patch * n_sp,
                   lambda = c(100, 50, 10)),
             nrow = n_sp,
             ncol = n_patch)
  N
}
