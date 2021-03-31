#' Calculate K based on environmental value
#'
#' @param env Vector of environmental values.
#' @param k_base Mean carrying capacity K.
#'
#' @details This function re-scales a vector of environmental values for patches to carrying capacities. the minimum carrying capacity = `0.5 * k_base` and the maximum carrying capacity = `1.5 * k_base`.
#'
#'  This function is designed to take `brnet()$df_patch$environment` as the `env` argument. See the `brnet()` function to control the headwater heterogeneity and spatial auto-correlation of environmental values between patches.
#'
#' @return Vector of carrying capacities
#'
#' @export
#'
#' @examples
#' env_to_k(env = c(0.5, 1, 1.5), k_base = 100)
#'
env_to_k <- function(env, k_base){
  if(any(is.na(c(env, k_base)))){
    stop("all inputs to `env_to_k()` need to be defined, \n one or more input values is 'NA' or 'NaN'" )
  }
  if(length(env) ==1){
    k = k_base
    k
  } else{
  k_min = 0.5 * k_base
  k_max = 1.5 * k_base
  k = round(((env - min(env))/ (max(env) - min(env))) *
            (k_max - k_min) + k_min)
  k}
}
