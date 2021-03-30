#' Title
#'
#' @param N Abundance matrix with `nrow = n_sp = 3` and `ncol = n_patch`. Rows 1:3 correspond to species B, C, and P respectively.
#' @param v_p_dispersal A vector of dispersal probabilities of length = 3. in `igp_sim()`, this is controlled with the `p_dispersal` argument. Values should be from 0 to 1
#' @param theta Parameter controlling the distance decay function, of either length = 1 (same for all species) or length = 3 (vary by species). LArger values of theta result in a faster decline of dispersal distance.
#' @param dist_mat a distance matrix with `nrow = ncol = n_patch` and diagonal == 0.
#'
#' @details This function is used internally in `igp_sim()` to calculate the number of individuals of each dispersing from patch i to patch j. Emigrants from patch i are more likely to become Immigrants to patch j if the patches are closer together.
#'
#' @return A matrix of abundances after accounting for dispersal with `nrow = n_sp = 3` and `ncol = n_patch`. Values may be returned with decimal values, but will be converted to integers using `rpois()` within the simulation model.
#'
#' @export
#'
#' @examples
#' n_patch = 5
#' n_sp = 3
#' N = matrix(rpois(n_sp * n_patch, lambda = 10), nrow = n_sp, ncol = n_patch)
#' # sort coordinates so you can easily see increasing distances
#'  x_coord = sort(runif(5, 0, 5))
#'  y_coord = sort(runif(5, 0, 5))
#'  dist_mat = data.matrix(dist(cbind(x_coord, y_coord)))
#' disperal_n(N = N, v_p_dispersal = 0.25, theta = 1, dist_mat = dist_mat)
#'
#'
dispersal_n <- function(N, v_p_dispersal, theta, dist_mat){
  test1 = all(v_p_dispersal <=1) & all(v_p_dispersal>=0)
  if (test1 == FALSE){
    stop("in dispersal_n, v_p_dispersal must be between 0 and 1")
  }
  if(length(theta) != 1 & length(theta) !=3){
    stop("in dispersal_n(), Theta must be length 1 or 3")
  }
  if(length(theta) == 1){
    v_theta = rep(theta, 3)
  }
  if(length(theta) == 3){
    v_theta = theta
  }

  # dispersal matrices per species
  # species B
  m_b_dispersal <- data.matrix(exp(-v_theta[1] * dist_mat))
  diag(m_b_dispersal) <- 0
  # species C
  m_c_dispersal <- data.matrix(exp(-v_theta[2] * dist_mat))
  diag(m_c_dispersal) <- 0
  # species P
  m_p_dispersal <- data.matrix(exp(-v_theta[3] * dist_mat))
  diag(m_p_dispersal) <- 0

  # estimated Emmigrants
  m_e_hat <- N * v_p_dispersal
  v_e_sum <- rowSums(m_e_hat)
  # distribute Emigrants to new patches as Immigrants
  # Separate calculation for each species
  m_i_b_raw <- m_e_hat[1,] %*% m_b_dispersal
  m_i_c_raw <- m_e_hat[2,] %*% m_c_dispersal
  m_i_p_raw <- m_e_hat[3,] %*% m_p_dispersal
  # combine immigrant totals
  m_i_raw <- rbind(m_i_b_raw, m_i_c_raw, m_i_p_raw)


  v_i_sum <- rowSums(m_i_raw)
  v_i_sum[v_i_sum == 0] <- 1
  m_i_prob <- m_i_raw / v_i_sum
  m_i_hat <- m_i_prob * v_e_sum
  # abundance matrix, N, after dispersal occurs
  m_n_prime <- N + m_i_hat - m_e_hat

  # make sure all values are numeric and >=0
  m_n_prime[is.nan(m_n_prime)] <- 0
  m_n_prime[is.na(m_n_prime)] <- 0
  m_n_prime[m_n_prime <0 ] <- 0

  m_n_prime
}
