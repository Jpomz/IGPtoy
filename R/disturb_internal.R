#' Determine if a disturbance occured and reduce population abundances in a time step.
#'
#' @param N Abundance matrix with `nrow = n_sp = 3` and `ncol = n_patch`.
#' @param adjacency_matrix Square matrix (`n_patch * n_patch`) describing which patches are adjacent to each other. This is used for branching river networks, and is designed to take the `adjacency_matrix` output from `brnet()` function as an input.
#' @param dist_mat Distance matrix indicating the distance between patches. This is designed to take the `distance_matrix` output from `brnet()` as an input. If it is not supplied (`dist_mat = NULL`), `n_patch` are assumed to be randomly dispersed in a square, 2D landscape with dimensions `landscape_size * landscape_size`
#' @param environment_value Vector of numbers describing an environmental condition in patch `i`. This vector will be scaled from 0 to 1 using the inverse logit function. This is designed to take the output of `df_patch$environment` from `brnet()` as an input.
#' @param disturb_type Character vector of one of `c("point-source", "regional", or NULL)`. "point-source" applies disturbances randomly to individual patches. If it is a river network, the disturbance travels downstream and decays according to the `disturb_decay`  If it is a 2D network, the disturbance decays with increasing distance from disturbed patch, and is controlled by `disturb_rho` argument
#' @param disturb_p probability of disturbance occurring. If `disturb_type = "point-source"`, whether or not a disturbance occurs is determined independently for each patch via `rbinom(n = n_patch, size = 1, prob = disturb_p)`. If `disturb_type = "regional"`, disturbance for the entire meta-community is determined in each time step via `rbinom(n = 1, size = 1, prob = disturb_p)`.
#' @param disturb_mag numeric value (0 to 1) controlling the initial magnitude of impact given a disturbance occurs. 0.1, 0.75 = population abundances reduced by 10 or 75%, respectively. Assuming a single disturbance (i.e. adjcanet or upstream patches not also disturbed) this values indicates the maximum disturbance magnitude a patch will be subjected to in a single time step.
#' @param disturb_rho integer >= 0. controls how quickly disturbance diminish with distance: 0 = no decay, all patches experience equal impact; 10 impact rapidly decays, with adjacent patches hardly being affected.
#' @param disturb_decay numeric from 0 to 1. Indicates what percent of the disturbance remains at the next patch (i.e. 0.1, 0.5, 0.75, 1 = 10, 50, 75 and 100% of disturbance impacts in adjacent patch, respectively.
#'
#' @return List with two elements. `N` = abundance matrix after accounting for disturbances. Values don't necessarily have to be integers. Numeric-double values will be converted to integer with `rpois()` in the simulation function. `patch_extinction` is an integer vector of `length = n_patch` indiciating if a disturbance did happen (1) or did not happen (0).
#'
#' @export
#'
#' @examples
disturb_internal <- function(N, adjacency_matrix, dist_mat,
                             environment_value,
                             disturb_type, disturb_p, disturb_mag,
                             disturb_rho, disturb_decay){
  n_patch <-  ncol(N)
  # no disturbances
  if(is.null(disturb_type)){
    N <- N
    patch_extinction <-  rep(0, n_patch)
    } else {
  # point-source disturbance
  if(disturb_type == "point-source"){
    if(is.null(adjacency_matrix)){ # move message outside of for loop
      #message("No adjacency matrix supplied, assuming disturbance in 2D habitat")}
  # point source in 2D habitats
  patch_extinction <- rbinom(n = n_patch,
                             size = 1, prob = disturb_p)
  patch_disturb_id <- which(patch_extinction==1)
  if(length(patch_disturb_id) == 0){
    N <-  N
  }
  if(length(patch_disturb_id) != 0){
    disturb_m <- data.matrix(exp(-disturb_rho * dist_mat))
    disturb_m <- as.matrix(disturb_m[,patch_disturb_id])
    N <- N*(1 - rowSums(disturb_m)*disturb_mag)[col(N)]
  }
  # point source disturbance in river network
  if(!is.null(adjacency_matrix)){
    #message("adj matrix not null")
    # Disturbance "decay" down stream
    patch_extinction <- rbinom(n = n_patch,
                               size = 1,
                               prob = disturb_p)
    patch_disturb_id <- which(patch_extinction==1)

    m_adj_up <- adjacency_matrix
    m_adj_up[lower.tri(m_adj_up)] <- 0

    # vector to hold disturbance magnitudes
    disturb_dummy <- disturb_v <-  rep(0, n_patch)
    # "source" disturbance magnitude
    disturb_dummy[patch_disturb_id] <- disturb_mag

    for(np in n_patch:1){
      disturb_v <- disturb_v + disturb_dummy
      disturb_dummy <- m_adj_up %*% (disturb_decay * disturb_dummy)
    }

    disturb_v[disturb_v >=1] <- 1
    disturb_v <- as.vector(disturb_v)

    N <- N*(1 - disturb_v)[col(N)]
    #message("end of disturb loop")
  }
    }
  }
    if(disturb_type == "regional"){
      if(is.null(environment_value)){
        stop("disturb_type = regional but no `environment_value` supplied")
      }
      patch_extinction <- rbinom(n = 1, size = 1, prob = disturb_p)
      if( patch_extinction == 1){
        # scale environment_value to inverse logit scale (0 to 1)
        env_logit <- exp(environment_value) / (1 + exp(environment_value))
        disturb_v <- env_logit * disturb_mag
        N <- N*(1 - disturb_v)[col(N)]
        patch_extinction <- rep(patch_extinction, n_patch)
      }
    }
  }
  N[N<0] <- 0
  return(list(N = N, patch_extinction = patch_extinction))
}

