#' Simulate IGP dynamics in 3 species communities through space and time
#'
#' @param n_patch single numeric value. Total number of patches in simulation
#' @param n_0 numeric value length should be either 1 or 3. starting abundances for species in each patch. If only one value is supplied it is assumed to be the same for all three species. Can designate starting values for each species with a vector of length 3. If NULL starting abundances for B, C, and P are 0.8, 0.5, and 0.25 * k_base, respectively.
#' @param dist_mat square distance matrix describing distance of all patches. If it is NULL, assumes a square landscape with `length = landscape_size`
#' @param landscape_size single numeric value. Length of a landscape on a side.
#' @param adjacency_matrix square matrix describing which patches are adjacent to one another.
#' @param k_function character string, one of: `c("patches-upstream", "random", "environment", NULL)`. `"patches-upstream"` calculates carrying capacity based on the number of nodes upstream. `"random"` assigns an integer value to each patch from a Poisson distribution with `lambda = k_base`. `"environment"` calculates K based on `df_patch$environment` value returned from `brnet()` function. `NULL` assumes that `k = k_base` in all patches.
#' @param n_upstream numeric vector of `length = n_patch`. used if `k_function = "patches-upstream"`
#' @param environment_value numeric vector of `length = n_patch`. used if `k_function = "environment"`
#' @param k_base single numeric value. Used to calculate carrying capacity, K. see ?k_function_internal for more.
#' @param k_c single numeric value. Constant used to calculate K when `k_function = "patches-upstream"`. Default = 10
#' @param k_min_exponent single numeric value. Minimum value of a uniform distribution for the exponent used to calculate K when `k_function = "patches_upstream"`
#' @param k_max_exponent single numeric value. Maximum value of a uniform distribution for the exponent used to calculate K when `k_function = "patches_upstream"`
#' @param p_dispersal numeric value (length should be one or 3). Probability of dispersal.
#' @param theta numeric value of length 1 or 3. Dispersal parameter describing dispersal capability of species. If length(theta) = 1 it is the same for all three species. Can be set for each species with vector of length = 3.
#' @param r_max single numeric value. Maximum reproductive number for the basal species, B, of the Beverton-Holt model.
#' @param ebc single Numeric value. Conversion efficiency of turning B biomass into new C biomass
#' @param ebp single Numeric value. Conversion efficiency of turning B biomass into new P biomass
#' @param ecp single Numeric value. Conversion efficiency of turning C biomass into new P biomass
#' @param alphabc single numeric value. Parameter controlling the predation of resource B by consumer C
#' @param betabc single numeric value. Parameter controlling the predation of resource B by consumer C
#' @param alphap single numeric value. Parameter controlling the predation of both resources by consumer P
#' @param betap single numeric value. Parameter controlling the predation of both resources by consumer P
#' @param P_pref single numeric value between 0-1 or NULL. Predator preference of B over C. if `NULL` (maybe NA?) calculates predator preference based on conversion efficiency rates and relative abundances of both prey species, i.e., Holling's "type 3" response.
#' @param s0 numeric value length should be either 1 or 3. Base survival probability for all (`length = 1`) or each (`length = 3`) species
#' @param disturb_type character string, one of: `c("point-source", "regional", or NULL)`. `"point-source"` applies a disturbance event with magnitude = `disturb_mag` to each patch with probability = `disturb_p` for each time step, and the magnitude of the disturbance decays moving downstream according to `disturb_rho` argument. `"regional"` applies a disturbance to the entire network in a given time step with probability = `disturb_p`. The magnitude of the disturbance varies based on the `environment_value`. `NULL` does not apply any disturbance events.
#' @param disturb_p single numeric value between 0-1. The probability that a disturbance occurs at a given time step.
#' @param disturb_mag single numeric value between 0-1. This controls the initial magnitude of impact given a disturbance occurs. For example, if a disturbance does occur, a value of 0.1, or 0.75 results in population abundances for each species being reduced by 10 or 75%, respectively. Assuming a single disturbance (i.e. adjacent or upstream patches not also disturbed) this values indicates the maximum disturbance magnitude a patch will be subjected to in a single time step.
#' @param disturb_rho Single numeric value >= 0. controls how quickly disturbance diminish with distance: 0 = no decay, all patches experience equal impact; 10 impact rapidly decays, with adjacent patches hardly being affected.
#' @param disturb_decay Single numeric value from 0 to 1. Indicates what percent of the disturbance remains at the next patch (i.e. 0.1, 0.5, 0.75, 1 = 10, 50, 75 and 100% of disturbance impacts in adjacent patch, respectively.
#' @param t Single numeric value. The number of time-steps to be saved.
#' @param plot_disturbance logical. If `TRUE` plots all patches that underwent a disturbance event.
#' @param plot_fcl logical. If `TRUE` plots the proportion of patches with a given food chain length through all time-steps.
#' @param plot_patch_dynamics logical. If `TRUE` plots population abundances through time for 5 random patches.
#'
#' @return `sp_dynamics` a data frame containing simulated IGP community dynamics.
#' @return `fcl` a data frame containing the proportion of patches that have a given food chain length for each time step.
#' @return `sim_params` a data frame which contains the parameter values used in the simulation. This data frame is a single row, unless one of the following argument is of length = 3: `p_dispersal`, `s0`, or `theta`, in which case the data frame will be 3 rows, corresponding to species B, C, and P, respectively.
#'
#' @importFrom dplyr %>% bind_rows filter pull select
#' @importFrom tidyr pivot_longer
#' @importFrom ggplot2 ggplot aes facet_grid facet_wrap geom_line geom_point label_both labeller labs geom_step theme_bw
#'
#' @export
#'
#' @examples
#' igp_sim()
#'
igp_sim <- function(n_patch = 20,
                    n_0 = NULL,
                    dist_mat = NULL,
                    landscape_size = 10,
                    adjacency_matrix = NULL,
                    k_function = NULL,
                    #c("patches-upstream", "random", "environment")
                    n_upstream = NULL,
                    environment_value = NULL,
                    k_base = 150,
                    k_c = 10,
                    k_min_exponent = 1.10,
                    k_max_exponent = 1.35,
                    #k = 500,
                    p_dispersal = 0.1,
                    theta = 1,
                    r_max = 2.5,
                    ebc = 2,
                    ebp = 2,
                    ecp = 2,
                    alphabc = 4,
                    alphap = 4,
                    betabc = 20,
                    betap = 20,
                    P_pref = 0.25, # preference of B over C
                    s0 = 0.75,
                    disturb_type = NULL, #c(NULL, "point-source", "regional")
                    disturb_p = 1e-4,
                    disturb_mag = 0.5,
                    disturb_rho = 1, # 2D habitats
                    disturb_decay = 0.75, # downstream
                    n_burnin = 200,
                    n_timestep = 1000,
                    plot_disturbance = FALSE,
                    # add disturb_type = "regional" subroutine
                    plot_fcl = FALSE,
                    plot_patch_dynamics = FALSE){
  # need to import functions:
  # dplyr:: bind_rows %>%
  # tidyr:: pivot_longer
  # ggplot2
  #library(tidyverse)

  param_df <- data.frame(
    k_base = k_base,
    k_function = ifelse(is.null(k_function), NA, k_function),
    r_max = r_max,
    alphabc = alphabc, betabc = betabc, ebc = ebc,
    alphap = alphap, betap = betap, ebp = ebp, ecp =ecp,
    P_pref = ifelse(is.null(P_pref), NA, P_pref),
    p_dispersal = p_dispersal, theta = theta,
    s0 = s0,
    disturb_type = ifelse(is.null(disturb_type), NA, disturb_type),
    disturb_p = disturb_p, disturb_mag = disturb_mag,
    disturb_rho = disturb_rho, disturb_decay = disturb_decay)
  if(nrow(param_df) == 1){
    row.names(param_df) <- "all_sp"
  }
  if(nrow(param_df) == 3){
    row.names(param_df) <- c("B", "C", "P")
  }

  n_sp = 3

  ## distance matrix ####
  dist_structure <- dist_mat_internal(dist_mat = dist_mat,
                    landscape_size = landscape_size,
                    n_patch = n_patch)
  dist_mat = dist_structure$dist_mat
  river_network_structure = dist_structure$river_network_structure
  if(river_network_structure == FALSE){
    environment_value <- rnorm(n_patch, mean = 0, sd = 1)
  }

  if(length(theta) == 1){
    v_theta = rep(theta, 3)
  }
  if(length(theta) == 3){
    v_theta = theta
  }
  # probability of dispersal ####
  if(length(p_dispersal) == 1){
    v_p_dispersal <- rep(x = p_dispersal, times = n_sp)
    message(
      paste("1 value of p_dispersal supplied:",
            p_dispersal,
            "for all three species"))
  }else{
    if(length(p_dispersal)!= n_sp)
      stop("length of p_dispersal should be 1  or n_sp")
    v_p_dispersal <- p_dispersal
    message(
      paste("3 values of p_dispersal supplied: B_dispersal =",
            p_dispersal[1],
            "C_dispersal =", p_dispersal[2],
            "and P_dispersal =", p_dispersal[3]))
  }

  # survival probability s0####
  if (length(s0) == 1){
    v_s0 <- rep(x = s0, times = n_sp)
    message(
      paste("1 value of s0 supplied:",
            s0,
            "for all three species"))
  }else{
    if(length(s0)!= n_sp)
      stop("length of s0 should be 1  or n_sp = 3")
    v_s0 <- s0
    message(
      paste("3 values of s0 supplied: B_s0 =", s0[1],
            "C_s0 =", s0[2],
            "and P_s0 =", s0[3]))
  }

  # dispersal distance ####
  # distance decay of dispersal, theta
  if(length(theta) != 1 & length(theta) !=3){
    stop("in igp_sim(), Theta must be length 1 or 3")
  }
  if(length(theta) == 1){
    message(
      paste("1 value of theta supplied: B_theta = C_theta = P_theta =",
            theta))
  }
  if(length(theta) == 3){
    message(
      paste("3 values of theta supplied: B_theta =", v_theta[1],
            "C_theta =", v_theta[2],
            "and P_theta =", v_theta[3]))
  }

  # dispersal matrices by species ####
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

  ## carrying capacity wrapper function ####
  b_k_list <- k_function_internal(
    k_function = k_function,
    k_c = k_c,
    k_min_exponent = k_min_exponent,
    k_max_exponent = k_max_exponent,
    k_base = k_base,
    r_max = r_max,
    n_upstream = n_upstream,
    n_patch = n_patch,
    river_network_structure = river_network_structure,
    environment_value = environment_value)
  k = b_k_list$k
  b = b_k_list$b

  # P_preference ####
  # predator preference fixed or variable
  # NULL or value from 0-1
  if(is.null(P_pref)){ # is.null(P_pref)
    fixed_P_pref = FALSE
    message("Predator preference varies with resource abundance")
  }
  if(is.numeric(P_pref) & length(P_pref == 1)){
    fixed_P_pref = TRUE
    message(paste("Predator preference is set at", P_pref))
  }

  ## initial community ####
  if(is.null(n_0)){
    N <- matrix(rpois(n = n_sp * n_patch,
                      lambda = c(k*0.8, k*0.5, k*0.25)),
                nrow = n_sp,
                ncol = n_patch)
    message(
      paste(
        "initial community abundance for B, C and P = 0.8, 0.5 and 0.25 *",
        k_base,
        "respectively"))
  } else{
    if(length(n_0) != 1 & length(n_0) != 3)
      stop("argument n_0 should be length 1, 3, or NULL")
    if(length(n_0) == 1){
      N <- matrix(rpois(n = n_sp*n_patch,
                        lambda = c(n_0, n_0, n_0)),
                  nrow = n_sp, ncol = n_patch)
      message(
        paste(
          "initial community abundance in all patches for all 3 species is",
          n_0))
    }
    if(length(n_0) == 3){
      N <- matrix(rpois(n = n_sp*n_patch,
                        lambda = c(n_0[1], n_0[2], n_0[3])),
                  nrow = n_sp, ncol = n_patch)
      message(paste("initial community abundance in all
                  patches for B =", n_0[1],
                  ", C =", n_0[2],
                  "and P =", n_0[3]),
              sep = "")
    }}

  # if basal species is 0 in a patch
  # make abundance of C and P = 0
  # make this a function "B_extinct"
  N[,N[1,] == 0] <- 0

  # result output ####
  output <- list()
  fcl_list <- list()
  # record starting conditions
  #fcl = get_fcl(N = N)
  #out = cbind(1:n_patch, t(N), i = 1, k, patch_extinction = 0, fcl)
  #colnames(out) <- c("patch", "B", "C", "P",
  #                   "time", "basal_k", "disturbance", "fcl")
  #output[[1]] <- as.data.frame(out)

  counter = 1 # for recording output at end of simulation loop

  n_sim = n_burnin + n_timestep
  ## simulation ####
  for (i in seq_len(n_sim)){
    # dispersal at beginning of time step
    m_n_prime <- dispersal_n(N = N,
                             v_p_dispersal = v_p_dispersal,
                             v_theta = v_theta,
                             dist_mat = dist_mat,
                             m_b_dispersal = m_b_dispersal,
                             m_c_dispersal = m_c_dispersal,
                             m_p_dispersal = m_p_dispersal)

    # turn species abundances into integers
    N <- matrix(rpois(n = n_sp * n_patch,
                      lambda = m_n_prime),
                nrow = n_sp,
                ncol = n_patch)

    # pop_sim() ####
    N = pop_sim(N = N,
            P_pref = P_pref,
            fixed_P_pref = fixed_P_pref,
            alphabc = alphabc,
            betabc = betabc,
            ebc = ebc,
            alphap = alphap,
            betap = betap,
            ebp = ebp,
            ecp = ecp,
            v_s0 = v_s0,
            b = b,
            k = k,
            r_max = r_max)

    # disturbance ####
    disturb_result = disturb_internal(
      N = N,
      adjacency_matrix = adjacency_matrix,
      dist_mat = dist_mat,
      environment_value = environment_value,
      disturb_type = disturb_type,
      disturb_p = disturb_p,
      disturb_mag = disturb_mag,
      disturb_rho = disturb_rho,
      disturb_decay = disturb_decay,
      river_network_structure = river_network_structure)

    N = disturb_result$N
    patch_extinction = disturb_result$patch_extinction

    # end of cycle ####
    # number of individuals  in patch x at time t
    N <-  matrix(rpois(n = n_sp * n_patch,
                       lambda = N),
                 nrow = n_sp,
                 ncol = n_patch)

    # if basal species is extinct in a patch
    # make abundance of C and P = 0
    N[,N[1,] == 0] <- 0

    if(i > n_burnin){

      # food chain length state
      fcl = get_fcl(N = N)

      # path_dynamics output
      out = cbind(1:n_patch,
                  t(N),
                  counter,
                  k,
                  patch_extinction,
                  fcl)
      colnames(out) <- c("patch",
                         "B",
                         "C",
                         "P",
                         "time",
                         "basal_k",
                         "disturbance", "fcl")
      output[[counter]] <- as.data.frame(out)
      fcl_list[[counter]] <- fcl_prop(get_fcl_state(N))
      counter = counter + 1
      #print(paste("loop iteration ", i))
    }
  }

  # make sure that these functions are properly imported
  dat <- dplyr::bind_rows(output)
  dat <- tidyr::pivot_longer(dat, 2:4, names_to = "species")

  fcl_df <- data.frame(dplyr::bind_rows(fcl_list), time = 1:n_timestep)
  fcl_df <- tidyr::pivot_longer(fcl_df, 1:5, names_to = "fcl_state")

  # Plots -------------------------------------------------------------------
# plot FCL ####
  # plot of food chain length
  if(plot_fcl == TRUE){
    fcl_plot <-
      ggplot(fcl_df,
             aes(y = value, x = time, color = fcl_state)) +
      geom_step() +
      facet_wrap(.~fcl_state) +
      theme_bw() +
      labs(y = "proportion of patches",
           title = "Food chain length") +
      NULL
    print(fcl_plot)
  }

  # plot patch dynamics ####
  # plot of patch dynamics
  if(plot_patch_dynamics == TRUE){

    plot_dat <- dat %>% filter(patch %in% sample(1:n_patch, 5))
    dist_dat <- dplyr::filter(plot_dat, disturbance == 1)

    if(nrow(dist_dat) >0){
      patch_plot <- ggplot(plot_dat,
                           aes(y = value,
                               x = time,
                               color = species))+
        geom_line() +
        geom_point(inherit.aes = FALSE,
                   data = dist_dat,
                   mapping = aes(x = time,
                                 y = -15),
                   shape = 2,
                   size = 2) +
        facet_wrap(.~patch, labeller = label_both) +
        labs(y = "Abundance",
             title = "Dynamics for 5 random patches") +
        theme_bw() +
        NULL
      print(patch_plot)
    } else {
      patch_plot <- ggplot(plot_dat,
                           aes(y = value,
                               x = time,
                               color = species))+
        geom_line() +
        facet_wrap(.~patch, labeller = label_both) +
        labs(y = "Abundance",
             title = "Dynamics for 5 random patches") +
        theme_bw() +
        NULL
      print(patch_plot)
    }
  }
  # plot disturbance ####
  # plot of patches with disturbances

  # add disturb_type = "regional" subroutine
  # plot 5 random patches, but show time point when disturbance occurred
  # title = "5 random patches experiencing regional disturbance"
  if(plot_disturbance == TRUE){
    dist_patch <- dplyr::filter(dat, disturbance == 1) %>%
      select(patch) %>% unique() %>% pull()
    dist_dat <- filter(dat, patch %in% dist_patch)
    dist_point <- dplyr::filter(dat, disturbance == 1)
    if(nrow(dist_dat) >0){
      dist_plot <- ggplot(dist_dat,
                          aes(y = value,
                              x = time,
                              color = species))+
        geom_line() +
        geom_point(inherit.aes = FALSE,
                   data = dist_point,
                   mapping = aes(x = time,
                                 y = -15),
                   shape = 2,
                   size = 2) +
        facet_wrap(.~patch, labeller = label_both) +
        labs(y = "Abundance",
             title = "All patches experiencing disturbances") +
        theme_bw() +
        NULL
      print(dist_plot)
    } else {message("No disturbances occured")}
  }

  # return ####
  return(list(sp_dynamics = dat,
              fcl_state_prop = fcl_df,
              sim_params = param_df))
}
