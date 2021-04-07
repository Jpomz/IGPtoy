library(devtools)
library(profvis)
library(mcbrnet)
#library(microbenchmark)

load_all()

# microbenchmark(
#   igp_sim(),
#   times = 10
# )


profvis({
  igp_sim(plot_disturbance = FALSE,
          plot_fcl = FALSE,
          plot_patch_dynamics = FALSE)
})

profvis({
  igp_sim()
})

# n_patch = 20
# landscape_size = 10
# N = matrix(rpois(3 *n_patch, lambda = 100),
#            nrow = 3,
#            ncol = n_patch)
# x_coord = runif(n_patch, 0, landscape_size)
# y_coord = runif(n_patch, 0, landscape_size)
# dist_mat = data.matrix(dist(cbind(x_coord, y_coord),
#                             diag = TRUE, upper = TRUE))
#
# profvis({
#   dispersal_n(N = N,
#               v_p_dispersal = c(0.1, 0.1, 0.1),
#               v_theta = c(1, 1, 1),
#               dist_mat = dist_mat)
# })

# v_theta = c(1, 1, 1)
#
# m1_start <- Sys.time()
# m1 = data.matrix(exp(-v_theta[1] * dist_mat))
# diag(m1) <- 0
# m1_end <- Sys.time()
# m1_end - m1_start
#
# m2_start <- Sys.time()
# m2 = exp(-v_theta[1] * dist_mat)
# diag(m2) <- 0
# m2_end <- Sys.time()
# m2_end - m2_start
#
# m3_start <- Sys.time()
# m3 = matrix(exp(-v_theta[1] * dist_mat))
# diag(m3) <- 0
# m3_end <- Sys.time()
# m3_end - m3_start
#
# identical(m1, m2, m3)
# diag(m1)
# diag(m2)


profvis({
  igp_sim(plot_disturbance = FALSE,
          plot_fcl = FALSE,
          plot_patch_dynamics = FALSE,
          disturb_type = "regional",
          disturb_p = 1e-2,
          environment_value = rnorm(20))
})


profvis({
  igp_sim(plot_disturbance = FALSE,
          plot_fcl = FALSE,
          plot_patch_dynamics = FALSE,
          P_pref = NULL)
})

profvis({
  igp_sim(n_patch = 100,
          plot_disturbance = FALSE,
          plot_fcl = FALSE,
          plot_patch_dynamics = FALSE)
})

profvis({
  igp_sim(t = 10000,
          plot_disturbance = FALSE,
          plot_fcl = FALSE,
          plot_patch_dynamics = FALSE)
})

profvis({
  igp_sim(n_patch = 20,
          plot_disturbance = FALSE,
          plot_fcl = FALSE,
          plot_patch_dynamics = FALSE)
})

profvis({
  igp_sim(n_patch = 1000,
          plot_disturbance = FALSE,
          plot_fcl = FALSE,
          plot_patch_dynamics = FALSE)
})

profvis({
  mcsim(n_sp = 3, n_patch = 20)
})

profvis({
  mcsim(n_sp = 3, n_patch = 1000)
})
