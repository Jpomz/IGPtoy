k_n_upstream <- function(k_base = 150,
                         k_c = 10,
                         k_min_exponent = 1.10,
                         k_max_exponent = 1.35,
                         r_max,
                         n_upstream,
                         n_patch){
  # carrying capacity ####
  # power law with exponent between 1.1-1.35 (Koening et al. 2019)
  # k_base = "minimum" carrying capacity
  # k_c = constant multiplier for power law

  # k_exp = exponent for carrying capacity function
  k_exp = runif(n = n_patch,
                min = k_min_exponent,
                max = k_max_exponent)
  k <- round(
    k_base + k_c * n_upstream^k_exp)
  b = (r_max - 1) /k
  list(b = b, k = k)
}
