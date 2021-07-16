get_fcl <- function(N, obs_P_pref){
  if(!is.matrix(N)){
    stop("get_fcl() requires N to be a matrix")
  }
  if(ncol(N) != length(obs_P_pref)){
    stop("get_fcl() requires obs_P_pref to be same length as n_patch in N")
  }
  N[,N[1,] == 0] <- 0
  N[N>0] = 1
  fcl = colSums(N)
  delta  = obs_P_pref[which(colSums(N)==3)]
  fcl[which(colSums(N)==3)] = delta * 2 + (1 - delta) *3
  fcl
}
