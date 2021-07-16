get_fcl <- function(N, lambda_b){
  if(!is.matrix(N)){
    stop("get_fcl() requires N to be a matrix")
  }
  if(ncol(N) != length(lambda_b)){
    stop("get_fcl() requires obs_P_pref to be same length as n_patch in N")
  }
  N[,N[1,] == 0] <- 0
  N[N>0] = 1
  fcl = colSums(N)
  delta  = lambda_b[which(colSums(N)==3)]
  fcl[which(colSums(N)==3)] = delta * 2 + (1 - delta) *3
  fcl
}
