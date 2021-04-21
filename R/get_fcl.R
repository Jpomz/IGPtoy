get_fcl <- function(N){
  if(!is.matrix(N)){
    stop("get_fcl() requires N to be a matrix")
  }
  N[,N[1,] == 0] <- 0
  N[N>0] = 1
  colSums(N)
}
