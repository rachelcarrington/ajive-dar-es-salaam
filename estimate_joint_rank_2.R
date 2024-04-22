
estimate_joint_rank <- function(dat, initial_ranks, nsamples=1000){
  nblocks <- length(dat)
  n <- dim(dat[[1]])[1]
  initial_ranks_totals <- c(0, cumsum(initial_ranks))
  for ( iter in 1:nsamples ){
    U_iter <- matrix(NA, nrow=n, ncol=sum(initial_ranks))
    for ( k in 1:nblocks ){
      Z <- matrix(rnorm(n * initial_ranks[k]), nrow=n)
      U_iter[(initial_ranks_totals[k] + 1):initial_ranks_totals[k + 1], ] <- cbind(U_iter, svd(Z)$u)
    }
    vals[iter] <- svd(U_iter)$d[1]
  }
  threshold <- quantile(vals, 0.95)

  U <- numeric(0)
  for ( i in 1:nblocks ){
    U <- cbind(U, svd(dat[[1]])$u[,1:initial_ranks[i]])
  }
  joint_rank_est <- sum(svd(U)$d > threshold)
  
  return(joint_rank_est)
}
