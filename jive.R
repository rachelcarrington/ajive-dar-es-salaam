jive <- function( dat, joint_rank=1, individual_ranks=rep( 1, length(dat) ), maxiter=1000, tol=10^(-10) ){
  
  # Function to implement the JIVE algorithm

  ## Parameters:
  ## dat: list of data matrices
  ## joint_rank: joint rank (default is 1)
  ## individual_ranks: vector containing individual ranks for each data matrix (default is 1)
  ## maxiter: maximum number of iterations to carry out (default is 1000)
  ## tol: tolerance to use for convergence (default is 10^(-10))

  if ( length(dat) < 2 ){
    stop("At least 2 data matrices must be supplied.")
  }

  # Create concatenated data matrix X
  K <- length(dat) # number of data matrices
  X <- dat[[1]]
  p <- rep(NA, K)
  p[1] <- ncol(dat[[1]])
  for ( i in 2:K ){
    X <- cbind(X, dat[[i]])
    p[i] <- ncol(dat[[i]])
  }
  n <- nrow(X)

  # Initialise joint matrix J and individual matrix A
  X_svd <- svd(X)
  if ( joint_rank > 1 ){
    J_hat <- J_hat_old <- X_svd$u[, 1:joint_rank] %*% diag(X_svd$d[1:joint_rank]) %*% t(X_svd$v[, 1:joint_rank])
  } else {
    J_hat <- J_hat_old <- X_svd$d[1] * X_svd$u[, 1] %*% t(X_svd$v[, 1])
  }
  A_hat <- as.list(rep(NA, K))
  
  # Apply JIVE
  iter <- 0
  dif <- 100
  while ( iter < maxiter ){

    if ( dif > tol ){
      iter <- iter + 1
      
      # Estimate A
      X_star <- X - J_hat
      inds0 <- 1
      for ( i in 1:length(dat) ){
        inds <- inds0 + 0:(p[i] - 1)
        X_star_svd <- svd(X_star[,inds])
        if ( individual_ranks[i] == 1 ){
          A_hat[[i]] <- X_star_svd$d[1] * (X_star_svd$u[,1] %*% t(X_star_svd$v[,1]))
        } else {
          A_hat[[i]] <- X_star_svd$u[, 1:individual_ranks[i]] %*% diag(X_star_svd$d[1:individual_ranks[i]]) %*% t(X_star_svd$v[, 1:individual_ranks[i]])
        }
        inds0 <- inds0 + ncol(dat[[i]])
      }

      # Estimate J
      A <- A_hat[[1]]
      for ( i in 2:K ){
        A <- cbind(A, A_hat[[i]])
      }
      X_star <- X - A
      X_svd <- svd(X)
      if ( joint_rank > 1 ){
        J_hat <- X_svd$u[, 1:joint_rank] %*% diag(X_svd$d[1:joint_rank]) %*% t(X_svd$v[, 1:joint_rank])
      } else {
        J_hat <- X_svd$d[1] * X_svd$u[, 1:joint_rank] %*% t(X_svd$v[, 1:joint_rank])
      }

      # Calculate change
      dif <- sum((J_hat - J_hat_old)^2) / length(dat)
      J_hat_old <- J_hat

    } else {
       
       n.iter <- iter
       iter <- maxiter

    }
  }

  # Save results
  joint_decomp=list(full=J_hat, loadings=svd(J_hat)$v[, 1:joint_rank], d=svd(J_hat)$d[1:joint_rank], 
                     scores=svd(J_hat)$u[, 1:joint_rank])

  individual_decomps <- as.list(rep(NA, K))
  for ( i in 1:K ){
   individual_decomps[[i]] <- list(full=A_hat[[i]], loadings=svd(A_hat[[i]])$v[,1:individual_ranks[i]], 
                                     d=svd(A_hat[[i]])$d[1:individual_ranks[i]], scores=svd(A_hat[[i]])$u[,1:individual_ranks[i]])
  }
  names(individual_decomps) <- names(dat)

  return(list(joint=joint_decomp, individual=individual_decomps, num.iter=n.iter))
}
