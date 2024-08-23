# Estimation of joint rank
# In the AJIVE function there are 2 methods for doing this: wedin and random
# In each case, the method returns a threshold for the squared singular values of M (=(U1, ..., Uk)).
# Components are only kept if the square of the corresponding singular value is greater than both thresholds.
# In practice the wedin threshold always seems to be 0, so the random threshold is recommended.

# Random samples
# Essentially we simulate data matrices from N(0,1), calculate M based on these data matrices, and calculate singular values.
# For each random sample, we calculate the (square of the) largest singular value of the estimated M.
# The threshold is the 95% quantile of these values.

get_random_threshold <- function( dat, initial_signal_ranks, nsamples=1000, alpha=0.05 ){
  
  nblocks <- length(dat)
  n <- dim(dat[[1]])[1]
  
  samples <- rep(NA, nsamples)
  for ( iter in 1:nsamples ){
    M_random <- numeric(0)
    for ( i in 1:nblocks ){
      X_random <- matrix(rnorm(n * initial_signal_ranks[i]), nrow=n)
      U_random <- svd(X_random)$u
      M_random <- cbind(M_random, U_random)
    }
    samples[iter] <- svd(M_random)$d[1]^2    
  }
  
  threshold <- quantile(samples, 1 - alpha)
  return(threshold)
}

# Wedin samples
# this always seems to return a threshold of 0, unless alpha is very large (~0.5)

get_wedin_threshold <- function( dat, initial_signal_ranks, nsamples=1000, alpha=0.05 ){
  
  nblocks <- length(dat)
  
  wedin_bound_samples <- numeric(0)
  
  for ( i in 1:nblocks ){
    X <- dat[[i]]
    X_svd <- svd(X)
    m <- min(dim(X))
    r <- initial_signal_ranks[i]
    
    U_resampled_norms <- V_resampled_norms <- rep(NA, nsamples)
    
    # Take U from the SVD of X and remove the first r columns
    U_tilde <- X_svd$u[, (r + 1):m]
    
    for ( index in 1:nsamples ){
      # Resample columns
      col_indexes <- sample(1:dim(U_tilde)[2], dim(U_tilde)[2], replace=TRUE)
      U_tilde_resampled <- U_tilde[, col_indexes]

      # Project onto X
      U_resampled_projection <- t(U_tilde_resampled) %*% X
      
      # Calculate 2-norm of projection
      U_resampled_norms[index] <- norm(U_resampled_projection, type="2")
    }
    
    # Do the same for V
    V_tilde <- X_svd$v[, (r + 1):m]
    
    for ( index in 1:nsamples ){
      # Resample columns
      col_indexes <- sample(1:dim(V_tilde)[2], dim(V_tilde)[2], replace=TRUE)
      V_tilde_resampled <- V_tilde[, col_indexes]
      
      # Project onto X
      V_resampled_projection <- X %*% V_tilde_resampled
      
      # Calculate 2-norm of projection
      V_resampled_norms[index] <- norm(V_resampled_projection, type="2")
    }
    
    sigma_min <- X_svd$d[r] # largest singular value in signal
    
    wedin_bound_samples_block <- rep(NA, nsamples)
    for ( j in 1:nsamples ){
      wedin_bound_samples_block[j] <- min((max(U_resampled_norms[j], V_resampled_norms[j]) / sigma_min)^2, 1)
    }
    
    wedin_bound_samples <- rbind(wedin_bound_samples, wedin_bound_samples_block)
    
  }
  
  wedin_samples <- nblocks - colSums(wedin_bound_samples)
  wedin_threshold <- quantile(wedin_samples, alpha)
  
  return(wedin_threshold)

}
