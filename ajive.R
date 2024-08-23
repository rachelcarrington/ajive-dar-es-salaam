# I think this function should also allow for ranks to be estimated
# for initial ranks we use a scree plot, these should be supplied; but for joint and individual ranks they can be estimated within the function

ajive <- function( dat, joint_rank=NULL, individual_ranks=NULL, initial_ranks=NULL ){

  # Function to implement AJIVE
  # If joint_rank = 0, implements PCA on each data matrix individually
  # Joint and individual ranks default to 1, if not specified

  # dat: list of data blocks: rows correspond to individuals, columns to variables
  # Each data block should have the same number of rows
  # joint_rank: should be a positive integer, defaults to 1 if not set
  # individual_ranks: vector of positive integers, of length equal to number of blocks, ranks default to 1 if not set
  # initial_ranks: initial ranks for first step of the algorithm, length equal to number of blocks,
  # defaults to sum of joint and individual ranks

  nblocks <- length(dat)

  # Check there are at least 2 blocks, with the same number of rows in each.
  if ( nblocks < 2 ){
    stop("There must be at least 2 data matrices to implement AJIVE.")
  } else {
    nrows_block1 <- dim(dat[[1]])[1]
    for ( i in 2:nblocks ){
      if ( dim(dat[[i]])[1] != nrows_block1 ){
        stop("Data matrices must have the same number of rows.")
      }
    }
  }

  # Set joint and individual ranks to 1, if not specified.
  if( is.null(joint_rank) ){
    joint_rank <- 1
  } else if ( joint_rank < 0 ){
    stop("The joint rank should be a non-negative integer.")
  } else if ( joint_rank != round(joint_rank) ){
    stop("The joint rank should be a non-negative integer.")
  }

  if( is.null(individual_ranks) ){
    individual_ranks <- rep(1, nblocks)
  } else if ( sum(individual_ranks < 0) >= 1 ){
    stop("The individual ranks should be non-negative integers.")
  } else if ( sum(abs(individual_ranks - round(individual_ranks))) != 0 ){
    stop("The individual ranks should be non-negative integers.")
  }


  if ( joint_rank == 0 ){
    # do PCA on individual matrices

    joint_decomp <- NA
    individual_decomps <- as.list(rep(NA, length(individual_ranks)))
    for ( i in 1:length(individual_ranks) ){
      if ( individual_ranks[i] != 0 ){
        svd_dat <- svd(dat[[i]])
        U <- svd_dat$u[, 1:individual_ranks[i], drop=FALSE]
        V <- svd_dat$v[, 1:individual_ranks[i], drop=FALSE]
        Sigma <- ifelse(individual_ranks[i] == 1, svd_dat$d[1], diag(svd_dat$d[1:individual_ranks[i]]))
        individual_decomps[[i]] <- list(full=U %*% Sigma %*% t(V), loadings=U, d=Sigma, scores=V)
      }
    }

  } else {

    # Set initial ranks for algorithm to sum of joint and individual ranks, unless otherwise specified.
    if ( is.null(initial_ranks) ){
      initial_ranks <- individual_ranks + joint_rank
      for ( i in 1:length(initial_ranks) ){
        initial_ranks[i] <- min(initial_ranks[i], min(dim(dat[[i]])))
      }
    }

    # Create concatenated data matrix X
    X <- dat[[1]]
    for ( i in 2:nblocks ){
      X <- cbind(X, dat[[i]])
    }
    n <- dim(X)[1]

    ############ AJIVE #####################################################################################################

    # Stage 1. Initial SVD, create M.
    M <- numeric(0)
    for ( i in 1:nblocks ){
      M <- cbind(M, svd(dat[[i]])$u[,1:initial_ranks[i]])
    }

    # 2. Find U_J.
    UJ <- svd(M)$u[,1:joint_rank]

    # 3a. Find J.
    J_hat <- UJ %*% t(UJ) %*% X

    svd_J <- svd(J_hat)

    if ( joint_rank == 1 ){
      joint_decomp=list(full=J_hat, loadings=matrix(svd_J$v[,1], ncol=1), d=svd_J$d[1], scores=matrix(svd_J$u[,1], ncol=1))
    } else {
      joint_decomp=list(full=J_hat, loadings=svd_J$v[,1:joint_rank], d=svd_J$d[1:joint_rank], scores=svd_J$u[,1:joint_rank] )
    }
    rownames(joint_decomp$full) <- rownames(joint_decomp$scores) <- rownames(dat[[1]])

    # 3b. Find A.
    individual_decomps <- as.list(rep(0, nblocks))
    for ( i in 1:nblocks ){
      A_i_tilde <- (diag(n) - UJ %*% t(UJ)) %*% dat[[i]]
      svd_A_i_tilde <- svd(A_i_tilde)
      if ( individual_ranks[i] == 1 ){
        A_i <- svd_A_i_tilde$d[1] * (svd_A_i_tilde$u[,1] %*% t(svd_A_i_tilde$v[,1]))
        svd_A_i <- svd(A_i)
        individual_decomps[[i]] <- list(full=A_i, loadings=matrix(svd_A_i$v[,1], ncol=1), d=svd_A_i$d[1], scores=matrix(svd_A_i$u[,1], ncol=1))
        rownames(individual_decomps[[i]]$full) <- rownames(individual_decomps[[i]]$scores) <- rownames(dat[[i]])
        colnames(individual_decomps[[i]]$full) <- rownames(individual_decomps[[i]]$loadings) <- colnames(dat[[i]])
      } else if ( individual_ranks[i] > 1 ){
        A_i <- svd_A_i_tilde$u[,1:individual_ranks[i]] %*% diag(svd_A_i_tilde$d[1:individual_ranks[i]]) %*% t(svd_A_i_tilde$v[,1:individual_ranks[i]])
        svd_A_i <- svd(A_i)
        individual_decomps[[i]] <- list(full=A_i, loadings=svd_A_i$v[,1:individual_ranks[i]], d=svd_A_i$d[1:individual_ranks[i]], scores=svd_A_i$u[,1:individual_ranks[i]])
        rownames(individual_decomps[[i]]$full) <- rownames(individual_decomps[[i]]$scores) <- rownames(dat[[i]])
        colnames(individual_decomps[[i]]$full) <- rownames(individual_decomps[[i]]$loadings) <- colnames(dat[[i]])
      } else {
        individual_decomps[[i]] <- list(full=NA, loadings=NA, d=NA, scores=NA)
      }

    }

  }

  names(individual_decomps) <- names(dat)

  return(list(joint=joint_decomp, individual=individual_decomps))
}
