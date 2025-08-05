
#' @title Perform the test of Bei (2024) for a given point
#' 
#' @description This function performs the unconditional moment restriction test 
#' as described in Bei (2024).
#' 
#' @param r Result of the projection for which the test should be carried out.
#' @param c The projection matrix. For now, c is restricted to being an
#' elementary vector, i.e. c = (0, ...,0, 1, 0, ..., 0).
#' @param t The time point at which to evaluate theta. 
#' @param par.space Matrix containing 2 columns and d_\theta rows, where d_\theta
#' is the dimension of the parameter space. The first column represents the lower
#' left corner of the parameter space, the second column represents the upper
#' right corner. At least for the time being, only rectangular parameter spaces
#' are allowed.
#' @param data Data frame on which to base the test.
#' @param hp List of hyperparameters needed.
#' @param verbose Boolean variable indicating whether to print updates of the
#' estimation process to the console.
#' @param inst.func.evals Matrix of precomputed instrumental function
#' evaluations for each observation in the data set. Used to speed up the 
#' simulations. If \code{NULL}, the evaluations will be computed during
#' execution of this function. Default is \code{inst.func.evals = NULL}.
#' @param alpha The significance level at which to perform the test. Default is
#' \code{alpha = 0.95}.
#' @param parallel Flag for whether or not parallel computing should be used.
#' Default is \code{parallel = FALSE}.
#' 
#' #' @references Bei, X. (2024). Local linearization based subvector inference in
#' moment inequality models. Journal of Econometrics, 238(1),
#' 105549-. https://doi.org/10.1016/j.jeconom.2023.10554
#' 
test.point_Bei <- function(r, c, t, par.space, data, hp, verbose = FALSE,
                           inst.func.evals = NULL, alpha = 0.95,
                           parallel = FALSE) {
  
  #### 0. Set some useful variables ####
  
  K.bar <- hp[["K.bar"]]
  n.cov <- sum(grepl("X[[:digit:]]+", colnames(data))) - 1
  n.param <- n.cov + 1
  delta.n <- hp[["delta.n"]]
  
  #### 1. Estimate the test statistic ####
  
  if (verbose) {
    message("\t Estimating test statistic...")
  }
  
  # Define initial value
  beta.init <- rowMeans(matrix(par.space[which(c == 0),], nrow = sum(c == 0)))
  
  # Get the test statistic
  out <- get.test.statistic(beta.init, data, par.space, t, hp, c, r,
                            inst.func.evals)
  Jnrh <- out[[1]]
  beta.hat <- out[[2]]
  
  # Define BetaI.r
  BetaI.r <- as.matrix(beta.hat)
  
  #### 2. Calculate the critical value for rn ####
  
  if (verbose) {
    message("\t Computing the initial critical value...")
  }
  
  cvLLn <- get.cvLLn(BetaI.r, data, t, hp, c, r, par.space, inst.func.evals,
                     alpha)
  
  #### 3. Obtain K.bar test statistics, starting from different values ####
  
  if (verbose & (Jnrh > cvLLn)) {
    message("\t Refinement step for critical value is skipped.")
  }
  
  if (Jnrh <= cvLLn) {
    if (verbose) {
      message(sprintf("\t Computing the K.bar (= %s) test statistics...", K.bar))
    }
    
    # Determine a set of K.bar initial values
    initial.values <- matrix(r, nrow = n.param, ncol = K.bar)
    for (idx in which(c == 0)) {
      initial.values[idx,] <- runif(K.bar, min = par.space[idx, 1], max = par.space[idx, 2])
    }
    
    # Initialize an object that will store all test statistic evaluations
    t.stat.evals <- matrix(nrow = K.bar, ncol = n.param + 1)
    colnames(t.stat.evals) <- c(paste0("X", 0:n.cov), "val")
    
    if (parallel) { # Using parallel computing
      
      t.stat.evals[] <-
        foreach(col.idx = 1:ncol(initial.values), .combine = 'rbind') %dopar% {
          
          source("lowLevelFunctions.R")
          
          # Select the parameter corresponding to this iteration
          beta.start <- initial.values[, col.idx]
          
          # Compute the test statistic
          out <- get.test.statistic(beta.start, data, par.space, t, hp, c, r,
                                    inst.func.evals)
          
          # Store the results
          c(out[[2]], out[[1]])
        }
      
    } else { # Using sequential computing
      
      # Compute \tilde{Beta}_I(r)
      for (col.idx in 1:ncol(initial.values)) {
        
        # Select the parameter corresponding to this iteration
        beta.start <- initial.values[, col.idx]
        
        # Compute the test statistic
        out <- get.test.statistic(beta.start, data, par.space, t, hp, c, r,
                                  inst.func.evals)
        
        # Store the results
        t.stat.evals[col.idx, ] <- c(out[[2]], out[[1]])
      }
      
    }
    
    # Reduce it to only the unique parameter vectors
    t.stat.evals <- t.stat.evals[which.unique(t.stat.evals[, 1:n.param]), , drop = FALSE]
    
    # Compute \hat{\Beta}_I(r)
    t.stat.evals <- t.stat.evals[which(t.stat.evals[, n.param + 1] <= Jnrh + delta.n), , drop = FALSE]
  }
  
  #### 4. For each, calculate the critical value ####
  
  if (Jnrh <= cvLLn) {
    
    # If in the previous step no eligible values were found, this step is
    # skipped.
    if (nrow(t.stat.evals) == 0) {
      if (verbose) {
        message("\t No eligible minimizers found in refinement step")
      }
      
    # If eligible values were found, their critical value is computed.
    } else {
      if (verbose) {
        message("\t Recomputing the critical value...")
      }
      
      # Obtain the matrix of covariates in the correct format
      BetaI.r <- t(matrix(t.stat.evals[, 1:n.param], ncol = n.param))
      
      # Compute the critical value
      cvLLn <- get.cvLLn(BetaI.r, data, t, hp, c, r, par.space, inst.func.evals,
                         alpha)
    }
  }
  
  #### 5. Return rh, Jnrh, and cvLLn ####
  
  if (verbose) {
    message("\t Returning the results...")
  }
  
  # Note: 'as.numeric(cvLLn)' gets rid of the name of cvLLn, as it of type
  #       'named num'.
  c("theta" = r, "t.stat" = Jnrh, "crit.val" = as.numeric(cvLLn))
}

#' @title Perform the test of Bei (2024) simultaneously for multiple time
#' points.
#' 
#' @description This function performs the unconditional moment restriction test 
#' as described in Bei (2024). This function directly extends
#' \code{test.point_Bei} by allowing for pairs of moment restrictions over a
#' grid of time points.
#' 
#' @param r Result of the projection for which the test should be carried out.
#' @param c The projection matrix. For now, c is restricted to being an
#' elementary vector, i.e. c = (0, ...,0, 1, 0, ..., 0).
#' @param t The time point at which to evaluate theta. Also allowed to
#' be a vector of time points (used in estimating the model under assumed time-
#' independent coefficients).
#' @param par.space Matrix containing 2 columns and d_\theta rows, where d_\theta
#' is the dimension of the parameter space. The first column represents the lower
#' left corner of the parameter space, the second column represents the upper
#' right corner. At least for the time being, only rectangular parameter spaces
#' are allowed.
#' @param data Data frame on which to base the test.
#' @param hp List of hyperparameters needed.
#' @param verbose Boolean variable indicating whether to print updates of the
#' estimation process to the console.
#' @param inst.func.evals Matrix of precomputed instrumental function
#' evaluations for each observation in the data set. Used to speed up the 
#' simulations. If \code{NULL}, the evaluations will be computed during
#' execution of this function. Default is \code{inst.func.evals = NULL}.
#' @param alpha The significance level at which to perform the test. Default is
#' \code{alpha = 0.95}.
#' @param parallel Flag for whether or not parallel computing should be used.
#' Default is \code{parallel = FALSE}.
#' 
#' #' @references Bei, X. (2024). Local linearization based subvector inference in
#' moment inequality models. Journal of Econometrics, 238(1),
#' 105549-. https://doi.org/10.1016/j.jeconom.2023.10554
#' 
test.point_Bei_MT <- function(r, c, t, par.space, data, hp, verbose = FALSE,
                              inst.func.evals = NULL, alpha = 0.95,
                              parallel = FALSE) {
  
  #### 0. Set some useful variables ####
  
  # Hyperparameters
  K.bar <- hp[["K.bar"]]
  n.cov <- sum(grepl("X[[:digit:]]+", colnames(data))) - 1
  n.param <- n.cov + length(t)
  delta.n <- hp[["delta.n"]]
  
  #### 1. Estimate the test statistic ####
  
  if (verbose) {
    message("\t Estimating test statistic...")
  }
  
  # Define initial value
  beta.init <- rowMeans(matrix(par.space[which(c == 0),], nrow = sum(c == 0)))
  if (length(t) > 1) {
    beta.init[2:length(t)] <- par.space[2:length(t), 1] + 1
  } 
  
  # Get the test statistic
  out <- get.test.statistic(beta.init, data, par.space, t, hp, c, r,
                            inst.func.evals)
  Jnrh <- out[[1]]
  beta.hat <- out[[2]]
  
  # Define BetaI.r
  BetaI.r <- as.matrix(beta.hat)
  
  #### 2. Calculate the critical value for rn ####
  
  if (verbose) {
    message("\t Computing the initial critical value...")
  }
  
  cvLLn <- get.cvLLn(BetaI.r, data, t, hp, c, r, par.space, inst.func.evals,
                     alpha)
  
  #### 3. Obtain K.bar test statistics, starting from different values ####
  
  if (verbose & (Jnrh > cvLLn)) {
    message("\t Refinement step for critical value is skipped.")
  }
  
  if (Jnrh <= cvLLn) {
    if (verbose) {
      message(sprintf("\t Computing the K.bar (= %s) test statistics...", K.bar))
    }
    
    # Determine a set of K.bar initial values
    initial.values <- matrix(r, nrow = n.param, ncol = K.bar)
    for (idx in which(c == 0)) {
      initial.values[idx,] <- runif(K.bar, min = par.space[idx, 1], max = par.space[idx, 2])
    }
    
    # Initialize an object that will store all test statistic evaluations
    t.stat.evals <- matrix(nrow = K.bar, ncol = n.param + 1)
    if (length(t) == 1) {
      colnames(t.stat.evals) <- c(paste0("X", 0:n.cov), "val")
    } else {
      colnames(t.stat.evals) <- c(paste0("X0.", t),
                                  paste0("X", 1:n.cov),
                                  "val")
    }
    
    if (parallel) { # Using parallel computing
      
      t.stat.evals[] <-
        foreach(col.idx = 1:ncol(initial.values), .combine = 'rbind') %dopar% {
          
          source("lowLevelFunctions.R")
          
          # Select the parameter corresponding to this iteration
          beta.start <- initial.values[, col.idx]
          
          # Compute the test statistic
          out <- get.test.statistic(beta.start, data, par.space, t, hp, c, r,
                                    inst.func.evals)
          
          # Store the results
          c(out[[2]], out[[1]])
        }
      
    } else { # Using sequential computing
      
      # Compute \tilde{Beta}_I(r)
      for (col.idx in 1:ncol(initial.values)) {
        
        # Select the parameter corresponding to this iteration
        beta.start <- initial.values[, col.idx]
        
        # Compute the test statistic
        out <- get.test.statistic(beta.start, data, par.space, t, hp, c, r,
                                  inst.func.evals)
        
        # Store the results
        t.stat.evals[col.idx, ] <- c(out[[2]], out[[1]])
      }
      
    }
    
    # Reduce it to only the unique parameter vectors
    t.stat.evals <- t.stat.evals[which.unique(t.stat.evals[, 1:n.param]), , drop = FALSE]
    
    # Compute \hat{\Beta}_I(r)
    t.stat.evals <- t.stat.evals[which(t.stat.evals[, n.param + 1] <= Jnrh + delta.n), , drop = FALSE]
  }
  
  #### 4. For each, calculate the critical value ####
  
  if (Jnrh <= cvLLn) {
    
    # If in the previous step no eligible values were found, this step is
    # skipped.
    if (nrow(t.stat.evals) == 0) {
      if (verbose) {
        message("\t No eligible minimizers found in refinement step")
      }
      
      # If eligible values were found, their critical value is computed.
    } else {
      if (verbose) {
        message("\t Recomputing the critical value...")
      }
      
      # Obtain the matrix of covariates in the correct format
      BetaI.r <- t(matrix(t.stat.evals[, 1:n.param], ncol = n.param))
      
      # Compute the critical value
      cvLLn <- get.cvLLn(BetaI.r, data, t, hp, c, r, par.space, inst.func.evals,
                         alpha)
    }
  }
  
  #### 5. Return rh, Jnrh, and cvLLn ####
  
  if (verbose) {
    message("\t Returning the results...")
  }
  
  # Note: 'as.numeric(cvLLn)' gets rid of the name of cvLLn, as it of type
  #       'named num'.
  c("theta" = r, "t.stat" = Jnrh, "crit.val" = as.numeric(cvLLn))
}
