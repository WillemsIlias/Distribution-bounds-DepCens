
#### Import dependencies ####

require("EnvStats")
require("splines2")
require("copula")
require("nloptr")

#### General utility functions ####

#' @title Tests whether two parameter vectors are approximately equal
#' 
#' @description
#' Function to check whether the given parameters are approximately equal, where
#' approximate equality is defined as the norm of the differences being smaller
#' than the given tolerance.
#' 
#' @param par1 First parameter
#' @param par2 Second parameter
#' @param tol Tolerance. Differences between \code{par1} and \code{par2} smaller
#' than \code{tol} will be neglected. Default value is \code{tol = 10^(-4)}.
#' 
#' @noRd
#' 
par.equal <- function(par1, par2, tol = 10^(-4)) {
  distance <- sqrt(sum((par1 - par2)^2))
  distance < tol
}

#' @title Returns the indices of the (approximately) unique parameter vectors 
#' 
#' @description
#' Returns the indices of the approximately unique parameter vectors in a 
#' matrix, where unicity is defined as there not existing another parameter 
#' vector that is equal up until a difference in norm less than the given
#' tolerance.
#' 
#' @param par.mat Matrix of parameter vectors.
#' @param tol Tolerance. Default value is \code{tol = 10^(-4)}.
#' 
#' @noRd
#' 
which.unique <- function(par.mat, tol = 10^(-4)) {
  
  # If par.mat is a vector, return 1
  if (!is.matrix(par.mat)) {
    return(c(1))
  }
  
  # Initialize a vector that will store all unique column indices
  idx.unique <- c()
  
  # For all but the last row, check whether there are duplicate rows that
  # follow.
  if (nrow(par.mat) > 1) {
    for (i in 1:(nrow(par.mat) - 1)) {
      dupe <- FALSE
      for (j in (i+1):nrow(par.mat)) {
        if (par.equal(par.mat[i,], par.mat[j,], tol = tol)) {
          dupe <- TRUE
        }
      }
      
      if (!dupe) {
        idx.unique <- c(idx.unique, i)
      }
    }
  }
  
  # The last element is always unique
  idx.unique <- c(idx.unique, nrow(par.mat))
  
  # Return the result
  idx.unique
}

#' @title Checks if a directory exists and, if necessary, creates it.
#' 
#' @description
#' Function that checks whether a given directory exists and, if necessary,
#' creates it.
#' 
#' @param dir.name Name of the directory whose existence is to be checked.
#' @param path Path to the parent directory in which the given directory should
#' exist. Default value is \code{path = NULL}, in which case the working
#' directory is used.
#' 
#' @noRd
#' 
check_create.dir <- function(dir.name, path = NULL) {
  
  # Path to directory
  if (!is.null(path)) {
    dir.path <- paste(c(path, dir.name), collapse = "/")
  } else {
    dir.path <- dir.name
  }
  
  # If the directory does not exist yet, create it
  if (!dir.exists(dir.path)) {
    dir.create(dir.path)
  }
}

#### Data generation function ####

#' @title Compute the average censoring in a simulation design.
#' 
#' @description This function computes the average censoring in a simulation
#' design as defined in generateData_simMain.R. This function assumes that
#' parameters 'beta.true', 'n', 'n.cov' and 'options' are already devined in the
#' global environment.
#' 
#' @param DGP The index of the DGP for for which to compute the average
#' censoring.
#' @param DGPset {Value indicating which DGP function should be called. Default
#' is \code{DGPset = "Main"}.
#'
#' @note
#' This function is deprecated and should not be used.
#' 
#' @noRd
#'
get.average.censoring <- function(DGP, DGPset = "Main") {
  
  # NOTE: DEPRECATED
  warning("The function 'get.average.censoring' is deprecated!")
  
  # Make sure the code in the environment is up-to-date
  source("lowLevelFunctions.R")
  
  # Generate K data sets and compute the censoring percentage in each
  K <- 250
  cens.per.vec <- rep(0, K)
  for (i in 1:K) {
    if (DGPset == "Main") {
      
      # Set correct DGP
      options.copy <- options
      options.copy[["DGP"]] <- DGP
      
      # Generate data
      data <- generateData_simMain(beta.true, n, n.cov, options.copy, H0.inv,
                                   plot.data = i == 1)
    } else if (DGPset == "Miss") {
      
      # Set correct DGP
      options.data.gen.copy <- options.data.gen
      options.data.gen.copy[["DGP"]] <- DGP
      
      # Generate data
      data <- generateData_simMiss(beta.true, n, n.cov, options.data.gen.copy, H0.inv,
                                   plot.data = i == 1)
    }
    
    cens.per.vec[i] <- 1 - sum(data$Delta)/nrow(data)
  }
  
  # Return the result
  mean(cens.per.vec)
}

#' @title Test whether covariates lie inside {1}^d.d x [0.5, 1.5]^d.c
#' 
#' @description This function tests for each covariate vector in the given
#' matrix of covariates whether or not it lies inside the region {1}^d.d x
#' [0.5, 1.5]^d.c. This function is used in generating data according to some
#' DGPs in the function 'generateData.R'.
#' 
#' @param X Either a matrix containing the covariates (and intercept), or a data
#' frame containing the covariates, named X1, X2, etc.
#' @param type.cov Vector containing the type of covariates.
#' 
#' @noRd
#'
inRegion1 <- function (X, type.cov) {
  if (is.data.frame(X)) {
    X <- X[, grepl("X[[:digit:]]*", colnames(X))]
  }
  check.in.region <- function(row) {
    all(0.5 <= row[type.cov == "c"]) & all(row[type.cov == "c"] <= 1.5) &
      all(row[type.cov == "d"] == 1)
  }
  apply(X[, -1, drop = FALSE], 1, check.in.region)
}

#' @title Test whether covariates lie inside {1}^d.d x [-Inf, -1]^d
#' 
#' @description This function tests for each covariate vector in the given
#' matrix of covariates whether or not it lies inside the region {1}^d.d x 
#' [-Inf, 1]^d. This function is used in generating data according to some DGPs
#' in the function 'generateData.R'.
#' 
#' @param X Either a matrix containing the covariates (and intercept), or a data
#' frame containing the covariates, named X1, X2, etc.
#' @param type.cov Vector containing the type of covariates.
#' 
#' @noRd
#'
inRegion2 <- function (X, type.cov) {
  if (is.data.frame(X)) {
    X <- X[, grepl("X[[:digit:]]*", colnames(X))]
  }
  check.in.region <- function(row) {
    all(row[type.cov == "c"] <= -1) & all(row[type.cov == "d"] == 1)
  }
  apply(X[, -1, drop = FALSE], 1, check.in.region)
}

#' @title Test whether covariates lie inside {0}^d.d x ([-1, 0] x [0, 1])^(d/2)
#' 
#' @description This function tests for each covariate vector in the given
#' matrix of covariates whether or not it lies inside the region {0}^d.d x 
#' ([-1, 0] x [0, 1])^(d/2). This function is used in generating data according
#' to some DGPs in the function 'generateData.R'. This function is also used in 
#' the implementation of \code{'inRegion4.R'}.
#' 
#' @param X Either a matrix containing the covariates (and intercept), or a data
#' frame containing the covariates, named X1, X2, etc.
#' @param type.cov Vector containing the type of covariates.
#' 
#' @noRd
#'
inRegion3 <- function (X, type.cov) {
  if (is.data.frame(X)) {
    X <- X[, grepl("X[[:digit:]]*", colnames(X))]
  }
  check.in.region <- function(row) {
    lbs <- rep(c(-1, 0), sum(type.cov == "c"))[1:sum(type.cov == "c")]
    ubs <- rep(c(0, 1), sum(type.cov == "c"))[1:sum(type.cov == "c")]
    equal.d <- all(row[type.cov != "c"] == 1)
    equal.c <- all(lbs <= row[type.cov == "c"]) & all(row[type.cov == "c"] <= ubs)
    equal.d & equal.c
  }
  apply(X[, -1, drop = FALSE], 1, check.in.region)
}

#' @title Test whether covariates lie inside {0}^d.d x ([0, 1] x [-1, 0])^(d/2)
#' 
#' @description This function tests for each covariate vector in the given
#' matrix of covariates whether or not it lies inside the region {0}^d.d x 
#' ([0, 1] x [-1, 0])^(d/2). This function is used in generating data according
#' to some DGPs in the function 'generateData.R'.
#' 
#' @param X Either a matrix containing the covariates (and intercept), or a data
#' frame containing the covariates, named X1, X2, etc.
#' @param type.cov Vector containing the type of covariates.
#' 
#' @noRd
#'
inRegion4 <- function (X, type.cov) {
  if (is.data.frame(X)) {
    X <- X[, grepl("X[[:digit:]]*", colnames(X))]
  }
  inRegion3(-X, type.cov)
}

#' @title Test whether covariates lie inside {0}^d.d x ([-0.5, 0.5] x
#' [-2, 2])^(d/2)
#' 
#' @description This function tests for each covariate vector in the given
#' matrix of covariates whether or not it lies inside the region {0}^d.d x 
#' ([-0.5, 0.5] x [-2, 2])^(d/2). This function is used in generating data
#' according to some DGPs in the function 'generateData_add.R'.
#' 
#' @param X Either a matrix containing the covariates (and intercept), or a data
#' frame containing the covariates, named X1, X2, etc.
#' @param type.cov Vector containing the type of covariates.
#' 
#' @noRd
#'
inRegion5 <- function (X, type.cov) {
  if (is.data.frame(X)) {
    X <- X[, grepl("X[[:digit:]]*", colnames(X))]
  }
  check.in.region <- function(row) {
    lbs <- rep(c(-0.5, -2), sum(type.cov == "c"))[1:sum(type.cov == "c")]
    ubs <- rep(c(0.5, 2), sum(type.cov == "c"))[1:sum(type.cov == "c")]
    equal.d <- all(row[type.cov != "c"] == 1)
    equal.c <- all(lbs <= row[type.cov == "c"]) & all(row[type.cov == "c"] <= ubs)
    equal.d & equal.c
  }
  apply(X[, -1, drop = FALSE], 1, check.in.region)
}

#' @title Test whether covariates lie inside {0}^d.d x ([-2, 2] x
#' [-0.5, 0.5])^(d/2)
#' 
#' @description This function tests for each covariate vector in the given
#' matrix of covariates whether or not it lies inside the region {0}^d.d x 
#' ([-2, 2] x [-0.5, 0.5])^(d/2). This function is used in generating data
#' according to some DGPs in the function 'generateData_add.R'.
#' 
#' @param X Either a matrix containing the covariates (and intercept), or a data
#' frame containing the covariates, named X1, X2, etc.
#' @param type.cov Vector containing the type of covariates.
#' 
#' @noRd
#'
inRegion6 <- function (X, type.cov) {
  if (is.data.frame(X)) {
    X <- X[, grepl("X[[:digit:]]*", colnames(X))]
  }
  check.in.region <- function(row) {
    lbs <- rep(c(-2, -0.5), sum(type.cov == "c"))[1:sum(type.cov == "c")]
    ubs <- rep(c(2, 0.5), sum(type.cov == "c"))[1:sum(type.cov == "c")]
    equal.d <- all(row[type.cov != "c"] == 1)
    equal.c <- all(lbs <= row[type.cov == "c"]) & all(row[type.cov == "c"] <= ubs)
    equal.d & equal.c
  }
  apply(X[, -1, drop = FALSE], 1, check.in.region)
}

#' @title Generates a data set according to the specified arguments.
#' 
#' @description This function generates a data set according to the specified
#' arguments.
#' 
#' @param beta.true True covariate vector, as a function of time.
#' @param n Sample size.
#' @param n.cov Number of covariates.
#' @param options List of additional arguments.
#' @param plot.data Boolean value indicating whether or not to plot the
#' generated data set. Default value is \code{plot.data = FALSE}.
#' 
#' @import stats copula
#' 
#' @noRd
#' 
generateData <- function(beta.true, n, n.cov, options, plot.data = FALSE) {
  
  # Extract the necessary hyperparameters
  if (options[["link.function"]] == "AFT_ll") {
    inv.Lambda <- Lambda_inverse_AFT_ll
  } else if (options[["link.function"]] == "Cox_wb") {
    inv.Lambda <- Lambda_inverse_Cox_wb
  }
  DGP <- options[["DGP"]]
  
  # Set the types of the covariates to be used based on the selected DGP.
  #
  # NOTE TO SELF: If you change these thresholds (20, 40, 60), also change them
  #               accordingly in 'simulate1D.R', 'set.hyperparameters.R', 
  #               'simFuncWrapper.R' and 'simulate1D.CCDF.R'
  if (DGP <= 20) {
    type.cov <- rep("c", n.cov)
  } else if (20 < DGP & DGP <= 40) {
    type.cov <- rep("c", n.cov)
    if ((ceiling(n.cov/2) + 1) <= length(type.cov)) {
      type.cov[(ceiling(n.cov/2) + 1):length(type.cov)] <- "b"
    }
  } else {
    type.cov <- rep("b", n.cov)
  }
  
  # For the data generation, we can make the following shortcut.
  beta.true <- beta.true(0)
  
  # Subset beta.true to the correct number of parameters
  beta.true <- beta.true[1:(n.cov + 1)]
  
  # Generate the intercept and the covariates
  X <- rep(1, n)
  for (i in 1:n.cov) {
    if (type.cov[i] == "c") {
      X <- cbind(X, rnorm(n))
    } else if (type.cov[i] == "b") {
      X <- cbind(X, rbinom(n, 1, 0.5))
    } else {
      stop("Invalid value for type of covariate")
    }
  }
  colnames(X) <- paste0("X", 0:n.cov)
  
  # Determine which covariate vectors lie in the predefined regions.
  idxs.region1 <- which(inRegion1(X, type.cov))
  idxs.region2 <- which(inRegion2(X, type.cov))
  idxs.region3 <- which(inRegion3(X, type.cov))
  idxs.region4 <- which(inRegion4(X, type.cov))
  
  # Set dependence structure based on selected DGP
  if (DGP %% 20 == 1) { # Independence, C ~ Unif
    
    # Latent event time T
    T <- inv.Lambda(runif(n)) - X %*% beta.true
    
    # Latent censoring time C
    C <- pmax(runif(n, min(T), max(T)), runif(n, min(T), max(T)))
    
  } else if (DGP %% 20 == 2) { # Positive dependence, C ~ Unif
    
    # Specify copula
    cop_to_use <- frankCopula(param = 6, dim = 2)
    
    # Generate data from copula
    U <- rCopula(n, cop_to_use)
    u1 <- U[,1]
    u2 <- U[,2]
    
    # Latent event time T
    T <- inv.Lambda(u1) - X %*% beta.true
    
    # Define quantile function of C
    FC_given_X_inv <- function(u2) {u2 * (max(T) - min(T)) + min(T) + 2}
    C <- FC_given_X_inv(u2)
    
  } else if (DGP %% 20 == 3) { # Negative dependence, C ~ Unif
    
    # Specify copula
    cop_to_use <- frankCopula(param = -6, dim = 2)
    
    # Generate data from copula
    U <- rCopula(n, cop_to_use)
    u1 <- U[,1]
    u2 <- U[,2]
    
    # Latent event time T
    T <- inv.Lambda(u1) - X %*% beta.true
    
    # Define quantile function of C
    FC_given_X_inv <- function(u2) {u2 * (max(T) - min(T)) + min(T) + 4}
    C <- FC_given_X_inv(u2)
    
  } else if (DGP %% 20 == 4) { # Independence, C ~ AFT_ll
    
    # Latent event time T
    T <- inv.Lambda(runif(n)) - X %*% beta.true
    
    # Latent censoring time C
    gamma <- beta.true
    gamma[1] <- gamma[1] - 1.3
    C <- Lambda_inverse_AFT_ll(runif(n)) - X %*% gamma
    
  } else if (DGP %% 20 == 5) { # Positive dependence, C ~ AFT_ll
    
    # Specify copula
    cop_to_use <- frankCopula(param = 6, dim = 2)
    
    # Generate data from copula
    U <- rCopula(n, cop_to_use)
    u1 <- U[,1]
    u2 <- U[,2]
    
    # Latent event time T
    T <- inv.Lambda(u1) - X %*% beta.true
    
    # Latent censoring time C
    gamma <- beta.true
    gamma[1] <- gamma[1] - 1.3
    C <- Lambda_inverse_AFT_ll(u2) - X %*% gamma
    
  } else if (DGP %% 20 == 6) { # Negative dependence, C ~ AFT_ll
    
    # Specify copula
    cop_to_use <- frankCopula(param = -6, dim = 2)
    
    # Generate data from copula
    U <- rCopula(n, cop_to_use)
    u1 <- U[,1]
    u2 <- U[,2]
    
    # Latent event time T
    T <- inv.Lambda(u1) - X %*% beta.true
    
    # Latent censoring time C
    gamma <- beta.true
    gamma[1] <- gamma[1] - 1.3
    C <- Lambda_inverse_AFT_ll(u2) - X %*% gamma
    
  } else if (DGP %% 20 == 7) { # Independence, C ~ AFT_ll, high cens
    
    # Latent event time T
    T <- inv.Lambda(runif(n)) - X %*% beta.true
    
    # Latent censoring time C
    gamma <- beta.true
    C <- Lambda_inverse_AFT_ll(runif(n)) - X %*% gamma
    
  } else if (DGP %% 20 == 8) { # Positive dependence, C ~ AFT_ll, high cens
    
    # Specify copula
    cop_to_use <- frankCopula(param = 6, dim = 2)
    
    # Generate data from copula
    U <- rCopula(n, cop_to_use)
    u1 <- U[,1]
    u2 <- U[,2]
    
    # Latent event time T
    T <- inv.Lambda(u1) - X %*% beta.true
    
    # Latent censoring time C
    gamma <- beta.true
    C <- Lambda_inverse_AFT_ll(u2) - X %*% gamma
    
  } else if (DGP %% 20 == 9) { # Negative dependence, C ~ AFT_ll, high cens
    
    # Specify copula
    cop_to_use <- frankCopula(param = -6, dim = 2)
    
    # Generate data from copula
    U <- rCopula(n, cop_to_use)
    u1 <- U[,1]
    u2 <- U[,2]
    
    # Latent event time T
    T <- inv.Lambda(u1) - X %*% beta.true
    
    # Latent censoring time C
    gamma <- beta.true
    C <- Lambda_inverse_AFT_ll(u2) - X %*% gamma
    
  } else if (DGP %% 20 == 10) { # Independence, C ~ AFT_ll, high cens, regions with less cens
    
    # Latent event time T
    T <- inv.Lambda(runif(n)) - X %*% beta.true
    
    # Latent censoring time C
    gamma <- beta.true
    gamma.region1 <- beta.true
    gamma.region2 <- beta.true
    gamma.region1[1] <- beta.true[1] - 3
    gamma.region2[1] <- beta.true[1] - 2
    
    idxs.outside <- which(!inRegion1(X, type.cov) & !inRegion2(X, type.cov))
    C <- rep(0, n)
    C[idxs.outside] <- Lambda_inverse_AFT_ll(runif(length(idxs.outside))) - X[idxs.outside, ] %*% gamma
    C[idxs.region1] <- Lambda_inverse_AFT_ll(runif(length(idxs.region1))) - X[idxs.region1, ] %*% gamma.region1
    C[idxs.region2] <- Lambda_inverse_AFT_ll(runif(length(idxs.region2))) - X[idxs.region2, ] %*% gamma.region2
  
  } else if (DGP %% 20 == 11) { # Pos. dep., C ~ AFT_ll, high cens, regions with less cens
      
    # Specify copula
    cop_to_use <- frankCopula(param = 6, dim = 2)
    
    # Generate data from copula
    U <- rCopula(n, cop_to_use)
    u1 <- U[,1]
    u2 <- U[,2]
    
    # Latent event time T
    T <- inv.Lambda(u1) - X %*% beta.true
    
    # Latent censoring time C
    gamma <- beta.true
    gamma.region1 <- beta.true
    gamma.region2 <- beta.true
    gamma.region1[1] <- beta.true[1] - 3
    gamma.region2[1] <- beta.true[1] - 2
    
    idxs.outside <- which(!inRegion1(X, type.cov) & !inRegion2(X, type.cov))
    C <- rep(0, n)
    C[idxs.outside] <- Lambda_inverse_AFT_ll(u2[idxs.outside]) - X[idxs.outside, ] %*% gamma
    C[idxs.region1] <- Lambda_inverse_AFT_ll(u2[idxs.region1]) - X[idxs.region1, ] %*% gamma.region1
    C[idxs.region2] <- Lambda_inverse_AFT_ll(u2[idxs.region2]) - X[idxs.region2, ] %*% gamma.region2
  
  } else if (DGP %% 20 == 12) { # Neg. dep., C ~ AFT_ll, high cens, regions with less cens
    
    # Specify copula
    cop_to_use <- frankCopula(param = -6, dim = 2)
    
    # Generate data from copula
    U <- rCopula(n, cop_to_use)
    u1 <- U[,1]
    u2 <- U[,2]
    
    # Latent event time T
    T <- inv.Lambda(u1) - X %*% beta.true
    
    # Latent censoring time C
    gamma <- beta.true
    gamma.region1 <- beta.true
    gamma.region2 <- beta.true
    gamma.region1[1] <- beta.true[1] - 3
    gamma.region2[1] <- beta.true[1] - 2
    
    idxs.outside <- which(!inRegion1(X, type.cov) & !inRegion2(X, type.cov))
    C <- rep(0, n)
    C[idxs.outside] <- Lambda_inverse_AFT_ll(u2[idxs.outside]) - X[idxs.outside, ] %*% gamma
    C[idxs.region1] <- Lambda_inverse_AFT_ll(u2[idxs.region1]) - X[idxs.region1, ] %*% gamma.region1
    C[idxs.region2] <- Lambda_inverse_AFT_ll(u2[idxs.region2]) - X[idxs.region2, ] %*% gamma.region2
    
  } else if (DGP %% 20 == 13) { # Positive dependence, cens comparable to DGP = 11
    
    # Specify copula
    cop_to_use <- frankCopula(param = 6, dim = 2)
    
    # Generate data from copula
    U <- rCopula(n, cop_to_use)
    u1 <- U[,1]
    u2 <- U[,2]
    
    # Latent event time T
    T <- inv.Lambda(u1) - X %*% beta.true
    
    # Latent censoring time C
    gamma <- beta.true
    gamma[1] <- beta.true[1] - 0.115
    C <- Lambda_inverse_AFT_ll(u2) - X %*% gamma
    
  } else if (DGP %% 20 == 14) { # Negative dependence, cens comparable to DGP = 12
    
    # Specify copula
    cop_to_use <- frankCopula(param = -6, dim = 2)
    
    # Generate data from copula
    U <- rCopula(n, cop_to_use)
    u1 <- U[,1]
    u2 <- U[,2]
    
    # Latent event time T
    T <- inv.Lambda(u1) - X %*% beta.true
    
    # Latent censoring time C
    gamma <- beta.true
    gamma[1] <- beta.true[1] - 0.2
    C <- Lambda_inverse_AFT_ll(u2) - X %*% gamma
    
  } else if (DGP %% 20 == 15) { # Independence, C ~ AFT_ll, high cens, regions with less cens 2
    
    # Latent event time T
    T <- inv.Lambda(runif(n)) - X %*% beta.true
    
    # Latent censoring time C
    gamma <- beta.true
    gamma.region1 <- beta.true
    gamma.region4 <- beta.true
    gamma[1] <- beta.true[1] + 0.43
    gamma.region1[1] <- beta.true[1] - 3
    gamma.region4[1] <- beta.true[1] - 2
    
    idxs.outside <- which(!inRegion1(X, type.cov) & !inRegion4(X, type.cov))
    C <- rep(0, n)
    C[idxs.outside] <- Lambda_inverse_AFT_ll(runif(length(idxs.outside))) - X[idxs.outside, ] %*% gamma
    C[idxs.region1] <- Lambda_inverse_AFT_ll(runif(length(idxs.region1))) - X[idxs.region1, ] %*% gamma.region1
    C[idxs.region4] <- Lambda_inverse_AFT_ll(runif(length(idxs.region4))) - X[idxs.region4, ] %*% gamma.region4
    
  } else if (DGP %% 20 == 16) { # Pos. dep., C ~ AFT_ll, high cens, regions with less cens 2
    
    # Specify copula
    cop_to_use <- frankCopula(param = 6, dim = 2)
    
    # Generate data from copula
    U <- rCopula(n, cop_to_use)
    u1 <- U[,1]
    u2 <- U[,2]
    
    # Latent event time T
    T <- inv.Lambda(u1) - X %*% beta.true
    
    # Latent censoring time C
    gamma <- beta.true
    gamma.region1 <- beta.true
    gamma.region4 <- beta.true
    gamma[1] <- beta.true[1] + 0.15
    gamma.region1[1] <- beta.true[1] - 3
    gamma.region4[1] <- beta.true[1] - 2
    
    idxs.outside <- which(!inRegion1(X, type.cov) & !inRegion4(X, type.cov))
    C <- rep(0, n)
    C[idxs.outside] <- Lambda_inverse_AFT_ll(u2[idxs.outside]) - X[idxs.outside, ] %*% gamma
    C[idxs.region1] <- Lambda_inverse_AFT_ll(u2[idxs.region1]) - X[idxs.region1, ] %*% gamma.region1
    C[idxs.region4] <- Lambda_inverse_AFT_ll(u2[idxs.region4]) - X[idxs.region4, ] %*% gamma.region4
    
  } else if (DGP %% 20 == 17) { # Neg. dep., C ~ AFT_ll, high cens, regions with less cens 2
    
    # Specify copula
    cop_to_use <- frankCopula(param = -6, dim = 2)
    
    # Generate data from copula
    U <- rCopula(n, cop_to_use)
    u1 <- U[,1]
    u2 <- U[,2]
    
    # Latent event time T
    T <- inv.Lambda(u1) - X %*% beta.true
    
    # Latent censoring time C
    gamma <- beta.true
    gamma.region1 <- beta.true
    gamma.region4 <- beta.true
    gamma[1] <- beta.true[1] + 0.2
    gamma.region1[1] <- beta.true[1] - 3
    gamma.region4[1] <- beta.true[1] - 2
    
    idxs.outside <- which(!inRegion1(X, type.cov) & !inRegion4(X, type.cov))
    C <- rep(0, n)
    C[idxs.outside] <- Lambda_inverse_AFT_ll(u2[idxs.outside]) - X[idxs.outside, ] %*% gamma
    C[idxs.region1] <- Lambda_inverse_AFT_ll(u2[idxs.region1]) - X[idxs.region1, ] %*% gamma.region1
    C[idxs.region4] <- Lambda_inverse_AFT_ll(u2[idxs.region4]) - X[idxs.region4, ] %*% gamma.region4
    
  } else if (DGP %% 20 == 18) { # Independence, C ~ AFT_ll, high cens, regions with less cens 3
    
    # Latent event time T
    T <- inv.Lambda(runif(n)) - X %*% beta.true
    
    # Latent censoring time C
    gamma <- beta.true
    gamma.region3 <- beta.true
    gamma.region4 <- beta.true
    gamma[1] <- beta.true[1] + 0.58
    gamma.region3[1] <- beta.true[1] - 3
    gamma.region4[1] <- beta.true[1] - 2
    
    idxs.outside <- which(!inRegion3(X, type.cov) & !inRegion4(X, type.cov))
    C <- rep(0, n)
    C[idxs.outside] <- Lambda_inverse_AFT_ll(runif(length(idxs.outside))) - X[idxs.outside, ] %*% gamma
    C[idxs.region3] <- Lambda_inverse_AFT_ll(runif(length(idxs.region3))) - X[idxs.region3, ] %*% gamma.region3
    C[idxs.region4] <- Lambda_inverse_AFT_ll(runif(length(idxs.region4))) - X[idxs.region4, ] %*% gamma.region4
    
  } else if (DGP %% 20 == 19) { # Pos. dep., C ~ AFT_ll, high cens, regions with less cens 3
    
    # Specify copula
    cop_to_use <- frankCopula(param = 6, dim = 2)
    
    # Generate data from copula
    U <- rCopula(n, cop_to_use)
    u1 <- U[,1]
    u2 <- U[,2]
    
    # Latent event time T
    T <- inv.Lambda(u1) - X %*% beta.true
    
    # Latent censoring time C
    gamma <- beta.true
    gamma.region3 <- beta.true
    gamma.region4 <- beta.true
    gamma[1] <- beta.true[1] + 0.25
    gamma.region3[1] <- beta.true[1] - 3
    gamma.region4[1] <- beta.true[1] - 2
    
    idxs.outside <- which(!inRegion3(X, type.cov) & !inRegion4(X, type.cov))
    C <- rep(0, n)
    C[idxs.outside] <- Lambda_inverse_AFT_ll(u2[idxs.outside]) - X[idxs.outside, ] %*% gamma
    C[idxs.region3] <- Lambda_inverse_AFT_ll(u2[idxs.region3]) - X[idxs.region3, ] %*% gamma.region3
    C[idxs.region4] <- Lambda_inverse_AFT_ll(u2[idxs.region4]) - X[idxs.region4, ] %*% gamma.region4
    
  } 
  
  # Observed identified minimum
  Y <- pmin(T, C)
  Delta <- as.numeric(Y == T)
  
  # Collect all variables
  data <- data.frame(Y, Delta, X)
  colnames(data) <- c("Y", "Delta", colnames(X))
  
  # Histogram of the observed times
  if (plot.data) {
    print(sprintf("Percentage of censored observations: %.2f%%",
                  100*(1 - sum(data$Delta)/n)))
    hist(Y)
  }
  
  # Return the results
  data
}

#' @title Additional data generating function.
#' 
#' @description This function generates a data set according to the specified
#' arguments, like 'generateData.R' above. It differs from the aforementioned
#' function in that some DGP's are slightly different
#' 
#' @param beta.true True covariate vector, as a function of time.
#' @param n Sample size.
#' @param n.cov Number of covariates.
#' @param options List of additional arguments.
#' @param plot.data Boolean value indicating whether or not to plot the
#' generated data set. Default value is \code{plot.data = FALSE}.
#'
#' @import stats copula
#' 
#' @noRd
#' 
generateData_add <- function(beta.true, n, n.cov, options, plot.data = FALSE) {
  
  # Extract the necessary hyperparameters
  if (options[["link.function"]] == "AFT_ll") {
    inv.Lambda <- Lambda_inverse_AFT_ll
  } else if (options[["link.function"]] == "Cox_wb") {
    inv.Lambda <- Lambda_inverse_Cox_wb
  }
  DGP <- options[["DGP"]]
  
  # Set the types of the covariates to be used based on the selected DGP.
  #
  # NOTE TO SELF: If you change these thresholds (20, 40, 60), also change them
  #               accordingly in 'simulate1D.R', 'set.hyperparameters.R' and
  #               'simFuncWrapper.R'
  if (DGP <= 20) {
    type.cov <- rep("c", n.cov)
  } else if (20 < DGP & DGP <= 40) {
    type.cov <- rep("c", n.cov)
    if ((ceiling(n.cov/2) + 1) <= length(type.cov)) {
      type.cov[(ceiling(n.cov/2) + 1):length(type.cov)] <- "b"
    }
  } else {
    type.cov <- rep("b", n.cov)
  }
  
  # Because, at least for now, beta.true only depends on t in its first element
  # (corresponding to the intercept) and moreover, this dependence is of the
  # form beta_0(t) = t + a, we can make the following short-cut.
  beta.true <- beta.true(0)
  
  # Subset beta.true to the correct number of parameters
  beta.true <- beta.true[1:(n.cov + 1)]
  
  # Generate the intercept and the covariates
  X <- rep(1, n)
  for (i in 1:n.cov) {
    if (type.cov[i] == "c") {
      X <- cbind(X, runif(n, -3, 3))
    } else if (type.cov[i] == "b") {
      X <- cbind(X, rbinom(n, 1, 0.5))
    } else {
      stop("Invalid value for type of covariate")
    }
  }
  colnames(X) <- paste0("X", 0:n.cov)
  
  # Determine which covariate vectors lie in the predefined regions.
  idxs.region1 <- which(inRegion1(X, type.cov))
  idxs.region2 <- which(inRegion2(X, type.cov))
  idxs.region3 <- which(inRegion3(X, type.cov))
  idxs.region4 <- which(inRegion4(X, type.cov))
  idxs.region5 <- which(inRegion5(X, type.cov))
  idxs.region6 <- which(inRegion6(X, type.cov))
  
  # Set dependence structure based on selected DGP
  if (DGP %% 20 == 1) { # Independence, region 5
    
    # Latent event time T
    T <- inv.Lambda(runif(n)) - X %*% beta.true
    
    # Latent censoring time C
    gamma <- beta.true
    gamma.region5 <- beta.true
    gamma[1] <- beta.true[1] + 0.4
    gamma.region5[1] <- beta.true[1] - 5
    
    idxs.outside <- which(!inRegion5(X, type.cov))
    C <- rep(0, n)
    C[idxs.outside] <- Lambda_inverse_AFT_ll(runif(length(idxs.outside))) - X[idxs.outside, ] %*% gamma
    C[idxs.region5] <- Lambda_inverse_AFT_ll(runif(length(idxs.region5))) - X[idxs.region5, ] %*% gamma.region5

  } else if (DGP %% 20 == 2) { # Positive dependence, region 5
    
    # Specify copula
    cop_to_use <- frankCopula(param = 6, dim = 2)
    
    # Generate data from copula
    U <- rCopula(n, cop_to_use)
    u1 <- U[,1]
    u2 <- U[,2]
    
    # Latent event time T
    T <- inv.Lambda(u1) - X %*% beta.true
    
    # Latent censoring time C
    gamma <- beta.true
    gamma.region5 <- beta.true
    gamma[1] <- beta.true[1] + 0
    gamma.region5[1] <- beta.true[1] - 5
    
    idxs.outside <- which(!inRegion5(X, type.cov))
    C <- rep(0, n)
    C[idxs.outside] <- Lambda_inverse_AFT_ll(u2[idxs.outside]) - X[idxs.outside, ] %*% gamma
    C[idxs.region5] <- Lambda_inverse_AFT_ll(u2[idxs.region5]) - X[idxs.region5, ] %*% gamma.region5
    
  } else if (DGP %% 20 == 3) { # Negative dependence, region 5
    
    # Specify copula
    cop_to_use <- frankCopula(param = -6, dim = 2)
    
    # Generate data from copula
    U <- rCopula(n, cop_to_use)
    u1 <- U[,1]
    u2 <- U[,2]
    
    # Latent event time T
    T <- inv.Lambda(u1) - X %*% beta.true
    
    # Latent censoring time C
    gamma <- beta.true
    gamma.region5 <- beta.true
    gamma[1] <- beta.true[1] + 0.2
    gamma.region5[1] <- beta.true[1] - 5
    
    idxs.outside <- which(!inRegion5(X, type.cov))
    C <- rep(0, n)
    C[idxs.outside] <- Lambda_inverse_AFT_ll(u2[idxs.outside]) - X[idxs.outside, ] %*% gamma
    C[idxs.region5] <- Lambda_inverse_AFT_ll(u2[idxs.region5]) - X[idxs.region5, ] %*% gamma.region5
    
  } else if (DGP %% 20 == 4) { # Independence, region 6
    
    # Latent event time T
    T <- inv.Lambda(runif(n)) - X %*% beta.true
    
    # Latent censoring time C
    gamma <- beta.true
    gamma.region6 <- beta.true
    gamma[1] <- beta.true[1] + 0.4
    gamma.region6[1] <- beta.true[1] - 5
    
    idxs.outside <- which(!inRegion6(X, type.cov))
    C <- rep(0, n)
    C[idxs.outside] <- Lambda_inverse_AFT_ll(runif(length(idxs.outside))) - X[idxs.outside, ] %*% gamma
    C[idxs.region6] <- Lambda_inverse_AFT_ll(runif(length(idxs.region6))) - X[idxs.region6, ] %*% gamma.region6
    
  } else if (DGP %% 20 == 5) { # Positive dependence, region 6
    
    # Specify copula
    cop_to_use <- frankCopula(param = 6, dim = 2)
    
    # Generate data from copula
    U <- rCopula(n, cop_to_use)
    u1 <- U[,1]
    u2 <- U[,2]
    
    # Latent event time T
    T <- inv.Lambda(u1) - X %*% beta.true
    
    # Latent censoring time C
    gamma <- beta.true
    gamma.region6 <- beta.true
    gamma[1] <- beta.true[1] + 0.05
    gamma.region6[1] <- beta.true[1] - 5
    
    idxs.outside <- which(!inRegion6(X, type.cov))
    C <- rep(0, n)
    C[idxs.outside] <- Lambda_inverse_AFT_ll(u2[idxs.outside]) - X[idxs.outside, ] %*% gamma
    C[idxs.region6] <- Lambda_inverse_AFT_ll(u2[idxs.region6]) - X[idxs.region6, ] %*% gamma.region6
    
  } else if (DGP %% 20 == 6) { # Negative dependence, region 6
    
    # Specify copula
    cop_to_use <- frankCopula(param = -6, dim = 2)
    
    # Generate data from copula
    U <- rCopula(n, cop_to_use)
    u1 <- U[,1]
    u2 <- U[,2]
    
    # Latent event time T
    T <- inv.Lambda(u1) - X %*% beta.true
    
    # Latent censoring time C
    gamma <- beta.true
    gamma.region6 <- beta.true
    gamma[1] <- beta.true[1] + 0.2
    gamma.region6[1] <- beta.true[1] - 5
    
    idxs.outside <- which(!inRegion6(X, type.cov))
    C <- rep(0, n)
    C[idxs.outside] <- Lambda_inverse_AFT_ll(u2[idxs.outside]) - X[idxs.outside, ] %*% gamma
    C[idxs.region6] <- Lambda_inverse_AFT_ll(u2[idxs.region6]) - X[idxs.region6, ] %*% gamma.region6
    
  } else if (DGP %% 20 == 7) { # Independence, regions 5 and 6
    
    # Latent event time T
    T <- inv.Lambda(runif(n)) - X %*% beta.true
    
    # Latent censoring time C
    gamma <- beta.true
    gamma.region5 <- gamma.region6 <- beta.true
    gamma[1] <- beta.true[1] + 0.8
    gamma.region5[1] <- gamma.region6[1] <- beta.true[1] - 5
    
    idxs.outside <- which(!inRegion5(X, type.cov) & !inRegion6(X, type.cov))
    C <- rep(0, n)
    C[idxs.outside] <- Lambda_inverse_AFT_ll(runif(length(idxs.outside))) - X[idxs.outside, ] %*% gamma
    C[idxs.region5] <- Lambda_inverse_AFT_ll(runif(length(idxs.region5))) - X[idxs.region5, ] %*% gamma.region5
    C[idxs.region6] <- Lambda_inverse_AFT_ll(runif(length(idxs.region6))) - X[idxs.region6, ] %*% gamma.region6
    
  } else if (DGP %% 20 == 8) { # Positive dependence, regions 5 and 6
    
    # Specify copula
    cop_to_use <- frankCopula(param = 6, dim = 2)
    
    # Generate data from copula
    U <- rCopula(n, cop_to_use)
    u1 <- U[,1]
    u2 <- U[,2]
    
    # Latent event time T
    T <- inv.Lambda(u1) - X %*% beta.true
    
    # Latent censoring time C
    gamma <- beta.true
    gamma.region5 <- gamma.region6 <- beta.true
    gamma[1] <- beta.true[1] + 0.2
    gamma.region5[1] <- gamma.region6[1] <- beta.true[1] - 5
    
    idxs.outside <- which(!inRegion5(X, type.cov) & !inRegion6(X, type.cov))
    C <- rep(0, n)
    C[idxs.outside] <- Lambda_inverse_AFT_ll(u2[idxs.outside]) - X[idxs.outside, ] %*% gamma
    C[idxs.region5] <- Lambda_inverse_AFT_ll(u2[idxs.region5]) - X[idxs.region5, ] %*% gamma.region5
    C[idxs.region6] <- Lambda_inverse_AFT_ll(u2[idxs.region6]) - X[idxs.region6, ] %*% gamma.region6
    
  } else if (DGP %% 20 == 9) { # Negative dependence, regions 5 and 6
    
    # Specify copula
    cop_to_use <- frankCopula(param = -6, dim = 2)
    
    # Generate data from copula
    U <- rCopula(n, cop_to_use)
    u1 <- U[,1]
    u2 <- U[,2]
    
    # Latent event time T
    T <- inv.Lambda(u1) - X %*% beta.true
    
    # Latent censoring time C
    gamma <- beta.true
    gamma.region5 <- gamma.region6 <- beta.true
    gamma[1] <- beta.true[1] + 0.7
    gamma.region5[1] <- gamma.region6[1] <- beta.true[1] - 5
    
    idxs.outside <- which(!inRegion5(X, type.cov) & !inRegion6(X, type.cov))
    C <- rep(0, n)
    C[idxs.outside] <- Lambda_inverse_AFT_ll(u2[idxs.outside]) - X[idxs.outside, ] %*% gamma
    C[idxs.region5] <- Lambda_inverse_AFT_ll(u2[idxs.region5]) - X[idxs.region5, ] %*% gamma.region5
    C[idxs.region6] <- Lambda_inverse_AFT_ll(u2[idxs.region6]) - X[idxs.region6, ] %*% gamma.region6
    
  } else if (DGP %% 20 == 10) { # Independence, regions 1 + 4
    
    # Latent event time T
    T <- inv.Lambda(runif(n)) - X %*% beta.true
    
    # Latent censoring time C
    gamma <- beta.true
    gamma.region1 <- beta.true
    gamma.region4 <- beta.true
    gamma[1] <- beta.true[1] + 0.2
    gamma.region1[1] <- beta.true[1] - 5
    gamma.region4[1] <- beta.true[1] - 6
    
    idxs.outside <- which(!inRegion1(X, type.cov) & !inRegion4(X, type.cov))
    C <- rep(0, n)
    C[idxs.outside] <- Lambda_inverse_AFT_ll(runif(length(idxs.outside))) - X[idxs.outside, ] %*% gamma
    C[idxs.region1] <- Lambda_inverse_AFT_ll(runif(length(idxs.region1))) - X[idxs.region1, ] %*% gamma.region1
    C[idxs.region4] <- Lambda_inverse_AFT_ll(runif(length(idxs.region4))) - X[idxs.region4, ] %*% gamma.region4
    
  } else if (DGP %% 20 == 11) { # Positive dependence, regions 1 + 4
    
    # Specify copula
    cop_to_use <- frankCopula(param = 6, dim = 2)
    
    # Generate data from copula
    U <- rCopula(n, cop_to_use)
    u1 <- U[,1]
    u2 <- U[,2]
    
    # Latent event time T
    T <- inv.Lambda(u1) - X %*% beta.true
    
    # Latent censoring time C
    gamma <- beta.true
    gamma.region1 <- beta.true
    gamma.region4 <- beta.true
    gamma[1] <- beta.true[1] 
    gamma.region1[1] <- beta.true[1] - 4
    gamma.region4[1] <- beta.true[1] - 5
    
    idxs.outside <- which(!inRegion1(X, type.cov) & !inRegion4(X, type.cov))
    C <- rep(0, n)
    C[idxs.outside] <- Lambda_inverse_AFT_ll(u2[idxs.outside]) - X[idxs.outside, ] %*% gamma
    C[idxs.region1] <- Lambda_inverse_AFT_ll(u2[idxs.region1]) - X[idxs.region1, ] %*% gamma.region1
    C[idxs.region4] <- Lambda_inverse_AFT_ll(u2[idxs.region4]) - X[idxs.region4, ] %*% gamma.region4
    
  } else if (DGP %% 20 == 12) { # Negative dependence, regions 1 + 4
    
    # Specify copula
    cop_to_use <- frankCopula(param = -6, dim = 2)
    
    # Generate data from copula
    U <- rCopula(n, cop_to_use)
    u1 <- U[,1]
    u2 <- U[,2]
    
    # Latent event time T
    T <- inv.Lambda(u1) - X %*% beta.true
    
    # Latent censoring time C
    gamma <- beta.true
    gamma.region1 <- beta.true
    gamma.region4 <- beta.true
    gamma[1] <- beta.true[1] + 0.1
    gamma.region1[1] <- beta.true[1] - 6
    gamma.region4[1] <- beta.true[1] - 7
    
    idxs.outside <- which(!inRegion1(X, type.cov) & !inRegion4(X, type.cov))
    C <- rep(0, n)
    C[idxs.outside] <- Lambda_inverse_AFT_ll(u2[idxs.outside]) - X[idxs.outside, ] %*% gamma
    C[idxs.region1] <- Lambda_inverse_AFT_ll(u2[idxs.region1]) - X[idxs.region1, ] %*% gamma.region1
    C[idxs.region4] <- Lambda_inverse_AFT_ll(u2[idxs.region4]) - X[idxs.region4, ] %*% gamma.region4
    
  } else if (DGP %% 20 == 13) { # Independence, regions 3 + 4
    
    # Latent event time T
    T <- inv.Lambda(runif(n)) - X %*% beta.true
    
    # Latent censoring time C
    gamma <- beta.true
    gamma.region3 <- beta.true
    gamma.region4 <- beta.true
    gamma[1] <- beta.true[1] + 0.2
    gamma.region3[1] <- beta.true[1] - 5
    gamma.region4[1] <- beta.true[1] - 6
    
    idxs.outside <- which(!inRegion3(X, type.cov) & !inRegion4(X, type.cov))
    C <- rep(0, n)
    C[idxs.outside] <- Lambda_inverse_AFT_ll(runif(length(idxs.outside))) - X[idxs.outside, ] %*% gamma
    C[idxs.region3] <- Lambda_inverse_AFT_ll(runif(length(idxs.region3))) - X[idxs.region3, ] %*% gamma.region3
    C[idxs.region4] <- Lambda_inverse_AFT_ll(runif(length(idxs.region4))) - X[idxs.region4, ] %*% gamma.region4
    
  } else if (DGP %% 20 == 14) { # Positive dependence, regions 3 + 4
    
    # Specify copula
    cop_to_use <- frankCopula(param = 6, dim = 2)
    
    # Generate data from copula
    U <- rCopula(n, cop_to_use)
    u1 <- U[,1]
    u2 <- U[,2]
    
    # Latent event time T
    T <- inv.Lambda(u1) - X %*% beta.true
    
    # Latent censoring time C
    gamma <- beta.true
    gamma.region3 <- beta.true
    gamma.region4 <- beta.true
    gamma[1] <- beta.true[1]
    gamma.region3[1] <- beta.true[1] - 5
    gamma.region4[1] <- beta.true[1] - 5
    
    idxs.outside <- which(!inRegion3(X, type.cov) & !inRegion4(X, type.cov))
    C <- rep(0, n)
    C[idxs.outside] <- Lambda_inverse_AFT_ll(u2[idxs.outside]) - X[idxs.outside, ] %*% gamma
    C[idxs.region3] <- Lambda_inverse_AFT_ll(u2[idxs.region3]) - X[idxs.region3, ] %*% gamma.region3
    C[idxs.region4] <- Lambda_inverse_AFT_ll(u2[idxs.region4]) - X[idxs.region4, ] %*% gamma.region4
    
  } else if (DGP %% 20 == 15) { # Negative dependence, regions 3 + 4
    
    # Specify copula
    cop_to_use <- frankCopula(param = -6, dim = 2)
    
    # Generate data from copula
    U <- rCopula(n, cop_to_use)
    u1 <- U[,1]
    u2 <- U[,2]
    
    # Latent event time T
    T <- inv.Lambda(u1) - X %*% beta.true
    
    # Latent censoring time C
    gamma <- beta.true
    gamma.region3 <- beta.true
    gamma.region4 <- beta.true
    gamma[1] <- beta.true[1]
    gamma.region3[1] <- beta.true[1] - 5
    gamma.region4[1] <- beta.true[1] - 5
    
    idxs.outside <- which(!inRegion3(X, type.cov) & !inRegion4(X, type.cov))
    C <- rep(0, n)
    C[idxs.outside] <- Lambda_inverse_AFT_ll(u2[idxs.outside]) - X[idxs.outside, ] %*% gamma
    C[idxs.region3] <- Lambda_inverse_AFT_ll(u2[idxs.region3]) - X[idxs.region3, ] %*% gamma.region3
    C[idxs.region4] <- Lambda_inverse_AFT_ll(u2[idxs.region4]) - X[idxs.region4, ] %*% gamma.region4
    
  } else if (DGP %% 20 == 16) { # Independence, C ~ AFT_ll, high cens. ~ DGP 7 in generateData
    
    # Latent event time T
    T <- inv.Lambda(runif(n)) - X %*% beta.true
    
    # Latent censoring time C
    gamma <- beta.true
    C <- Lambda_inverse_AFT_ll(runif(n)) - X %*% gamma
    
  } else if (DGP %% 20 == 17) { # Pos. dep., C ~ AFT_ll, high cens. ~ DGP 13 in generateData
    
    # Specify copula
    cop_to_use <- frankCopula(param = 6, dim = 2)
    
    # Generate data from copula
    U <- rCopula(n, cop_to_use)
    u1 <- U[,1]
    u2 <- U[,2]
    
    # Latent event time T
    T <- inv.Lambda(u1) - X %*% beta.true
    
    # Latent censoring time C
    gamma <- beta.true
    gamma[1] <- beta.true[1] - 0.115
    C <- Lambda_inverse_AFT_ll(u2) - X %*% gamma
    
  } else if (DGP %% 20 == 18) { # Neg. dep., C ~ AFT_ll, high cens. ~ DGP 14 in generateData
    
    # Specify copula
    cop_to_use <- frankCopula(param = -6, dim = 2)
    
    # Generate data from copula
    U <- rCopula(n, cop_to_use)
    u1 <- U[,1]
    u2 <- U[,2]
    
    # Latent event time T
    T <- inv.Lambda(u1) - X %*% beta.true
    
    # Latent censoring time C
    gamma <- beta.true
    gamma[1] <- beta.true[1] - 0.2
    C <- Lambda_inverse_AFT_ll(u2) - X %*% gamma
    
  }
  
  # Observed identified minimum
  Y <- pmin(T, C)
  Delta <- as.numeric(Y == T)
  
  # Collect all variables
  data <- data.frame(Y, Delta, X)
  colnames(data) <- c("Y", "Delta", colnames(X))
  
  # Histogram of the observed times
  if (plot.data) {
    print(sprintf("Percentage of censored observations: %.2f%%",
                  100*(1 - sum(data$Delta)/n)))
    hist(Y)
  }
  
  # Return the results
  data
}

#' @title Data generation function for the main simulation.
#' 
#' @description This function generates a data set according to the specified
#' arguments.
#' 
#' @param beta.true True covariate vector, as a function of time.
#' @param n Sample size.
#' @param n.cov Number of covariates.
#' @param options List of additional arguments.
#' @param H0.inv Inverse of the intercept (function of time).
#' @param plot.data Boolean value indicating whether or not to plot the
#' generated data set. Default value is \code{plot.data = FALSE}.
#' 
#' @import stats copula
#' 
#' @noRd
#' 
generateData_simMain <- function(beta.true, n, n.cov, options, H0.inv,
                                 plot.data = FALSE) {
  
  # Extract the necessary hyperparameters
  if (options[["link.function"]] == "AFT_ll") {
    inv.Lambda <- Lambda_inverse_AFT_ll
  } else if (options[["link.function"]] == "Cox_wb") {
    inv.Lambda <- Lambda_inverse_Cox_wb
  }
  DGP <- options[["DGP"]]
  
  # Subset the parameter vector to the first n.cov covariate effects.
  if (class(beta.true) == "function") {
    beta.true <- beta.true(0)[2:(n.cov + 1)]
  }
  
  # Set the types of the covariates to be used based on the selected DGP.
  #
  # NOTE TO SELF: If you change these thresholds (20, 40, 60), also change them
  #               accordingly in 'simulate1D.R', 'set.hyperparameters.R', 
  #               'simFuncWrapper.R' and 'simulate1D.CCDF.R'
  if (DGP <= 20) {
    type.cov <- rep("c", n.cov)
  } else if (20 < DGP & DGP <= 40) {
    type.cov <- rep("c", n.cov)
    if ((ceiling(n.cov/2) + 1) <= length(type.cov)) {
      type.cov[(ceiling(n.cov/2) + 1):length(type.cov)] <- "b"
    }
  } else {
    type.cov <- rep("b", n.cov)
  }
  
  # Generate the intercept and the covariates
  X <- rep(1, n)
  for (i in 1:n.cov) {
    if (type.cov[i] == "c") {
      X <- cbind(X, rnorm(n))
    } else if (type.cov[i] == "b") {
      X <- cbind(X, rbinom(n, 1, 0.5))
    } else {
      stop("Invalid value for type of covariate")
    }
  }
  colnames(X) <- paste0("X", 0:n.cov)
  
  # Get matrix of just the covariates (no intercept)
  X.noint <- X[, -1]
  
  # Set dependence structure based on selected DGP
  if (DGP %% 20 == 1) { # For AFT_ll: Independence, C ~ Exp, ~25% cens.
    
    # Latent event time T
    T <- H0.inv(exp(inv.Lambda(runif(n)) - X.noint %*% beta.true))
    
    # Latent censoring time C
    lambda.c <- 0.1
    C <- (-1/lambda.c) * log(1 - runif(n))
    
  } else if (DGP %% 20 == 2) { # For AFT_ll: Pos. dep., C ~ Exp, ~25% cens.
    
    # Specify copula
    cop_to_use <- frankCopula(param = 6, dim = 2)
    
    # Generate data from copula
    U <- rCopula(n, cop_to_use)
    u1 <- U[,1]
    u2 <- U[,2]
    
    # Latent event time T
    T <- H0.inv(exp(inv.Lambda(u1) - X.noint %*% beta.true))
    
    # Latent censoring time C
    lambda.c <- 0.15
    C <- (-1/lambda.c) * log(1 - u2)
    
  } else if (DGP %% 20 == 3) { # For AFT_ll: Neg. dep., C ~ Exp, ~25% cens.
    
    # Specify copula
    cop_to_use <- frankCopula(param = -6, dim = 2)
    
    # Generate data from copula
    U <- rCopula(n, cop_to_use)
    u1 <- U[,1]
    u2 <- U[,2]
    
    # Latent event time T
    T <- H0.inv(exp(inv.Lambda(u1) - X.noint %*% beta.true))
    
    # Latent censoring time C
    lambda.c <- 0.07
    C <- (-1/lambda.c) * log(1 - u2)
    
  } else if (DGP %% 20 == 4) { # For AFT_ll: Independence, C ~ Exp, ~65% cens.
    
    # Latent event time T
    T <- H0.inv(exp(inv.Lambda(runif(n)) - X.noint %*% beta.true))
    
    # Latent censoring time C
    lambda.c <- 0.8
    C <- (-1/lambda.c) * log(1 - runif(n))
    
  } else if (DGP %% 20 == 5) { # For AFT_ll: Pos. dep., C ~ Exp, ~65% cens.
    
    # Specify copula
    cop_to_use <- frankCopula(param = 6, dim = 2)
    
    # Generate data from copula
    U <- rCopula(n, cop_to_use)
    u1 <- U[,1]
    u2 <- U[,2]
    
    # Latent event time T
    T <- H0.inv(exp(inv.Lambda(u1) - X.noint %*% beta.true))
    
    # Latent censoring time C
    lambda.c <- 0.65
    C <- (-1/lambda.c) * log(1 - u2)
    
  } else if (DGP %% 20 == 6) { # For AFT_ll: Neg. dep., C ~ Exp, ~65% cens.
    
    # Specify copula
    cop_to_use <- frankCopula(param = -6, dim = 2)
    
    # Generate data from copula
    U <- rCopula(n, cop_to_use)
    u1 <- U[,1]
    u2 <- U[,2]
    
    # Latent event time T
    T <- H0.inv(exp(inv.Lambda(u1) - X.noint %*% beta.true))
    
    # Latent censoring time C
    lambda.c <- 1.2
    C <- (-1/lambda.c) * log(1 - u2)
    
  } else if (DGP %% 20 == 7) { # For Cox_wb: Independence, C ~ Exp, ~25% cens.
    
    # Latent event time T
    T <- H0.inv(exp(inv.Lambda(runif(n)) - X.noint %*% beta.true))
    
    # Latent censoring time C
    lambda.c <- 0.22
    C <- (-1/lambda.c) * log(1 - runif(n))
    
  } else if (DGP %% 20 == 8) { # For Cox_wb: Pos. dep., C ~ Exp, ~25% cens.
    
    # Specify copula
    cop_to_use <- frankCopula(param = 6, dim = 2)
    
    # Generate data from copula
    U <- rCopula(n, cop_to_use)
    u1 <- U[,1]
    u2 <- U[,2]
    
    # Latent event time T
    T <- H0.inv(exp(inv.Lambda(u1) - X.noint %*% beta.true))
    
    # Latent censoring time C
    lambda.c <- 0.3
    C <- (-1/lambda.c) * log(1 - u2)
    
  } else if (DGP %% 20 == 9) { # For Cox_wb: Neg. dep., C ~ Exp, ~25% cens.
    
    # Specify copula
    cop_to_use <- frankCopula(param = -6, dim = 2)
    
    # Generate data from copula
    U <- rCopula(n, cop_to_use)
    u1 <- U[,1]
    u2 <- U[,2]
    
    # Latent event time T
    T <- H0.inv(exp(inv.Lambda(u1) - X.noint %*% beta.true))
    
    # Latent censoring time C
    lambda.c <- 0.17
    C <- (-1/lambda.c) * log(1 - u2)
    
  } else if (DGP %% 20 == 10) { # For Cox_wb: Independence, C ~ Exp, ~65% cens.
    
    # Latent event time T
    T <- H0.inv(exp(inv.Lambda(runif(n)) - X.noint %*% beta.true))
    
    # Latent censoring time C
    lambda.c <- 1.3
    C <- (-1/lambda.c) * log(1 - runif(n))
    
  } else if (DGP %% 20 == 11) { # For Cox_wb: Pos. dep., C ~ Exp, ~65% cens.
    
    # Specify copula
    cop_to_use <- frankCopula(param = 6, dim = 2)
    
    # Generate data from copula
    U <- rCopula(n, cop_to_use)
    u1 <- U[,1]
    u2 <- U[,2]
    
    # Latent event time T
    T <- H0.inv(exp(inv.Lambda(u1) - X.noint %*% beta.true))
    
    # Latent censoring time C
    lambda.c <- 1
    C <- (-1/lambda.c) * log(1 - u2)
    
  } else if (DGP %% 20 == 12) { # For Cox_wb: Neg. dep., C ~ Exp, ~65% cens.
    
    # Specify copula
    cop_to_use <- frankCopula(param = -6, dim = 2)
    
    # Generate data from copula
    U <- rCopula(n, cop_to_use)
    u1 <- U[,1]
    u2 <- U[,2]
    
    # Latent event time T
    T <- H0.inv(exp(inv.Lambda(u1) - X.noint %*% beta.true))
    
    # Latent censoring time C
    lambda.c <- 1.6
    C <- (-1/lambda.c) * log(1 - u2)
    
  }  else if (DGP %% 20 == 13) { # For AFT_ll: Independence, C ~ Exp, ~2% cens.
    
    # Latent event time T
    T <- H0.inv(exp(inv.Lambda(runif(n)) - X.noint %*% beta.true))
    
    # Latent censoring time C
    lambda.c <- 0.002
    C <- (-1/lambda.c) * log(1 - runif(n))
    
  } else if (DGP %% 20 == 14) { # For AFT_ll: Pos. dep., C ~ Exp, ~2% cens.
    
    # Specify copula
    cop_to_use <- frankCopula(param = 6, dim = 2)
    
    # Generate data from copula
    U <- rCopula(n, cop_to_use)
    u1 <- U[,1]
    u2 <- U[,2]
    
    # Latent event time T
    T <- H0.inv(exp(inv.Lambda(u1) - X.noint %*% beta.true))
    
    # Latent censoring time C
    lambda.c <- 0.006
    C <- (-1/lambda.c) * log(1 - u2)
    
  } else if (DGP %% 20 == 15) { # For AFT_ll: Neg. dep., C ~ Exp, ~2% cens.
    
    # Specify copula
    cop_to_use <- frankCopula(param = -6, dim = 2)
    
    # Generate data from copula
    U <- rCopula(n, cop_to_use)
    u1 <- U[,1]
    u2 <- U[,2]
    
    # Latent event time T
    T <- H0.inv(exp(inv.Lambda(u1) - X.noint %*% beta.true))
    
    # Latent censoring time C
    lambda.c <- 0.0005
    C <- (-1/lambda.c) * log(1 - u2)
    
  } else if (DGP %% 20 == 16) { # For Cox_wb: Independence, C ~ Exp, ~2% cens.
    
    # Latent event time T
    T <- H0.inv(exp(inv.Lambda(runif(n)) - X.noint %*% beta.true))
    
    # Latent censoring time C
    lambda.c <- 0.005
    C <- (-1/lambda.c) * log(1 - runif(n))
    
  } else if (DGP %% 20 == 17) { # For Cox_wb: Pos. dep., C ~ Exp, ~2% cens.
    
    # Specify copula
    cop_to_use <- frankCopula(param = 6, dim = 2)
    
    # Generate data from copula
    U <- rCopula(n, cop_to_use)
    u1 <- U[,1]
    u2 <- U[,2]
    
    # Latent event time T
    T <- H0.inv(exp(inv.Lambda(u1) - X.noint %*% beta.true))
    
    # Latent censoring time C
    lambda.c <- 0.03
    C <- (-1/lambda.c) * log(1 - u2)
    
  } else if (DGP %% 20 == 18) { # For Cox_wb: Neg. dep., C ~ Exp, ~2% cens.
    
    # Specify copula
    cop_to_use <- frankCopula(param = -6, dim = 2)
    
    # Generate data from copula
    U <- rCopula(n, cop_to_use)
    u1 <- U[,1]
    u2 <- U[,2]
    
    # Latent event time T
    T <- H0.inv(exp(inv.Lambda(u1) - X.noint %*% beta.true))
    
    # Latent censoring time C
    lambda.c <- 0.003
    C <- (-1/lambda.c) * log(1 - u2)
    
  }
  
  # Observed identified minimum
  Y <- pmin(T, C)
  Delta <- as.numeric(Y == T)
  
  # Collect all variables
  data <- data.frame(Y, Delta, X)
  colnames(data) <- c("Y", "Delta", colnames(X))
  
  # Histogram of the observed times
  if (plot.data) {
    print(sprintf("Percentage of censored observations: %.2f%%",
                  100*(1 - sum(data$Delta)/n)))
    hist(Y)
  }
  
  # Return the results
  data
}

#' @title Data generation function for the additional simulation.
#' 
#' @description This function generates a data set according to the specified
#' arguments. As opposed to \code{generateData_simMain.R}, it generates the
#' covariates in a way that they are depedendent. (Achieved through the use of
#' a copula).
#' 
#' @param beta.true True covariate vector, as a function of time.
#' @param n Sample size.
#' @param n.cov Number of covariates.
#' @param options List of additional arguments.
#' @param H0.inv Inverse of the intercept (function of time).
#' @param plot.data Boolean value indicating whether or not to plot the
#' generated data set. Default value is \code{plot.data = FALSE}.
#' 
#' @import stats copula
#' 
#' @noRd
#' 
generateData_simAdd <- function(beta.true, n, n.cov, options, H0.inv,
                                plot.data = FALSE) {
  
  # Extract the necessary hyperparameters
  if (options[["link.function"]] == "AFT_ll") {
    inv.Lambda <- Lambda_inverse_AFT_ll
  } else if (options[["link.function"]] == "Cox_wb") {
    inv.Lambda <- Lambda_inverse_Cox_wb
  }
  DGP <- options[["DGP"]]
  
  # Subset the parameter vector to the first n.cov covariate effects.
  if (class(beta.true) == "function") {
    beta.true <- beta.true(0)[2:(n.cov + 1)]
  }
  
  # Set the types of the covariates to be used based on the selected DGP.
  #
  # NOTE TO SELF: If you change these thresholds (20, 40, 60), also change them
  #               accordingly in 'simulate1D.R', 'set.hyperparameters.R', 
  #               'simFuncWrapper.R' and 'simulate1D.CCDF.R'
  if (DGP <= 20) {
    type.cov <- rep("c", n.cov)
  } else if (20 < DGP & DGP <= 40) {
    type.cov <- rep("c", n.cov)
    if ((ceiling(n.cov/2) + 1) <= length(type.cov)) {
      type.cov[(ceiling(n.cov/2) + 1):length(type.cov)] <- "b"
    }
  } else {
    type.cov <- rep("b", n.cov)
  }
  
  # Throw an error if specified number of covariates is not equal to two. These
  # cases are not implemented in this function.
  if (n.cov != 2) {
    stop("The specified number of covariates should equal 2.")
  }
  
  # Generate the intercept and the covariates
  X <- rep(1, n)
  X.U <- rCopula(n, normalCopula(0.8))
  for (i in 1:n.cov) {
    if (type.cov[i] == "c") {
      X <- cbind(X, qnorm(X.U[,i]))
    } else if (type.cov[i] == "b") {
      X <- cbind(X, as.numeric(X.U[,i] >= 0.5))
    } else {
      stop("Invalid value for type of covariate")
    }
  }
  colnames(X) <- paste0("X", 0:n.cov)
  
  # Get matrix of just the covariates (no intercept)
  X.noint <- X[, -1]
  
  # Set dependence structure based on selected DGP
  if (DGP %% 20 == 1) { # For AFT_ll: Independence, C ~ Exp, ~30% cens.
    
    # Latent event time T
    T <- H0.inv(exp(inv.Lambda(runif(n)) - X.noint %*% beta.true))
    
    # Latent censoring time C
    lambda.c <- 0.13
    C <- (-1/lambda.c) * log(1 - runif(n))
    
  } else if (DGP %% 20 == 2) { # For AFT_ll: Pos. dep., C ~ Exp, ~25% cens.
    
    # Specify copula
    cop_to_use <- frankCopula(param = 6, dim = 2)
    
    # Generate data from copula
    U <- rCopula(n, cop_to_use)
    u1 <- U[,1]
    u2 <- U[,2]
    
    # Latent event time T
    T <- H0.inv(exp(inv.Lambda(u1) - X.noint %*% beta.true))
    
    # Latent censoring time C
    lambda.c <- 0.17
    C <- (-1/lambda.c) * log(1 - u2)
    
  } else if (DGP %% 20 == 3) { # For AFT_ll: Neg. dep., C ~ Exp, ~25% cens.
    
    # Specify copula
    cop_to_use <- frankCopula(param = -6, dim = 2)
    
    # Generate data from copula
    U <- rCopula(n, cop_to_use)
    u1 <- U[,1]
    u2 <- U[,2]
    
    # Latent event time T
    T <- H0.inv(exp(inv.Lambda(u1) - X.noint %*% beta.true))
    
    # Latent censoring time C
    lambda.c <- 0.10
    C <- (-1/lambda.c) * log(1 - u2)
    
  } else if (DGP %% 20 == 4) { # For AFT_ll: Independence, C ~ Exp, ~65% cens.
    
    # Latent event time T
    T <- H0.inv(exp(inv.Lambda(runif(n)) - X.noint %*% beta.true))
    
    # Latent censoring time C
    lambda.c <- 1
    C <- (-1/lambda.c) * log(1 - runif(n))
    
  } else if (DGP %% 20 == 5) { # For AFT_ll: Pos. dep., C ~ Exp, ~65% cens.
    
    # Specify copula
    cop_to_use <- frankCopula(param = 6, dim = 2)
    
    # Generate data from copula
    U <- rCopula(n, cop_to_use)
    u1 <- U[,1]
    u2 <- U[,2]
    
    # Latent event time T
    T <- H0.inv(exp(inv.Lambda(u1) - X.noint %*% beta.true))
    
    # Latent censoring time C
    lambda.c <- 0.6
    C <- (-1/lambda.c) * log(1 - u2)
    
  } else if (DGP %% 20 == 6) { # For AFT_ll: Neg. dep., C ~ Exp, ~65% cens.
    
    # Specify copula
    cop_to_use <- frankCopula(param = -6, dim = 2)
    
    # Generate data from copula
    U <- rCopula(n, cop_to_use)
    u1 <- U[,1]
    u2 <- U[,2]
    
    # Latent event time T
    T <- H0.inv(exp(inv.Lambda(u1) - X.noint %*% beta.true))
    
    # Latent censoring time C
    lambda.c <- 1.3
    C <- (-1/lambda.c) * log(1 - u2)
    
  } else if (DGP %% 20 == 7) { # For Cox_wb: Independence, C ~ Exp, ~25% cens.
    
    # Latent event time T
    T <- H0.inv(exp(inv.Lambda(runif(n)) - X.noint %*% beta.true))
    
    # Latent censoring time C
    lambda.c <- 0.25
    C <- (-1/lambda.c) * log(1 - runif(n))
    
  } else if (DGP %% 20 == 8) { # For Cox_wb: Pos. dep., C ~ Exp, ~25% cens.
    
    # Specify copula
    cop_to_use <- frankCopula(param = 6, dim = 2)
    
    # Generate data from copula
    U <- rCopula(n, cop_to_use)
    u1 <- U[,1]
    u2 <- U[,2]
    
    # Latent event time T
    T <- H0.inv(exp(inv.Lambda(u1) - X.noint %*% beta.true))
    
    # Latent censoring time C
    lambda.c <- 0.31
    C <- (-1/lambda.c) * log(1 - u2)
    
  } else if (DGP %% 20 == 9) { # For Cox_wb: Neg. dep., C ~ Exp, ~25% cens.
    
    # Specify copula
    cop_to_use <- frankCopula(param = -6, dim = 2)
    
    # Generate data from copula
    U <- rCopula(n, cop_to_use)
    u1 <- U[,1]
    u2 <- U[,2]
    
    # Latent event time T
    T <- H0.inv(exp(inv.Lambda(u1) - X.noint %*% beta.true))
    
    # Latent censoring time C
    lambda.c <- 0.19
    C <- (-1/lambda.c) * log(1 - u2)
    
  } else if (DGP %% 20 == 10) { # For Cox_wb: Independence, C ~ Exp, ~65% cens.
    
    # Latent event time T
    T <- H0.inv(exp(inv.Lambda(runif(n)) - X.noint %*% beta.true))
    
    # Latent censoring time C
    lambda.c <- 1.3
    C <- (-1/lambda.c) * log(1 - runif(n))
    
  } else if (DGP %% 20 == 11) { # For Cox_wb: Pos. dep., C ~ Exp, ~65% cens.
    
    # Specify copula
    cop_to_use <- frankCopula(param = 6, dim = 2)
    
    # Generate data from copula
    U <- rCopula(n, cop_to_use)
    u1 <- U[,1]
    u2 <- U[,2]
    
    # Latent event time T
    T <- H0.inv(exp(inv.Lambda(u1) - X.noint %*% beta.true))
    
    # Latent censoring time C
    lambda.c <- 0.9
    C <- (-1/lambda.c) * log(1 - u2)
    
  } else if (DGP %% 20 == 12) { # For Cox_wb: Neg. dep., C ~ Exp, ~65% cens.
    
    # Specify copula
    cop_to_use <- frankCopula(param = -6, dim = 2)
    
    # Generate data from copula
    U <- rCopula(n, cop_to_use)
    u1 <- U[,1]
    u2 <- U[,2]
    
    # Latent event time T
    T <- H0.inv(exp(inv.Lambda(u1) - X.noint %*% beta.true))
    
    # Latent censoring time C
    lambda.c <- 1.5
    C <- (-1/lambda.c) * log(1 - u2)
    
  }  else if (DGP %% 20 == 13) { # For AFT_ll: Independence, C ~ Exp, ~2% cens.
    
    # Latent event time T
    T <- H0.inv(exp(inv.Lambda(runif(n)) - X.noint %*% beta.true))
    
    # Latent censoring time C
    lambda.c <- 0.002
    C <- (-1/lambda.c) * log(1 - runif(n))
    
  } else if (DGP %% 20 == 14) { # For AFT_ll: Pos. dep., C ~ Exp, ~2% cens.
    
    # Specify copula
    cop_to_use <- frankCopula(param = 6, dim = 2)
    
    # Generate data from copula
    U <- rCopula(n, cop_to_use)
    u1 <- U[,1]
    u2 <- U[,2]
    
    # Latent event time T
    T <- H0.inv(exp(inv.Lambda(u1) - X.noint %*% beta.true))
    
    # Latent censoring time C
    lambda.c <- 0.01
    C <- (-1/lambda.c) * log(1 - u2)
    
  } else if (DGP %% 20 == 15) { # For AFT_ll: Neg. dep., C ~ Exp, ~2% cens.
    
    # Specify copula
    cop_to_use <- frankCopula(param = -6, dim = 2)
    
    # Generate data from copula
    U <- rCopula(n, cop_to_use)
    u1 <- U[,1]
    u2 <- U[,2]
    
    # Latent event time T
    T <- H0.inv(exp(inv.Lambda(u1) - X.noint %*% beta.true))
    
    # Latent censoring time C
    lambda.c <- 0.04
    C <- (-1/lambda.c) * log(1 - u2)
    
  } else if (DGP %% 20 == 16) { # For Cox_wb: Independence, C ~ Exp, ~2% cens.
    
    # Latent event time T
    T <- H0.inv(exp(inv.Lambda(runif(n)) - X.noint %*% beta.true))
    
    # Latent censoring time C
    lambda.c <- 0.01
    C <- (-1/lambda.c) * log(1 - runif(n))
    
  } else if (DGP %% 20 == 17) { # For Cox_wb: Pos. dep., C ~ Exp, ~2% cens.
    
    # Specify copula
    cop_to_use <- frankCopula(param = 6, dim = 2)
    
    # Generate data from copula
    U <- rCopula(n, cop_to_use)
    u1 <- U[,1]
    u2 <- U[,2]
    
    # Latent event time T
    T <- H0.inv(exp(inv.Lambda(u1) - X.noint %*% beta.true))
    
    # Latent censoring time C
    lambda.c <- 0.03
    C <- (-1/lambda.c) * log(1 - u2)
    
  } else if (DGP %% 20 == 18) { # For Cox_wb: Neg. dep., C ~ Exp, ~2% cens.
    
    # Specify copula
    cop_to_use <- frankCopula(param = -6, dim = 2)
    
    # Generate data from copula
    U <- rCopula(n, cop_to_use)
    u1 <- U[,1]
    u2 <- U[,2]
    
    # Latent event time T
    T <- H0.inv(exp(inv.Lambda(u1) - X.noint %*% beta.true))
    
    # Latent censoring time C
    lambda.c <- 0.005
    C <- (-1/lambda.c) * log(1 - u2)
    
  }
  
  # Observed identified minimum
  Y <- pmin(T, C)
  Delta <- as.numeric(Y == T)
  
  # Collect all variables
  data <- data.frame(Y, Delta, X)
  colnames(data) <- c("Y", "Delta", colnames(X))
  
  # Histogram of the observed times
  if (plot.data) {
    print(sprintf("Percentage of censored observations: %.2f%%",
                  100*(1 - sum(data$Delta)/n)))
    hist(Y)
  }
  
  # Return the results
  data
}

#' @title Data generation function for the simulation regarding misspecification.
#' 
#' @description This function generates a data set according to the specified
#' arguments. It is mostly a copy-paste from the function generateData_simMain.R
#' above, but the censoring distribution slightly addapted in order to control
#' the percentage of censored observations.
#' 
#' @param beta.true True covariate vector, as a function of time.
#' @param n Sample size.
#' @param n.cov Number of covariates.
#' @param options.data.gen List of additional arguments.
#' @param H0.inv Inverse of the intercept (function of time).
#' @param plot.data Boolean value indicating whether or not to plot the
#' generated data set. Default value is \code{plot.data = FALSE}.
#' 
#' @import stats copula
#' 
#' @noRd
#' 
generateData_simMiss <- function(beta.true, n, n.cov, options.data.gen, H0.inv,
                                 plot.data = FALSE) {
  
  # Parameter used for development of this function. Set to TRUE MANUALLY if
  # desired.
  # NOTE: When the function is run as is with 'test.mode = TRUE', it will fail!
  test.mode <- FALSE
  
  # Extract the necessary hyperparameters
  if (options.data.gen[["link.function"]] == "AFT_ll") {
    inv.Lambda <- Lambda_inverse_AFT_ll
  } else if (options.data.gen[["link.function"]] == "Cox_wb") {
    inv.Lambda <- Lambda_inverse_Cox_wb
  }
  DGP <- options.data.gen[["DGP"]]
  
  # Subset the parameter vector to the first n.cov covariate effects.
  if (class(beta.true) == "function") {
    beta.true <- beta.true(0)[2:(n.cov + 1)]
  }
  
  # Set the types of the covariates to be used based on the selected DGP.
  #
  # NOTE TO SELF: If you change these thresholds (20, 40, 60), also change them
  #               accordingly in 'simulate1D.R', 'set.hyperparameters.R', 
  #               'simFuncWrapper.R' and 'simulate1D.CCDF.R'
  if (DGP <= 20) {
    type.cov <- rep("c", n.cov)
  } else if (20 < DGP & DGP <= 40) {
    type.cov <- rep("c", n.cov)
    if ((ceiling(n.cov/2) + 1) <= length(type.cov)) {
      type.cov[(ceiling(n.cov/2) + 1):length(type.cov)] <- "b"
    }
  } else {
    type.cov <- rep("b", n.cov)
  }
  
  # Generate the intercept and the covariates
  X <- rep(1, n)
  for (i in 1:n.cov) {
    if (type.cov[i] == "c") {
      X <- cbind(X, rnorm(n))
    } else if (type.cov[i] == "b") {
      X <- cbind(X, rbinom(n, 1, 0.5))
    } else {
      stop("Invalid value for type of covariate")
    }
  }
  colnames(X) <- paste0("X", 0:n.cov)
  
  # Get matrix of just the covariates (no intercept)
  X.noint <- X[, -1]
  
  # Set dependence structure based on selected DGP
  if (DGP %% 20 == 1) { # For AFT_ll: Independence, C ~ Exp, ~25% cens.
    
    # Latent event time T
    T <- H0.inv(exp(inv.Lambda(runif(n)) - X.noint %*% beta.true))
    
    # Latent censoring time C
    lambda.c <- 0.12
    C <- (-1/lambda.c) * log(1 - runif(n))
    
    if (test.mode) {
      get.average.censoring(1, DGPset = "Miss")
    }
    
  } else if (DGP %% 20 == 2) { # For AFT_ll: Pos. dep., C ~ Unif, ~25% cens.
    
    # Specify copula
    cop_to_use <- frankCopula(param = 6, dim = 2)
    
    # Generate data from copula
    U <- rCopula(n, cop_to_use)
    u1 <- U[,1]
    u2 <- U[,2]
    
    # Latent event time T
    T <- H0.inv(exp(inv.Lambda(u1) - X.noint %*% beta.true))
    
    # Latent censoring time C
    lambda.c <- 0.18
    C <- (-1/lambda.c) * log(1 - u2)
    
    if (test.mode) {
      get.average.censoring(2, DGPset = "Miss")
    }
    
  } else if (DGP %% 20 == 4) { # For AFT_ll: Independence, C ~ Exp, ~65% cens.
    
    # Latent event time T
    T <- H0.inv(exp(inv.Lambda(runif(n)) - X.noint %*% beta.true))
    
    # Latent censoring time C
    lambda.c <- 1
    C <- (-1/lambda.c) * log(1 - runif(n))
    
    if (test.mode) {
      get.average.censoring(4, DGPset = "Miss")
    }
    
  } else if (DGP %% 20 == 5) { # For AFT_ll: Pos. dep., C ~ Exp, ~65% cens.
    
    # Specify copula
    cop_to_use <- frankCopula(param = 6, dim = 2)
    
    # Generate data from copula
    U <- rCopula(n, cop_to_use)
    u1 <- U[,1]
    u2 <- U[,2]
    
    # Latent event time T
    T <- H0.inv(exp(inv.Lambda(u1) - X.noint %*% beta.true))
    
    # Latent censoring time C
    lambda.c <- 0.75
    C <- (-1/lambda.c) * log(1 - u2)
    
    if (test.mode) {
      get.average.censoring(5, DGPset = "Miss")
    }
    
  }
  
  # Observed identified minimum
  Y <- pmin(T, C)
  Delta <- as.numeric(Y == T)
  
  # Collect all variables
  data <- data.frame(Y, Delta, X)
  colnames(data) <- c("Y", "Delta", colnames(X))
  
  # Histogram of the observed times
  if (plot.data) {
    print(sprintf("Percentage of censored observations: %.2f%%",
                  100*(1 - sum(data$Delta)/n)))
    hist(Y)
  }
  
  # Return the results
  data
}

#' @title Data generation function for the simulations under many covariates.
#' 
#' @description This function generates a data set according to the specified
#' arguments.
#' 
#' @param beta.true True covariate vector, as a function of time.
#' @param n Sample size.
#' @param n.cov Number of covariates.
#' @param options List of additional arguments.
#' @param H0.inv Inverse of the intercept (function of time).
#' @param plot.data Boolean value indicating whether or not to plot the
#' generated data set. Default value is \code{plot.data = FALSE}.
#' 
#' @import stats copula
#' 
#' @noRd
#' 
generateData_simManyCov <- function(beta.true, n, n.cov, options, H0.inv,
                                    plot.data = FALSE) {
  
  # Extract the necessary hyperparameters
  if (options[["link.function"]] == "AFT_ll") {
    inv.Lambda <- Lambda_inverse_AFT_ll
  } else if (options[["link.function"]] == "Cox_wb") {
    inv.Lambda <- Lambda_inverse_Cox_wb
  }
  DGP <- options[["DGP"]]
  
  # Subset the parameter vector to the first n.cov covariate effects.
  if (class(beta.true) == "function") {
    beta.true <- beta.true(0)[2:(n.cov + 1)]
  }
  
  # Set the types of the covariates to be used based on the selected DGP.
  #
  # NOTE TO SELF: If you change these thresholds (20, 40, 60), also change them
  #               accordingly in 'simulate1D.R', 'set.hyperparameters.R', 
  #               'simFuncWrapper.R' and 'simulate1D.CCDF.R'
  if (DGP <= 20) {
    type.cov <- rep("c", n.cov)
  } else if (20 < DGP & DGP <= 40) {
    type.cov <- rep("c", n.cov)
    if ((ceiling(n.cov/2) + 1) <= length(type.cov)) {
      type.cov[(ceiling(n.cov/2) + 1):length(type.cov)] <- "b"
    }
  } else {
    type.cov <- rep("b", n.cov)
  }
  
  # Generate the intercept and the covariates
  X <- rep(1, n)
  for (i in 1:n.cov) {
    if (type.cov[i] == "c") {
      X <- cbind(X, rnorm(n))
    } else if (type.cov[i] == "b") {
      X <- cbind(X, rbinom(n, 1, 0.5))
    } else {
      stop("Invalid value for type of covariate")
    }
  }
  colnames(X) <- paste0("X", 0:n.cov)
  
  # Get matrix of just the covariates (no intercept)
  X.noint <- X[, -1]
  
  # Set dependence structure based on selected DGP
  if (DGP %% 20 == 1) { # For AFT_ll: Independence, C ~ Exp, ~25% cens.
    
    # Latent event time T
    T <- H0.inv(exp(inv.Lambda(runif(n)) - X.noint %*% beta.true))
    
    # Latent censoring time C
    lambda.c <- 0.35
    C <- (-1/lambda.c) * log(1 - runif(n))
    
  } else if (DGP %% 20 == 2) { # For AFT_ll: Pos. dep., C ~ Exp, ~25% cens.
    
    # Specify copula
    cop_to_use <- frankCopula(param = 6, dim = 2)
    
    # Generate data from copula
    U <- rCopula(n, cop_to_use)
    u1 <- U[,1]
    u2 <- U[,2]
    
    # Latent event time T
    T <- H0.inv(exp(inv.Lambda(u1) - X.noint %*% beta.true))
    
    # Latent censoring time C
    lambda.c <- 0.45
    C <- (-1/lambda.c) * log(1 - u2)
    
  } else if (DGP %% 20 == 3) { # For AFT_ll: Neg. dep., C ~ Exp, ~25% cens.
    
    # Specify copula
    cop_to_use <- frankCopula(param = -6, dim = 2)
    
    # Generate data from copula
    U <- rCopula(n, cop_to_use)
    u1 <- U[,1]
    u2 <- U[,2]
    
    # Latent event time T
    T <- H0.inv(exp(inv.Lambda(u1) - X.noint %*% beta.true))
    
    # Latent censoring time C
    lambda.c <- 0.25
    C <- (-1/lambda.c) * log(1 - u2)
    
  } else if (DGP %% 20 == 4) { # For AFT_ll: Independence, C ~ Exp, ~65% cens.
    
    # Latent event time T
    T <- H0.inv(exp(inv.Lambda(runif(n)) - X.noint %*% beta.true))
    
    # Latent censoring time C
    lambda.c <- 5
    C <- (-1/lambda.c) * log(1 - runif(n))
    
  } else if (DGP %% 20 == 5) { # For AFT_ll: Pos. dep., C ~ Exp, ~65% cens.
    
    # Specify copula
    cop_to_use <- frankCopula(param = 6, dim = 2)
    
    # Generate data from copula
    U <- rCopula(n, cop_to_use)
    u1 <- U[,1]
    u2 <- U[,2]
    
    # Latent event time T
    T <- H0.inv(exp(inv.Lambda(u1) - X.noint %*% beta.true))
    
    # Latent censoring time C
    lambda.c <- 3.5
    C <- (-1/lambda.c) * log(1 - u2)
    
  } else if (DGP %% 20 == 6) { # For AFT_ll: Neg. dep., C ~ Exp, ~65% cens.
    
    # Specify copula
    cop_to_use <- frankCopula(param = -6, dim = 2)
    
    # Generate data from copula
    U <- rCopula(n, cop_to_use)
    u1 <- U[,1]
    u2 <- U[,2]
    
    # Latent event time T
    T <- H0.inv(exp(inv.Lambda(u1) - X.noint %*% beta.true))
    
    # Latent censoring time C
    lambda.c <- 6
    C <- (-1/lambda.c) * log(1 - u2)
    
  } else if (DGP %% 20 == 7) { # For Cox_wb: Independence, C ~ Exp, ~25% cens.
    
    # Latent event time T
    T <- H0.inv(exp(inv.Lambda(runif(n)) - X.noint %*% beta.true))
    
    # Latent censoring time C
    lambda.c <- 0.7
    C <- (-1/lambda.c) * log(1 - runif(n))
    
  } else if (DGP %% 20 == 8) { # For Cox_wb: Pos. dep., C ~ Exp, ~25% cens.
    
    # Specify copula
    cop_to_use <- frankCopula(param = 6, dim = 2)
    
    # Generate data from copula
    U <- rCopula(n, cop_to_use)
    u1 <- U[,1]
    u2 <- U[,2]
    
    # Latent event time T
    T <- H0.inv(exp(inv.Lambda(u1) - X.noint %*% beta.true))
    
    # Latent censoring time C
    lambda.c <- 0.9
    C <- (-1/lambda.c) * log(1 - u2)
    
  } else if (DGP %% 20 == 9) { # For Cox_wb: Neg. dep., C ~ Exp, ~25% cens.
    
    # Specify copula
    cop_to_use <- frankCopula(param = -6, dim = 2)
    
    # Generate data from copula
    U <- rCopula(n, cop_to_use)
    u1 <- U[,1]
    u2 <- U[,2]
    
    # Latent event time T
    T <- H0.inv(exp(inv.Lambda(u1) - X.noint %*% beta.true))
    
    # Latent censoring time C
    lambda.c <- 0.55
    C <- (-1/lambda.c) * log(1 - u2)
    
  } else if (DGP %% 20 == 10) { # For Cox_wb: Independence, C ~ Exp, ~65% cens.
    
    # Latent event time T
    T <- H0.inv(exp(inv.Lambda(runif(n)) - X.noint %*% beta.true))
    
    # Latent censoring time C
    lambda.c <- 8
    C <- (-1/lambda.c) * log(1 - runif(n))
    
  } else if (DGP %% 20 == 11) { # For Cox_wb: Pos. dep., C ~ Exp, ~65% cens.
    
    # Specify copula
    cop_to_use <- frankCopula(param = 6, dim = 2)
    
    # Generate data from copula
    U <- rCopula(n, cop_to_use)
    u1 <- U[,1]
    u2 <- U[,2]
    
    # Latent event time T
    T <- H0.inv(exp(inv.Lambda(u1) - X.noint %*% beta.true))
    
    # Latent censoring time C
    lambda.c <- 6
    C <- (-1/lambda.c) * log(1 - u2)
    
  } else if (DGP %% 20 == 12) { # For Cox_wb: Neg. dep., C ~ Exp, ~65% cens.
    
    # Specify copula
    cop_to_use <- frankCopula(param = -6, dim = 2)
    
    # Generate data from copula
    U <- rCopula(n, cop_to_use)
    u1 <- U[,1]
    u2 <- U[,2]
    
    # Latent event time T
    T <- H0.inv(exp(inv.Lambda(u1) - X.noint %*% beta.true))
    
    # Latent censoring time C
    lambda.c <- 9
    C <- (-1/lambda.c) * log(1 - u2)
    
  }
  
  # Observed identified minimum
  Y <- pmin(T, C)
  Delta <- as.numeric(Y == T)
  
  # Collect all variables
  data <- data.frame(Y, Delta, X)
  colnames(data) <- c("Y", "Delta", colnames(X))
  
  # Histogram of the observed times
  if (plot.data) {
    print(sprintf("Percentage of censored observations: %.2f%%",
                  100*(1 - sum(data$Delta)/n)))
    hist(Y)
  }
  
  # Return the results
  data
}

#### Link functions ####

#' @title Link function (Cox model)
#' 
#' @description
#' This function defines the Cox PH link function.
#' 
#' @param t time parameter.
#' 
Lambda_Cox_wb <- function(t) {
  1 - exp(-exp(t))
}

#' @title Derivative of link function (Cox model)
#' 
#' @description
#' This function defines the derivative of the Cox PH link function.
#' 
#' @param t time parameter.
#' 
dLambda_Cox_wb <- function(t) {
  exp(t - exp(t))
}

#' @title Inverse of link function (Cox model)
#' 
#' @description
#' This function defines the inverse of the Cox PH link function.
#' 
#' @param p probability.
#' 
Lambda_inverse_Cox_wb <- function(p) {
  log(-log(1-p))
}

#' @title Link function (AFT model)
#' 
#' @description
#' This function defines the AFT link function.
#' 
#' @param t time parameter.
#' 
Lambda_AFT_ll <- function(t) {
  1 - 1/(1 + exp(t))
}

#' @title Derivative of link function (AFT model)
#' 
#' @description
#' This function defines the derivative of the AFT link function.
#' 
#' @param t time parameter.
#' 
dLambda_AFT_ll <- function(t) {
  exp(t)/(1 + exp(t))^2
}

#' @title Inverse of link function (AFT model)
#' 
#' @description
#' This function defines the inverse of the AFT link function.
#' 
#' @param p probability.
#' 
Lambda_inverse_AFT_ll <- function(p) {
  log(p / (1 - p))
}

#### Family of instrumental functions ####

#' @title Normalize the covariates of a data set to lie in the unit interval by
#' scaling based on the ranges of the covariates.
#' 
#' @description This function normalized the covariates in the data to lie in 
#' the unit interval based on either the empirical or known ranges of the
#' covariates. It is useful to perform this step when defining the instrumental
#' functions later on. This function is used in \code{G.box.R}, \code{G.spline.R}
#' and by extension in \code{G.cd.R}.
#'
#' @param data (optional) Data set to be used to construct the normalizing
#' transformation. Default is \code{data = NULL}.
#' @param x (optional) Vector of covariates to be normalized alongside the data.
#' Default is \code{x = NULL}.
#' @param cov.ranges (optional) Matrix that specifies the range of each of the
#' covariates in the data set. Each column corresponds to a covariate. The first
#' row contains the lower bound, the second row contains the upper bound.
#' If not supplied, the data will be normalized based on the minimum and maximum
#' detected values. If supplied, the non data-dependent transformation function
#' listed in the appendix of Andrews, Shi 2013 will be used. Default is
#' \code{cov.ranges = NULL}.
#' @param norm.cov.out (optional) The output of a previous call to this function.
#' Can be used to speed up computation. If both \code{data} and
#' \code{norm.cov.out} are supplied to the function, this method will throw an
#' error. Default is \code{norm.cov.out = NULL}.
#' @param idxs.c (optional) Vector of indices of covariates that are continuous.
#' Note that that indices are relative to the covariate vector, not the full
#' data set. Default value is \code{idxs.c = "all"}, which indicates that all
#' elements should be regarded as continuous. If \code{idxs.c = NULL}, all
#' elements are regarded as discrete.
#' @param ... Allows easier interchangeability between covariate normalization
#' functions. All arguments specified under \code{...} will be ignored.
#' 
#' @references Andrews, D.W.K. and Shi, X. (2013). Inference based on
#' confitional moment inequalities. Econometrica. 81(2):609-666.
#' 
normalize.covariates <- function(data = NULL, x = NULL, cov.ranges = NULL,
                                 idxs.c = "all", norm.cov.out = NULL, ...) {
  
  # Precondition checks
  if (is.null(data) & is.null(norm.cov.out)) {
    stop("Either data or norm.cov.out should be supplied to this function.")
  }
  if (!is.null(data) & !is.null(norm.cov.out)) {
    stop("Ambiguous function arguments: both data and norm.cov.out are supplied.")
  }
  
  # Extract the covariates from the data set, if applicable. Else extract the
  # necessary parameters from the previous function call.
  if (!is.null(data)) {
    
    # Get all covariates names
    cov.idxs <- which(grepl("X[[:digit:]]+", colnames(data)))[-1]
    covariate.names <- colnames(data)[cov.idxs]
    
    # Only retain names of continuous covariates
    idxs.c <- if (all(idxs.c == "all")) {1:length(covariate.names)} else {idxs.c}
    cont.cov.names <- covariate.names[idxs.c]
  } else {
    covariate.names <- norm.cov.out$covariate.names
    cont.cov.names <- norm.cov.out$cont.cov.names
  }
  
  # If supplied, rename the entries of x.
  if (!is.null(x)) {
    names(x) <- covariate.names
  }
  
  # Initialize object that will store the data with normalized covariates
  normalized.data <- data
  x.c.norm <- NULL
  
  # Compute the minimum and maximum value of each continuous covariates based on
  # the data. If the output of a previous call to this function was provided,
  # this step can be skipped.
  X.ranges <- matrix(nrow = 2, ncol = length(cont.cov.names))
  colnames(X.ranges) <- cont.cov.names
  if (is.null(norm.cov.out)) {
    for (cov.name in cont.cov.names) {
      if (is.null(cov.ranges)) {
        X.ranges[1, cov.name] <- min(data[, cov.name])
        X.ranges[2, cov.name] <- max(data[, cov.name])
      } else {
        X.ranges[1, cov.name] <- cov.ranges[1, cov.name]
        X.ranges[2, cov.name] <- cov.ranges[2, cov.name]
      }
    }
  } else {
    X.ranges <- norm.cov.out$X.ranges
  }
  
  # For each covariate...
  for (cov.name in cont.cov.names) {
    
    # Get covariate of this iteration, alongside its min/max value
    if (!is.null(data)) {
      X <- data[, cov.name]
    }
    min.X <- X.ranges[1, cov.name]
    max.X <- X.ranges[2, cov.name]
    
    # Normalize the covariate based on the functions defined in the supplement
    # of Andrews, Shi (2013).
    if ((min.X > -Inf) & (max.X < Inf)) {
      if (is.null(norm.cov.out)) {
        X.normalized <- (X - min.X)/(max.X - min.X)
      }
      if (!is.null(x)) {
        x.elem <- x[cov.name]
        x.c.norm <- c(x.c.norm, (x.elem - min.X)/(max.X - min.X))
      }
    } else if ((min.X > -Inf) & (max.X == Inf)) {
      if (is.null(norm.cov.out)) {
        X.normalized <- (exp(X - min.X) - 1)/(1 + exp(X - min.X))
      }
      if (!is.null(x)) {
        x.elem <- x[cov.name]
        x.c.norm <- c(x.c.norm, (exp(x.elem - min.X) - 1)/(1 + exp(x.elem - min.X)))
      }
    } else if ((min.X == -Inf) & max.X < Inf) {
      if (is.null(norm.cov.out)) {
        X.normalized <- (2*exp(X - max.X))/(1 + exp(X - max.X))
      }
      if (!is.null(x)) {
        x.elem <- x[cov.name]
        x.c.norm <- c(x.c.norm, (2*exp(x.elem - max.X))/(1 + exp(x.elem - max.X)))
      }
    } else {
      if (is.null(norm.cov.out)) {
        X.normalized <- exp(X)/(1 - exp(X))
      }
      if (!is.null(x)) {
        x.elem <- x[cov.name]
        x.c.norm <- c(x.c.norm, exp(x.elem)/(1 - exp(x.elem)))
      }
    }
    
    # Transform x.norm to a named vector (in the first iteration, it will be a
    # list).
    x.c.norm <- unlist(x.c.norm)
    
    # Store the result
    if (is.null(norm.cov.out)) {
      normalized.data[, cov.name] <- X.normalized
    }
  }
  
  # If a previous output of this function was supplied, access the precomputed
  # normalized data
  if (!is.null(norm.cov.out)) {
    normalized.data <- norm.cov.out$normalized.data
  }
  
  # If x fell outside of the range of the covariates, the normalized values for 
  # x might be negative. Note however that during the execution of the test by
  # Bei (2024), this should never occur.
  x.c.norm <- pmax(pmin(x.c.norm, 1), 0)
  
  # Reconstruct the full vector for x, with normalized continuous elements
  x.norm <- x
  x.norm[cont.cov.names] <- x.c.norm
  
  # Return the normalized data
  return(list("normalized.x" = x.norm,
              "normalized.data" = normalized.data,
              "X.ranges" = X.ranges,
              "covariate.names" = covariate.names,
              "cont.cov.names" = cont.cov.names))
}

#' @title Normalize the covariates of a data set to lie in the unit interval by
#' transforming based on PCA.
#' 
#' @description This function normalized the covariates in the data to lie in 
#' the unit interval based on a principal component analysis. It is useful to
#' perform this step when defining the instrumental functions later on. This
#' function is used in \code{G.box}, \code{G.spline} and by extension \code{G.cd}.
#'
#' @param data (optional) Data set to be used to construct the normalizing
#' transformation. Default is \code{data = NULL}.
#' @param x (optional) Vector of covariates to be normalized alongside the data.
#' Default is \code{x = NULL}.
#' @param idxs.c (optional) Vector of indices of covariates that are continuous.
#' Note that that indices are relative to the covariate vector, not the full
#' data set. Default value is \code{idxs.c = "all"}, which indicates that all
#' elements should be regarded as continuous. If \code{idxs.c = NULL}, all
#' elements are regarded as discrete.
#' @param norm.cov.out (optional) The output of a previous call to this function.
#' Can be used to speed up computation. If both \code{data} and
#' \code{norm.cov.out} are supplied to the function, the function will throw an
#' error. Default is \code{norm.cov.out = NULL}
#' @param ... Allows easier interchangeability between covariate normalization
#' functions. All arguments specified under \code{...} will be ignored.
#' 
#' @import stats
#' 
normalize.covariates2 <- function(data = NULL, x = NULL, idxs.c = "all",
                                  norm.cov.out = NULL, ...) {
  
  # Precondition checks
  if (is.null(data) & is.null(norm.cov.out)) {
    stop("Either data or norm.cov.out should be supplied to this function.")
  }
  if (!is.null(data) & !is.null(norm.cov.out)) {
    stop("Ambiguous function arguments: both data and norm.cov.out are supplied.")
  }
  
  # Extract the covariates from the data set, if applicable. Else extract the
  # necessary parameters from the previous function call.
  if (!is.null(data)) {
    
    # Get all covariates names
    cov.idxs <- which(grepl("X[[:digit:]]+", colnames(data)))[-1]
    covariate.names <- colnames(data)[cov.idxs]
    
    # Only retain names of continuous covariates
    idxs.c <- if (class(idxs.c) == "character") {1:length(covariate.names)} else {idxs.c}
    cont.cov.names <- covariate.names[idxs.c]
    
    # Obtain matrix of continuous covariates
    covariate.matrix <- as.matrix(data[, cont.cov.names, drop = FALSE])
    
  } else {
    covariate.names <- norm.cov.out$covariate.names
    cont.cov.names <- norm.cov.out$cont.cov.names
  }
  
  # If supplied, rename the entries of x.
  if (!is.null(x)) {
    names(x) <- covariate.names
  }
  
  # Define some useful variables
  n <- nrow(data)
  
  # Define function to transform unit circle towards [-1, 1]^d. Used to make
  # observations more evenly distributed in [-1, 1]^d later on. Else, extract
  # them from the previous output.
  transform.scores <- function(scores) {
    sin((1/2)*pi*scores)
  }
  
  # If the output of a previous call to this function was not supplied, obtain
  # the transformation parameters to be used.
  if (is.null(norm.cov.out)) {
    
    # If there are continuous covariates in the data set to be normalized...
    if (length(idxs.c > 0)) {
      
      # Apply principal component analysis to the covariates and obtain scores
      ev <- eigen(var(covariate.matrix))$vectors
      scores <- covariate.matrix %*% ev
      
      # Scale every covariate to have range of length 2
      scale <- apply(scores, 2, max) - apply(scores, 2, min)
      scale.mat <- matrix(rep(2/scale, n), ncol = ncol(scores), byrow = TRUE)
      scores.scaled <- scale.mat * scores
      
      # Shift the covariates into [-1, 1]^d
      shift <- apply(scores.scaled, 2, min) + 1
      shift.mat <- matrix(rep(shift, n), ncol = ncol(scores.scaled), byrow = TRUE)
      scores.scaled.shifted <- scores.scaled - shift.mat
      
      # Transform covariates to be more uniformly distributed in [-1, 1]^d
      scores.transformed <- transform.scores(scores.scaled.shifted)
      
      # Shift and scale the covariates into [0, 1]^d
      covariates.norm <- (1/2)*scores.transformed + (1/2)
      
      # Create the data frame with normalized covariates
      data.norm <- data
      data.norm[, cont.cov.names] <- covariates.norm
      
    # If there are no continuous covariates in the data to be normalized...
    } else {
      
      ev <- NULL
      scale <- NULL
      shift <- NULL
      data.norm <- data
    }
    
  # If the argument for norm.cov.out was provided...
  } else {
    
    ev <- norm.cov.out$ev
    scale <- norm.cov.out$scale
    shift <- norm.cov.out$shift
    data.norm <- norm.cov.out$normalized.data
    
  }
  
  # If x was supplied, transform x in the same way
  x.norm <- NULL
  if (!is.null(x)) {
    
    if (length(idxs.c) > 0) {
      # Transform the continuous elements in x
      x.c <- x[cont.cov.names]
      x.c.norm <- (1/2)*transform.scores(((matrix(x.c, nrow = 1) %*% ev) * (2/scale)) - shift) + (1/2)
      x.c.norm <- as.numeric(x.c.norm)
      
      # Construct entire vector
      x.norm <- x
      x.norm[cont.cov.names] <- x.c.norm
    } else {
      x.norm <- x
    }
  }
  
  # Return the results
  return(list("normalized.x" = x.norm,
              "normalized.data" = data.norm,
              "ev" = ev,
              "scale" = scale,
              "shift" = shift,
              "covariate.names" = covariate.names,
              "cont.cov.names" = cont.cov.names))
}

#' @title Get anchor points on which to base the instrumental functions 
#' 
#' @description The points returned by this function can be used as corner
#' points in the family of box functions, or as knots in the family of B-spline
#' functions.
#'
#' @param data Data set.
#' @param n.if.per.cov Number of instrumental functions to use per continuous
#' covariate.
#' @param normalized Boolean value indicating whether the covariates in the
#' given data frame have been normalized. Default is \code{normalized = FALSE}.
#' 
get.anchor.points <- function(data, n.if.per.cov, normalized = FALSE) {
  
  # Get column indices of covariates in the data (excluding the intercept)
  cov.idxs <- which(grepl("X[[:digit:]]+", colnames(data)))[-1]
  
  # Get the number of covariates
  n.cov <- length(cov.idxs)
  
  # If the data should have normalized covariates, check that it is the case
  if (normalized) {
    for (cov.idx in cov.idxs) {
      if (!all((0 <= data[, cov.idx]) & (data[, cov.idx] <= 1))) {
        stop("Unnormalized data detected")
      }
    }
  }
  
  # Initialize object that will store the anchor points
  ap <- matrix(nrow = n.cov, ncol = n.if.per.cov + 1)
  rownames(ap) <- paste0("X", 1:n.cov)
  
  # For each covariate, determine an appropriate range for the boxes
  for (idx in 1:n.cov) {
    if (normalized) {
      
      # When the data is normalized, anchor points are evenly spaced in [0, 1]
      ap[idx, ] <- seq(0, 1, length.out = n.if.per.cov + 1)
      
    } else {
      
      # Select covariate of this iteration
      X <- data[, cov.idxs[idx]]
      
      # Determine the grid of corner points corresponding to this covariate for the
      # boxes.
      ap[idx, ] <- seq(min(X), max(X), length.out = n.if.per.cov + 1)
      
    }
  }
  
  # Return the results
  ap
}

#' @title Family of box functions
#' 
#' @description This function defined the class of box functions as defined in
#' Willems et al. (2024+).
#' 
#' @param x Vector of covariates to be normalized alongside the data. Default is
#' \code{x = NULL}.
#' @param g.idx Index of the instrumental function, in \{1, ..., n.inst.func\}.
#' @param data Data frame.
#' @param n.box.per.cov Number of box functions to consider per continuous 
#' covariate.
#' @param norm.func Function to be used to normalize the covariates.
#' @param cov.ranges Matrix of ranges of the covariates. Used for normalizing
#' the data to the unit interval before applying the instrumental functions.
#' Default is \code{cov.ranges = NULL}.
#' @param norm.cov.out Output of a preliminary call to the supplied covariate
#' normalization function.
#' @param ... Additional arguments will be ignored. Useful for allowing
#' compatibility with the implementations of other instrument function families.
#' Specifically, it allows to ignore the \code{degree} argument used in
#' 'G.spline.R' and 'G.cd.R'.
#' 
#' @importFrom EnvStats
#' 
G.box <- function(x, g.idx, data, n.box.per.cov, norm.func, cov.ranges = NULL, 
                  norm.cov.out = NULL, ...) {
  
  # Normalize the covariates to lie in the unit interval
  if (is.null(norm.cov.out)) {
    out.norm <- norm.func(data = data, x = x, cov.ranges = cov.ranges,
                          norm.cov.out = NULL)
  } else {
    out.norm <- norm.func(data = NULL, x = x, cov.ranges = cov.ranges,
                          norm.cov.out = norm.cov.out)
  }
  x.norm <- out.norm[["normalized.x"]]
  
  # Get column indices of covariates in the data (excluding the intercept)
  cov.idxs <- which(grepl("X[[:digit:]]+", colnames(data)))[-1]
  
  # Get the number of covariates
  n.cov <- length(cov.idxs)
  
  # Construct matrix of corner points of the boxes
  cp <- matrix(rep(seq(0, 1, length.out = n.box.per.cov + 1), n.cov),
               nrow = n.cov, byrow = TRUE)
  rownames(cp) <- paste0("X", 1:n.cov)
  
  # For each dimension, determine the range of the covariate corresponding to
  # this box function
  range.idxs <- base(g.idx - 1, base = n.box.per.cov, num.digits = n.cov) + 1
  ranges <- NULL
  for (idx in 1:n.cov) {
    ranges <- rbind(ranges, c(cp[idx, range.idxs[idx]], cp[idx, range.idxs[idx] + 1]))
  }
  
  # Return the results
  as.numeric(all(ranges[,1] <= x.norm) & all(x.norm <= ranges[,2]))
}

#' @title Evaluate the specified B-spline, defined on the unit interval
#' 
#' @description This function evaluates the specified B-spline defined on the
#' unit interval, when considering \code{n.if.per.cov} B-splines. Currently, the
#' implementation is based on the one in Andrews, Shi 2013 (supplementary
#' materials).
#' 
#' @param x value inside the unit interval at which to evaluate the spline.
#' @param spline.index Index of the spline to evaluate.
#' @param n.if.per.cov Number of B-splines to consider over the unit interval.
#' @param degree Degree of the B-splines. Default is \code{degree = 3}.
#' 
#' @import splines2
#' 
#' @references Andrews, D.W.K. and Shi, X. (2013). Inference based on
#' confitional moment inequalities. Econometrica. 81(2):609-666.
#' 
Bspline.unit.interval <- function(x, spline.index, n.if.per.cov, degree = 3) {
  
  # Precondition checks
  if (n.if.per.cov <= degree) {
    stop("n.if.per.cov must be larger than degree")
  }
  
  # Create vector of (boundary) knots to be used to construct the spline
  knots <- seq(0, 1, length.out = n.if.per.cov + 1 - degree)
  width <- diff(knots)[1]
  if (degree > 1) {
    knots <- c(-((degree - 1):1) * width, knots, 1 + 1:(degree - 1) * width)
  }
  if (degree != 0) {
    boundary.knots <- c(-degree * width, 1 + degree * width)
  } else {
    boundary.knots <- c(-width, 1)
    knots <- knots[-c(length(knots))]
  }
  
  # Obtain all spline function evaluations
  spline.evals <- bSpline(x, knots = knots, Boundary.knots = boundary.knots,
                          degree = degree)
  
  # Subset spline functions to the ones inside the unit interval
  spline.evals <- spline.evals[max(1, degree):(length(spline.evals) - degree)]
  
  # Return the requested spline function
  spline.evals[spline.index]
}

#' @title Family of spline instrumental functions
#' 
#' @description This function normalizes the covariates to lie in the unit
#' interval and then evaluates each B-spline at each observation, multiplying
#' together the results per observation.
#' 
#' @param x The vector of covariates at which to evaluate the B-splines
#' @param g.idx The index of the instrumental function. Note that g.idx ranges
#' between 1 and n.if.per.cov^n.cov, as an instrumental function is the product
#' of the appropriate B-spline evaluation for each element in the covariate
#' vector.
#' @param data Data frame containing the data.
#' @param n.if.per.cov Number of instrumental variables to be used per covariate.
#' @param norm.func Function to be used to normalize the covariates.
#' @param cov.ranges Matrix of ranges of the covariates. Used for normalizing
#' the covariates. If \code{cov.ranges = NULL}, the data will be normalized in a
#' data-dependent way. Default is \code{cov.ranges = NULL}.
#' @param norm.cov.out Output of a preliminary call to the given covariate
#' normalization function. Default is \code{norm.cov.out = NULL}.
#' @param degree Degree of B-splines to use. Default value is \code{degree = 3}.
#' 
G.spline <- function(x, g.idx, data, n.if.per.cov, norm.func, cov.ranges = NULL,
                     norm.cov.out = NULL, degree = 3) {
  
  # Get the number of covariates
  n.cov <- sum(grepl("X[1-9][[:digit:]]*", colnames(data)))
  
  # Normalize the covariates to lie in the unit interval
  if (is.null(norm.cov.out)) {
    out.norm <- norm.func(data = data, x = x, cov.ranges = cov.ranges,
                          norm.cov.out = NULL)
  } else {
    out.norm <- norm.func(data = NULL, x = x, cov.ranges = cov.ranges,
                          norm.cov.out = norm.cov.out)
  }
  x.norm <- out.norm[["normalized.x"]]
  
  # Get the index of the B-spline to be used for each of the covariates
  spline.idxs <- base(g.idx - 1, base = n.if.per.cov, num.digits = n.cov) + 1
  
  # Evaluate each element of the normalized covariate vector on the appropriate
  # B-spline. Multiply all results.
  spline.args <- cbind(x.norm, spline.idxs, n.if.per.cov, degree)
  spline.wrapper <- function(inp) {Bspline.unit.interval(inp[1], inp[2], inp[3], inp[4])}
  spline.evals <- apply(X = spline.args, MARGIN = 1, FUN = spline.wrapper)
  
  # Return the results
  prod(spline.evals)
}

#' @title Family of continuous/discrete instrumental function
#' 
#' @description The function normalizes the continuous covariates to lie in the
#' unit interval and then evaluates the subvector of continuous covariates on
#' the specified family of instrumental function. For the discrete elements,
#' indicator functions are used for each level.
#' 
#' @param x The vector of covariates at which to evaluate the B-splines
#' @param g.idx The index of the instrumental function.
#' @param data Data frame containing the data.
#' @param n.if.per.cov Number of instrumental functions per continuous covariate.
#' @param idxs.c Vector of indices of the continuous elements in the vector of
#' covariates.
#' @param G.c Family of instrumental functions to use for the subvector of
#' continuous covariates.
#' @param norm.func Function to be used to normalize the covariates.
#' @param discrete.covariate.levels Matrix containing as rows all possible
#' 'combined' levels of the discrete covariates. Default is
#' \code{discrete.covariate.levels = NULL}.
#' @param cov.ranges Matrix containing as its rows the lower and upper bounds
#' for each continuous covariate. Default is \code{cov.ranges = NULL}.
#' @param norm.cov.out Output of a preliminary call to a covariate normalization
#' function (defined above). This is used to speed up computations. Note that
#' this argument should only apply to continuous covariates!! Default is
#' \code{norm.cov.out = NULL}.
#' @param degree Degree of the spline functions to be used as instrumental
#' functions for the continuous covariates (if applicable). Default is
#' \code{degree = 3}.
#' 
G.cd <- function(x, g.idx, data, n.if.per.cov, idxs.c, G.c, norm.func,
                 discrete.covariate.levels = NULL, cov.ranges = NULL,
                 norm.cov.out = NULL, degree = 3) {
  
  # Get the number of covariates
  n.cov <- sum(grepl("X[1-9][[:digit:]]*", colnames(data)))
  
  # Obtain subvectors of continuous and discrete covariates
  x.c <- x[idxs.c]
  x.d <- x[setdiff(1:length(x), idxs.c)]
  
  # Subset data to continuous and discrete covariates (also catching the cases
  # where there are no continuous/discrete variables).
  names.cov.c <- setdiff(paste("X", idxs.c, sep = ""), "X")
  names.cov.d <- setdiff(paste("X", setdiff(1:n.cov, idxs.c), sep = ""), "X")
  data.c <- data[, !(colnames(data) %in% names.cov.d)]
  data.d <- data[, !(colnames(data) %in% names.cov.c)]
  
  # Obtain the indices for the discrete instrumental functions
  n.inst.func.d <- max(nrow(unique(data[, names.cov.d, drop = FALSE])), 1)
  g.idx.d <- ((g.idx - 1) %% n.inst.func.d) + 1
  
  # Obtain the indices for the continuous instrumental functions
  n.inst.func.c <- n.if.per.cov^length(x.c)
  g.idx.c <- ((g.idx - 1) %/% n.inst.func.d) + 1
  
  # Subset the matrix of covariate ranges to the continuous covariates
  cov.ranges.c <- cov.ranges[, colnames(cov.ranges) %in% names.cov.c, drop = FALSE]
  
  # Obtain the instrumental function evaluation for the continuous elements
  eval.c <- 1
  if (length(x.c) > 0) {
    if (length(norm.cov.out$covariate.names) == length(norm.cov.out$cont.cov.names)) {
      norm.cov.out.c <- norm.cov.out
    } else {
      stop("Provided norm.cov.out argument should only apply to the continuous cov.")
    }
    eval.c <- G.c(x.c, g.idx.c, data.c, n.if.per.cov, norm.func, degree = degree,
                  cov.ranges = cov.ranges.c, norm.cov.out = norm.cov.out.c)
  }
  
  # Obtain the instrumental function evaluations for the discrete elements
  eval.d <- 1
  if (length(x.d) > 0) {
    eval.d <- as.numeric(all(x.d == discrete.covariate.levels[g.idx.d,]))
  }
  
  # Return the result
  eval.c * eval.d
}

#' @title Family of discrete/continuous instrumental functions, in the case of
#' many covariates.
#' 
#' @description This function defines the family of discrete/continuous
#' instrumental functions in the case of many covariates. It does so by
#' considering a instrumental functions for each pair of entries in the given
#' covariate vector.
#' 
#' @param x The vector of covariates at which to evaluate the B-splines
#' @param g.idx The index of the instrumental function.
#' @param data Data frame containing the data.
#' @param n.if.per.cov Number of instrumental functions per continuous covariate.
#' @param idxs.c Vector of indices of the continuous elements in the vector of
#' covariates.
#' @param G.c Family of instrumental functions to use for the subvector of
#' continuous covariates.
#' @param norm.func Function to be used to normalize the covariates.
#' @param info.manycov Data frame containing some information about the global
#' structure of the instrumental functions of this class. If
#' \code{info.manycov = NULL}, it will be computed during execution. Default is
#' \code{info.manycov = NULL}.
#' @param cov.ranges Matrix containing as its rows the lower and upper bounds
#' for each continuous covariate. Default is \code{cov.ranges = NULL}.
#' @param degree Degree of the spline functions to be used as instrumental
#' functions for the continuous covariates (if applicable). Default is
#' \code{degree = 3}.
#' 
G.cd.mc <- function(x, g.idx, data, n.if.per.cov, idxs.c, G.c, norm.func,
                    info.manycov = NULL, cov.ranges = NULL,
                    norm.cov.out = NULL, degree = 3, ...) {
  
  #### Precompute/preset some necessary variables ####
  
  # Obtain vector of covariate names
  cov.names <- colnames(data)[grep("X[1-9][[:digit:]]*$", colnames(data))]
  
  # If the necessary information is not pre-supplied...
  if (is.null(info.manycov)) {
    
    # Obtain each pair of covariates in the data. For each, determine the amount
    # of instrumental functions.
    info.manycov <- data.frame(cov.pair = character(), n.if = numeric())
    for (cov.name.idx1 in 1:length(cov.names)) {
      for (cov.name.idx2 in 2:length(cov.names)) {
        if (cov.name.idx2 > cov.name.idx1) {
          
          # Name of covariates in the pair
          cov.name1 <- cov.names[cov.name.idx1]
          cov.name2 <- cov.names[cov.name.idx2]
          
          # Number of instrumental functions for each
          n.if1 <- ifelse(cov.name.idx1 %in% idxs.c, n.if.per.cov, length(unique(data[, cov.name1])))
          n.if2 <- ifelse(cov.name.idx2 %in% idxs.c, n.if.per.cov, length(unique(data[, cov.name2])))
          
          # Total number of instrumental functions
          n.if <- n.if1 * n.if2
          
          # Add to information data frame
          row <- list(cov.pair = sprintf("%s, %s", cov.name1, cov.name2),
                      n.if = n.if)
          info.manycov <- rbind(info.manycov, row)
        }
      }
    }
    
    # Add supplementary rows and columns
    info.manycov <- cbind(idx = 1:nrow(info.manycov),
                          info.manycov,
                          cumsum = cumsum(info.manycov$n.if))
    info.manycov <- rbind(list(idx = 0, cov.pair = "", n.if = 0, cumsum = 0),
                          info.manycov)
  }
  
  #### Select the relevant variables ####
  
  # Get pair of covariates corresponding to the given index of instrumental
  # function.
  cov.pair <- info.manycov[min(which(info.manycov$cumsum >= g.idx)), "cov.pair"]
  pair.vec <- strsplit(cov.pair, split = ", ")[[1]]
  
  # Get subset of data and covariate vector corresponding to the covariates
  data.sub <- data[, c("Y", "Delta", "X0", pair.vec)]
  x.sub <- x[which(cov.names %in% pair.vec)]
  g.idx.sub <- g.idx - info.manycov[min(which(info.manycov$cumsum >= g.idx)) - 1, "cumsum"]
  cov.ranges.sub <- cov.ranges[, pair.vec]
  
  #### Construct the class of instrumental functions for this pair ####
  
  # If both variables in the pair are continuous...
  if (all(which(cov.names %in% pair.vec) %in% idxs.c)) {
    eval <- G.c(x = x.sub, g.idx = g.idx.sub, data = data.sub,
                n.if.per.cov = n.if.per.cov, norm.func = norm.func,
                cov.ranges = cov.ranges.sub, degree = degree)
    
  # If both variables in the pair are binary...
  } else if (all(!(which(cov.names %in% pair.vec) %in% idxs.c))) {
    levels.var1 <- unique(data[, pair.vec[1]])
    levels.var2 <- unique(data[, pair.vec[2]])
    discrete.covariate.levels <- expand.grid(levels.var1, levels.var2)
    discrete.covariate.levels[order(discrete.covariate.levels[,1],
                                    discrete.covariate.levels[,2],
                                    decreasing = TRUE),]
    eval <- as.numeric(all(discrete.covariate.levels[g.idx.sub, ] == x.sub))
    
  # If one variable in the pair in continuous and the other one is binary...
  } else {
    cont.var.idx <- which(pair.vec %in% cov.names[idxs.c])
    disc.var.idx <- setdiff(1:2, cont.var.idx)
    levels.disc <- sort(unique(data[, pair.vec[disc.var.idx]]))
    n.levels.disc <- length(levels.disc)
    
    g.idx.sub.d <- ((g.idx.sub - 1) %% n.levels.disc) + 1
    g.idx.sub.c <- ((g.idx.sub - 1) %/% n.levels.disc) + 1
    data.sub.c <- data[, c("Y", "Delta", "X0", pair.vec[cont.var.idx])]
    cov.ranges.sub.c <- cov.ranges.sub[, pair.vec[cont.var.idx]]
    
    eval.d <- as.numeric(x.sub[disc.var.idx] == levels.disc[g.idx.sub.d])
    eval.c <- G.c(x = x.sub[cont.var.idx], g.idx = g.idx.sub.c,
                  data = data.sub.c, n.if.per.cov = n.if.per.cov,
                  norm.func = norm.func, cov.ranges = cov.ranges.sub.c,
                  degree = degree)
    
    eval <- eval.d * eval.c
  }
  
  #### Return the result ####
  
  eval
}

#' @title Evaluate each instrumental function at each of the observations.
#' 
#' @description Obtain the evaluations of each observation on each of the
#' instrumental functions. (Used in function get.mi.mat.R)
#' 
#' @param data Data frame.
#' @param hp List of hyperparameters. Notably, it contains the instrumental
#' function to be used in an element named \code{G}.
#' 
get.instrumental.function.evals <- function(data, hp) {
  
  # Unpack hyperparameters
  n.inst.func <- hp[["n.inst.func"]]
  G <- hp[["G"]]
  
  # Initialize matrix that will store the instrumental function evaluations
  inst.func.evals <- matrix(nrow = nrow(data), ncol = n.inst.func)
  
  for (i in 1:nrow(data)) {
    
    # Get the covariate values of the i-th observation. Leave out the intercept.
    X <- as.matrix(data[i, grepl("X[[:digit:]]+", colnames(data))])
    X.no_int <- X[-1]
    
    # For each instrumental function, evaluate it at the covariates values of 
    # the i-th observation.
    for (j in 1:n.inst.func) {
      inst.func.evals[i, j] <- G(X.no_int, j)
    }
  }
  
  # Return the results
  inst.func.evals
}

#### Moment functions + derivatives ####

#' @title [DEPRECATED] Component function of the vector of moment functions m.
#' 
#' @description THIS FUNCTION IS DEPRECATED AND WILL THROW A WARNING WHEN USED.
#' Use the faster function 'get.mi.mat.R' instead
#' 
#' @param i Index of observation
#' @param j Index of moment function.
#' @param data Data frame.
#' @param beta Vector of coefficients.
#' @param t Time point at which to compute the moment function
#' @param hp List of hyperparamerers.
#' 
#' @noRd
#' 
m.comp <- function(i, j, data, beta, t, hp) {
  
  # This function is deprecated
  warning("Attempted to use a deprecated function (m.comp.R)")
  
  # Unpack data
  Y <- data[i, "Y"]
  Delta <- data[i, "Delta"]
  X <- as.matrix(data[i, grepl("X[[:digit:]]+", colnames(data))])
  X.no_int <- matrix(X[-1], nrow = 1)
  
  # Unpack hyperparameters
  Lambda <- hp[["Lambda"]]
  G <- hp[["G"]]
  n.inst.func <- hp[["n.inst.func"]]
  
  # Compute moment function
  if (j <= n.inst.func) {
    (Lambda(X %*% beta) - as.numeric(Y <= t & Delta == 1)) * G(X.no_int, j)
  } else {
    (as.numeric(Y <= t) - Lambda(X %*% beta) ) * G(X.no_int, j - n.inst.func)
  }
}

#' @title [DEPRECATED] Vector of moment functions
#' 
#' @description THIS FUNCTION IS DEPRECATED AND WILL THROW A WARNING WHEN USED.
#' Use the faster function 'get.mi.mat.R' instead
#' 
#' @param i Index of observation
#' @param data Data frame.
#' @param beta Vector of coefficients.
#' @param t Time point at which to compute the moment function.
#' 
#' @noRd
#' 
m <- function(i, data, beta, t, hp) {
  
  # This function is deprecated
  warning("Attempted to use a deprecated function (m.R).")
  
  # Number of instrumental functions
  n.inst.func <- hp[["n.inst.func"]]
  
  # Create vector of moment functions evaluated at (data, theta)
  rtrn <- c()
  for (j in 1:(2*n.inst.func)) {
    rtrn <- c(rtrn, m.comp(i, j, data, beta, t, hp))
  }
  
  # Return the results
  rtrn
}

#' @title Compute the conditional moment evaluations
#' 
#' @description This function computes the 1(Y <= t) - Lambda(X^T beta(t)) and
#' Lambda(X^T beta(t)) - 1(Y <= t, Delta = 1) parts of the moment functions.
#' (Used in function get.mi.mat.R)
#' 
#' @param data Data frame.
#' @param beta Vector of coefficients.
#' @param t Time point of interest.
#' @param hp List of hyperparameters.
#' 
#' @returns A vector of 2n elements containing in the first n positions the
#' evaluations of 1(Y <= t) - Lambda(X^T beta(t)) and in the last n positions
#' the evaluations of Lambda(X^T beta(t)) - 1(Y <= t, Delta = 1).
#' 
get.cond.moment.evals <- function(data, beta, t, hp) {
  
  # Unpack hyperparameters
  Lambda <- hp[["Lambda"]]
  
  # Initialize matrix that will store the results
  evals <- matrix(nrow = nrow(data), ncol = 2)
  
  # For each observation, compute the evaluation of the two conditional moment
  # functions
  for (i in 1:nrow(data)) {
    
    # Get the values pertaining to the i-th observation
    Y <- data[i, "Y"]
    Delta <- data[i, "Delta"]
    X <- as.matrix(data[i, grepl("X[[:digit:]]+", colnames(data))])
    
    # Compute moment functions
    evals[i, 1] <- Lambda(X %*% beta) - as.numeric(Y <= t & Delta == 1)
    evals[i, 2] <- as.numeric(Y <= t) - Lambda(X %*% beta)
  }
  
  # Return the results
  evals
}

#' @title Faster implementation of vector of moment functions.
#' 
#' @description
#' This function obtains the moment function evaluations.
#' 
#' 
#' @param i Index of observation
#' @param data Data frame.
#' @param beta Vector of coefficients.
#' @param t Time point at which to compute the moment function. Also allowed to
#' be a vector of time points (used in estimating the model under assumed time-
#' independent coefficients).
#' 
get.mi.mat <- function(data, beta, t, hp, inst.func.evals = NULL) {
  
  # Extract hyperparameters
  n.inst.func <- hp[["n.inst.func"]]
  n.cov <- sum(grepl("[1-9][[:digit:]]*$", colnames(data)))
  
  # Get instrumental function evaluations
  if (is.null(inst.func.evals)) {
    inst.func.evals <- t(get.instrumental.function.evals(data, hp))
  }
  
  # Create matrix of replicates of instrumental function evaluation. I.e.
  # ife = [inst.func.evals
  #        inst.func.evals
  #        ...
  #        inst.func.evals].
  ife <- do.call(rbind, replicate(2*length(t), inst.func.evals, simplify=FALSE))
  
  # Get conditional moment evaluations at each time point
  cmfe <- NULL
  for (time.point in t) {
    if (class(beta) == "function") {
      beta.t <- beta(time.point)
    } else if (length(t) == 1) {
      beta.t <- beta
    } else {
      beta.t <- beta[c(which(t == time.point), (length(t) + 1):length(beta))]
    }
    cond.m.evals <- get.cond.moment.evals(data, beta.t, time.point, hp)
    cmfe1 <- matrix(rep(cond.m.evals[,1], n.inst.func), nrow = n.inst.func, byrow = TRUE)
    cmfe2 <- matrix(rep(cond.m.evals[,2], n.inst.func), nrow = n.inst.func, byrow = TRUE)
    cmfe <- rbind(cmfe, cmfe1, cmfe2)
  }
  
  # Combine the conditional moment function with the instrumental function
  # evaluations and return the result.
  cmfe * ife
}

#' @title Vector of sample average of each moment function (\bar{m}_n(\theta)).
#' 
#' @description This function obtains the vector of sample averages of each
#' moment function.
#' 
#' 
#' @param data Data frame.
#' @param beta Vector of coefficients.
#' @param t Time point at which to compute the moment functions. Also allowed to
#' be a vector of time points (used in estimating the model under assumed time-
#' independent coefficients).
#' @param hp List of hyperparameters.
#' @param mi.mat Matrix of moment function evaluations. Can be used to avoid
#' some computation. Default is \code{mi.mat = NULL}.
#' 
m.bar <- function(data, beta, t, hp, mi.mat = NULL) {
  
  # Number of instrumental functions
  n.inst.func <- hp[["n.inst.func"]]
  
  # Initialize vector that will contain the sum of all moment function
  # evaluations.
  m.evals.sum <- rep(0, 2*length(t)*n.inst.func)
  
  # Obtain the sum
  for (i in 1:nrow(data)) {
    if (is.null(mi.mat)) {
      mi <- NULL
      for (time.point in t) {
        mi <- c(mi, m(i, data, beta, time.point, hp))
      }
    } else {
      mi <- mi.mat[,i]
    }
    m.evals.sum <- m.evals.sum + mi
  }
  
  # Return the average
  m.evals.sum/nrow(data)
}

#' @title [DEPRECATED] Component function of the vector of derivatives of moment
#' functions (with respect to \beta).
#' 
#' @description This function obtains the vector of partial derivatives of a
#' moment function, evaluated at a specified observation. This function is
#' deprecated.
#' 
#' 
#' @param i Index of observation
#' @param j Index of moment function.
#' @param data Data frame.
#' @param beta Vector of coefficients.
#' @param t Time point at which to compute the derivative of the moment function
#' (not actually used in the implementation below)
#' @param hp List of hyperparameters
#' 
#' @returns A vector containing the partial derivatives of the selected moment
#' function, evaluated at the specified observation.
#' 
#' @noRd
#' 
dm.comp <- function(i, j, data, beta, t, hp) {
  
  # Unpack data
  Y <- data[i, "Y"]
  Delta <- data[i, "Delta"]
  X <- as.matrix(data[i, grepl("X[[:digit:]]+", colnames(data))])
  X.no_int <- matrix(X[-1], nrow = 1)
  
  # Unpack hyperparameters
  dLambda <- hp[["dLambda"]]
  G <- hp[["G"]]
  n.inst.func <- hp[["n.inst.func"]]
  
  # Compute vector derivative of moment function
  if (j <= n.inst.func) {
    G(X.no_int, j) * dLambda(as.numeric(X %*% beta)) * X
  } else {
    - G(X.no_int, j - n.inst.func) * dLambda(as.numeric(X %*% beta)) * X
  }
}

#' @title [DEPRECATED] Vector of derivatives of moment functions
#' 
#' @description This function returns a matrix containing the partial
#' derivatives of each moment function, evaluated at the specified observation.
#' 
#' @param i Index of observation
#' @param data Data frame.
#' @param beta Vector of coefficients.
#' @param t Time point at which to compute the derivative of the moment function
#' @param hp List of hyperparameters.
#' 
#' @returns A matrix containing the partial derivatives of each moment
#' function, evaluated at the specified observation. Each row corresponds to a
#' moment function, each column corresponds to a coefficient.
#' 
#' @noRd
#' 
dm <- function(i, data, beta, t, hp) {
  
  # Warn user that the function is deprecated and might not work
  warning("Using deprecated function 'dm'!")
  
  # Number of instrumental functions
  n.inst.func <- hp[["n.inst.func"]]
  
  # Create vector of moment functions evaluated at (data, theta)
  rtrn <- NULL
  for (j in 1:(2*n.inst.func)) {
    rtrn <- rbind(rtrn, dm.comp(i, j, data, beta, t, hp))
  }
  
  # Return the results
  rtrn
}

#' @title Matrix of derivatives of conditional moment functions
#' 
#' @description This function evaluates the derivatives of the conditional
#' moment function at each observation. Used in get.dmi.tens.R
#' 
#' @param data Data frame.
#' @param beta Parameter vector.
#' @param t Time point of interest.
#' @param hp List of hyperparameters.
#' 
get.deriv.mom.func <- function(data, beta, t, hp) {
  
  # Extract hyperparameters
  n.inst.func <- hp[["n.inst.func"]]
  dLambda <- hp[["dLambda"]]
  n.param <- length(beta)
  n <- nrow(data)
  
  # Initialize matrix that will store the results
  evals.m1 <- matrix(nrow = n, ncol = n.param)
  evals.m2 <- matrix(nrow = n, ncol = n.param)
  
  # For each observation, compute the evaluation of the two conditional moment
  # functions
  for (i in 1:nrow(data)) {
    
    # Get the values pertaining to the i-th observation
    Y <- data[i, "Y"]
    Delta <- data[i, "Delta"]
    X <- as.matrix(data[i, grepl("X[[:digit:]]+", colnames(data))])
    
    # Compute derivatives of moment functions
    evals.m1[i,] <- dLambda(as.numeric(X %*% beta)) * X
  }
  evals.m2 <- -evals.m1
  
  # Return the results
  evals <- array(dim = c(n, n.param, 2))
  evals[, , 1] <- evals.m1
  evals[, , 2] <- evals.m2
  evals
}

#' @title Faster implementation to obtain the tensor of the evaluations of the
#' derivatives of the moment functions at each observation.
#' 
#' @description This function provides a faster implementation of obtaining the
#' evaluations of the derivative of the moment functions at each observation
#' (wrt the previous implementation using 'dm.comp' and 'dm.R'). Used in the
#' function G.hat.R
#' 
#' @param data Data frame.
#' @param beta Vector of coefficients.
#' @param t Time point of interest. Also allowed to
#' be a vector of time points (used in estimating the model under assumed time-
#' independent coefficients).
#' @param hp List of hyperparameters.
#' @param inst.func.evals Precomputed matrix of instrumental function
#' evaluations. Defaults is \code{inst.func.evals = NULL}, in which case the
#' evaluations will be done inside this function.
#' 
get.dmi.tens <- function(data, beta, t, hp, inst.func.evals = NULL) {
  
  # Extract hyperparameters
  n.inst.func <- hp[["n.inst.func"]]
  n.param <- length(beta)
  n <- nrow(data)
  
  # Compute the instrumental function evaluations if necessary
  if (is.null(inst.func.evals)) {
    inst.func.evals <- t(get.instrumental.function.evals(data, hp))
  }
  
  # Repeat the matrix of instrumental function evaluations into a tensor of the
  # correct dimension (note that the instrumental function evaluations do not
  # depend on the covariates).
  inst.func.tens <- array(dim = c(2*n.inst.func*length(t), n.param, n))
  for (i in 1:n.param) {
    inst.func.tens.i <- inst.func.evals
    for (dummy in 2:(2*length(t))) {
      inst.func.tens.i <- rbind(inst.func.tens.i, inst.func.evals)
    }
    inst.func.tens[, i, ] <- inst.func.tens.i
  }
  
  # Compute matrix of derivatives of conditional moment functions, evaluated at
  # each observation.
  deriv.mom.evals.list <- list()
  for (time.point in t) {
    if (class(beta) == "function") {
      beta.t <- beta(time.point)
    } else if (length(t) == 1) {
      beta.t <- beta 
    } else {
      beta.t <- beta[c(which(t == time.point), (length(t)+1):length(beta))]
    }
    
    # Compute the derivatives of the moment functions at each time point wrt the
    # appropriate intercept parameter. Derivatives wrt intercept parameters at
    # other time points are zero
    deriv.mom.func <- array(0, dim = c(n, n.param, 2))
    deriv.mom.func[, c(which(t == time.point), (length(t)+1):n.param), ] <-
      get.deriv.mom.func(data, beta.t, time.point, hp)
    
    # Store the result
    deriv.mom.evals.list[[as.character(time.point)]] <- deriv.mom.func
  }
  
  # Create tensor of evaluations of the derivatives of the moment functions
  deriv.cond.mom.func <- array(dim = c(2*length(t)*n.inst.func, n.param, n))
  for (time.point.idx in 1:length(deriv.mom.evals.list)) {
    for (j in 1:n.inst.func) {
      deriv.cond.mom.func[2*n.inst.func*(time.point.idx - 1) + j, ,] <-
        t(deriv.mom.evals.list[[time.point.idx]][, , 1])
      deriv.cond.mom.func[2*n.inst.func*(time.point.idx - 1) + j + n.inst.func, ,] <-
        t(deriv.mom.evals.list[[time.point.idx]][, , 2])
    }
  }
  
  # Compute tensor of evaluations of derivatives of unconditional moment
  # functions.
  dmi.tens <- inst.func.tens * deriv.cond.mom.func
}

#' @title Vector of sample average of each moment function (\bar{m}_n(\theta)).
#' 
#' @description This function computes the matrix containing the sample average
#' of the partial derivatives of the moment functions.
#' 
#' @param data Data frame.
#' @param beta Vector of coefficients.
#' @param t Time point at which to compute the derivative of the moment
#' functions. Also allowed to
#' be a vector of time points (used in estimating the model under assumed time-
#' independent coefficients).
#' @param hp List of hyperparameters.
#' @param dmi.mat Matrix of derivative moment function evaluations. Can be used
#' to avoid some computation. Default is \code{mi.mat = NULL}.
#'
#' @returns A matrix containing the sample average of the partial derivatives of
#' the moment functions. Each row corresponds to a moment function, each column
#' corresponds to a coefficient.
dm.bar <- function(data, beta, t, hp, dmi.tens = NULL) {
  
  # Number of instrumental functions
  n.inst.func <- hp[["n.inst.func"]]
  
  # Number of covariates
  n.param <- length(beta)
  
  # Initialize vector that will contain the sum of all moment function
  # evaluations.
  dm.evals.sum <- matrix(0, nrow = 2*length(t)*n.inst.func, ncol = n.param)
  
  # Obtain the sum
  for (i in 1:nrow(data)) {
    if (is.null(dmi.tens)) {
      dmi <- dm(i, data, beta, t, hp)
    } else {
      dmi <- dmi.tens[,,i]
    }
    dm.evals.sum <- dm.evals.sum + dmi
  }
  
  # Return the average
  dm.evals.sum/nrow(data)
}

#### Variances-covariance/correlation of moment functions ####

#' @title Compute the variance-covariance matrix of the moment functions.
#' 
#' @description This function comptutes the empricical variance-covariance
#' matrix of the moment functions.
#' 
#' @param data Data frame.
#' @param beta Coefficient vector.
#' @param t Time point of interest.
#' @param hp List of hyperparameters.
#' @param m.avg A precomputed vector of the sample average of the moment
#' functions. If not supplied, this vector is computed. Default is
#' \code{m.avg = NULL}.
#' @param mi.mat A precomputed matrix of moment function evaluations at each
#' observation. If supplied, some computations can be skipped. Default is
#' \code{mi.mat = NULL}.
#' 
Sigma.hat <- function(data, beta, t, hp, m.avg = NULL, mi.mat = NULL) {
  
  # Number of instrumental functions
  n.inst.func <- hp[["n.inst.func"]]
  
  # Sample average of the moment functions
  if (is.null(m.avg)) {
    m.avg <- m.bar(data, beta, t, hp)
  }
  
  # Initialize matrix that will contain the sum of all outer products used in
  # obtaining the sample variance-covariance matrix
  sig.evals.sum <- matrix(0, nrow = 2*length(t)*n.inst.func, ncol = 2*length(t)*n.inst.func)
  
  # Obtain the sum
  for (i in 1:nrow(data)) {
    if (is.null(mi.mat)) {
      mi <- m(i, data, beta, t, hp)
    } else {
      mi <- mi.mat[,i]
    }
    sig.evals.sum <- sig.evals.sum + outer(mi - m.avg, mi - m.avg)
  }
  
  # Return the average
  sig.evals.sum/nrow(data)
}

#' @title Obtain the diagonal matrix of sample variances of moment functions
#' 
#' @description This function computes the diagonal matrix of the sample
#' variance-covariance matrix. 
#' 
#' @param input Can either be the variance-covariance matrix obtained from the
#' function Sigma.hat, or the data frame.
#' @param beta The coefficient vector. Only needs to be supplied when the
#' argument for \code{input} is the data frame.
#' @param t The time point of interest. Only needs to be supplied when the
#' argument for \code{input} is the data frame.
#' @param hp List of hyperparameters. Only needs to be supplied when the
#' argument for \code{input} is the data frame.
#' @param m.avg See documentation of \code{Sigma.hat}. Only needs to be supplied
#' when the argument for \code{input} is the data frama.
#' @param mi.mat See documentation of \code{Sigma.hat}. Only needs to be supplied
#' when the argument for \code{input} is the data frama.
#' 
D.hat <- function(input, beta = NULL, t = NULL, hp = NULL, m.avg = NULL,
                  mi.mat = NULL) {
  
  # If the given input is a variance-covariance matrix...
  if (is.matrix(input)) {
    Sigma <- input
    return(diag(diag(Sigma), nrow = nrow(Sigma)))
    
  # If the given input is a data frame...
  } else {
    
    data <- input
    
    # beta, t and hp should be specified in this case.
    if (any(is.null(c(beta, t, hp)))) {
      stop("When input is not a matrix, beta, t and hp should be specified.")
    }
    
    # Evaluations of the moment functions at each observation
    if (is.null(mi.mat)) {
      for (i in 1:nrow(data)) {
        mi.mat <- cbind(mi.mat, m(i, data, beta, t, hp))
      }
    }
    
    # Sample average of the moment functions
    if (is.null(m.avg)) {
      m.avg <- m.bar(data, beta, t, hp, mi.mat = mi.mat)
    }
    
    # Initialize variable that will store the result
    sum.sq.dev <- rep(0, length(m.avg))
    
    # Obtain the sum of squared deviations from the sample mean
    for (i in 1:nrow(data)) {
      mi <- mi.mat[,i]
      sum.sq.dev <- sum.sq.dev + (mi - m.avg)^2
    }
    
    # Return the resuls
    return(diag(sum.sq.dev/nrow(data), nrow = length(m.avg)))
  }
}

#' @title Obtain the correlation matrix of the moment functions
#' 
#' @description This function computes the correlation matrix corresponding to
#' the variance-covariance matrix as returned by \code{Sigma.hat.R}
#' 
#' @param Sigma The output of the function Sigma.hat
#' 
Omega.hat <- function(Sigma) {
  sqrt.D.hat.inv <- solve(sqrt(D.hat(Sigma)))
  sqrt.D.hat.inv %*% Sigma %*% sqrt.D.hat.inv
}

#' @title Obtain the matrix of partial derivatives of the sample variances.
#' 
#' @description This function computes the matrix of sample derivatives of the 
#' sample variances.
#' 
#' @param data Data frame.
#' @param beta Vector of coefficients.
#' @param t Time point at which to evaluate the (derivatives of) the moment
#' functions. Also allowed to
#' be a vector of time points (used in estimating the model under assumed time-
#' independent coefficients).
#' @param hp List of hyperparamerers.
#' @param mi.mat A precomputed matrix of moment function evaluations at each
#' observation. If supplied, some computations can be skipped. Default is
#' \code{mi.mat = NULL}.
#' @param m.avg A precomputed vector of the sample average of the moment
#' functions. If not supplied, this vector is computed. Default is
#' \code{m.avg = NULL}.
#' @param dm.avg Matrix of precomputed sample averages of the derivatives of the
#' moment functions. Default is \code{dm.avg = NULL}.
#' @param dmi.tens 3D tensor of precomputed evaluations of the derivatives of 
#' the moment functions. Default is \code{dmi.tens = NULL}.
#' 
#' @returns A matrix containing the partial derivatives of the variances of the
#' moment functions. Each row corresponds to a moment function, each column
#' corresponds to a covariate.
dD.hat <- function(data, beta, t, hp, mi.mat = NULL, m.avg = NULL,
                   dm.avg = NULL, dmi.tens = NULL) {
  
  # Define some useful variables
  n <- nrow(data)
  n.param <- length(beta)
  n.inst.func <- hp[["n.inst.func"]]
  
  # Evaluations of the moment functions at each observation
  if (is.null(mi.mat)) {
    for (i in 1:nrow(data)) {
      mi.mat <- cbind(mi.mat, m(i, data, beta, t, hp))
    }
  }
  
  # Sample average of the moment functions
  if (is.null(m.avg)) {
    m.avg <- m.bar(data, beta, t, hp, mi.mat = mi.mat)
  }
  
  # Evaluations of the derivatives of moment functions at each observation
  if (is.null(dmi.tens)) {
    dmi.tens <- array(dim = c(2*n.inst.func, n.param, n))
    for (i in 1:n) {
      dmi.tens[, ,i] <- dm(i, data, beta, t, hp)
    }
  }
  
  # Sample average of the derivatives of moment functions
  if (is.null(dm.avg)) {
    dm.avg <- dm.bar(data, beta, t, hp, dmi.tens = mi.tens)
  }
  
  # Compute dD.hat
  sum <- 0
  for (i in 1:n) {
    mi <- mi.mat[,i]
    dmi <- dmi.tens[,,i]
    sum <- sum + 2 * matrix(rep(mi - m.avg, n.param), ncol = n.param) * (dmi - dm.avg)
  }
  
  # Return the result
  sum/n
}

#' @title Compute the Gn matrix in step 3b of Bei (2024).
#' 
#' @param data Data frame.
#' @param beta Vector of coefficients.
#' @param t Time point at which to evaluate the (derivatives of) the moment
#' functions.
#' @param hp List of hyperparamerers.
#' @param mi.mat A precomputed matrix of moment function evaluations at each
#' observation. If supplied, some computations can be skipped. Default is
#' \code{mi.mat = NULL}.
#' @param m.avg A precomputed vector of the sample average of the moment
#' functions. If not supplied, this vector is computed. Default is
#' \code{m.avg = NULL}.
#' @param dm.avg Matrix of precomputed sample averages of the derivatives of the
#' moment functions. Default is \code{dm.avg = NULL}.
#' @param dmi.tens 3D tensor of precomputed evaluations of the derivatives of 
#' the moment functions. Default is \code{dmi.tens = NULL}.
#' @param D Diagonal of D-matrix.
#' 
#' @returns A matrix containing the partial derivatives of the variances of the
#' moment functions. Each row corresponds to a moment function, each column
#' corresponds to a covariate.
#' 
#' @references Bei, X. (2024). Local linearieation based subvector inference in
#' moment inequality models. Journal of Econometrics. 238:105549-
#' 
G.hat <- function(data, beta, t, hp, mi.mat = NULL, m.avg = NULL,
                  dm.avg = NULL, dmi.tens = NULL, D = NULL) {
  
  # Define some useful variables
  n <- nrow(data)
  n.param <- length(beta)
  n.inst.func <- hp[["n.inst.func"]]
  
  # Evaluations of the moment functions at each observation
  if (is.null(mi.mat)) {
    mi.mat <- get.mi.mat(data, beta, t, hp)
  }
  
  # Sample average of the moment functions
  if (is.null(m.avg)) {
    m.avg <- m.bar(data, beta, t, hp, mi.mat = mi.mat)
  }
  
  # Evaluations of the derivatives of moment functions at each observation
  if (is.null(dmi.tens)) {
    dmi.tens <- get.dmi.tens(data, beta, t, hp)
  }
  
  # Sample average of the derivatives of moment functions
  if (is.null(dm.avg)) {
    dm.avg <- dm.bar(data, beta, t, hp, dmi.tens = dmi.tens)
  }
  
  # Compute the diagonal of Dn
  if (is.null(D)) {
    Dn <- diag(D.hat(data, beta, t, hp, m.avg = m.avg, mi.mat = mi.mat))
  } else {
    Dn <- diag(D)
  }
  
  # Compute the derivative of Dn
  dDn <- dD.hat(data, beta, t, hp, mi.mat = mi.mat, m.avg = m.avg,
                dm.avg = dm.avg, dmi.tens = dmi.tens)
  
  # Compute Gn
  a <- dm.avg * matrix(rep(sqrt(Dn), n.param), ncol = n.param) 
  b <- matrix(rep(m.avg, n.param), ncol = n.param) * (1/2) * matrix(rep(Dn^(-1/2), n.param), ncol = n.param) * dDn
  
  (a - b)/matrix(rep(Dn, n.param), ncol = n.param)
}

#### S functions ####

#' @title S-function
#' 
#' @description This function computes the loss function at a given point.
#' 
#' @param m Vector of averages of moment functions.
#' @param Sigma Sample variance-covariance matrix of moment functions.
#' 
#' @returns S(m, Sigma).
#' 
S.func <- function(m, Sigma) {
  
  # Number of moment functions
  p <- length(m)
  
  # Initialize variable
  S <- 0
  for (j in 1:p) {
    S <- S + max(-m[j]/sqrt(Sigma[j, j]), 0)^2
  }
  
  # Return results
  S
}

#### Estimating the test statistic ####

#' @title 'Loss function' of the test statistic.
#' 
#' @description This function implements the loss function used in computing
#' the test statistic.
#' 
#' @param beta.sub Subvector of coefficient vector.
#' @param data Data frame.
#' @param t Time point of interest. Also allowed to
#' be a vector of time points (used in estimating the model under assumed time-
#' independent coefficients).
#' @param hp List of hyperparameters.
#' @param c Unit vector containing unity at the location of the parameter of
#' interest.
#' @param r Value of the parameter of interest that is tested.
#' @param inst.func.evals Pre-computed matrix of insturmental function
#' evaluations. If not supplied, it will be computed during execution of this
#' function.
#' 
#' @returns S-functions evaluation for the specified parameter vector.
#' 
lf.ts <- function(beta.sub, data, t, hp, c, r, inst.func.evals = NULL) {
  
  # Sample size
  n <- nrow(data)
  
  # Minimum variance (used for computational reasons)
  min.var <- hp[["min.var"]]
  
  # Make the completed parameter vector
  if (length(t) == 1) {
    beta <- rep(r, length(beta.sub) + 1)
    beta[which(c == 0)] <- beta.sub
  } else {
    beta <- function(time.point) {
      beta <- rep(r, length(beta.sub) + 1)
      beta[which(c == 0)] <- beta.sub
      beta[1:length(t)] <- cumsum(beta[1:length(t)])
      beta <- beta[-which(t != time.point)]
      beta
    }
  }
  
  # Matrix of moment function evaluations
  mi.mat <- get.mi.mat(data, beta, t, hp, inst.func.evals)
  
  # Sample average of the moment functions
  m.avg <- m.bar(data, beta, t, hp, mi.mat = mi.mat)
  
  # Sample variance-covariance matrix
  svc <- Sigma.hat(data, beta, t, hp, mi.mat = mi.mat, m.avg = m.avg)
  
  # Ensure the invertibility of the sample variance-covariance matrix
  svc <- svc + min.var * diag(ncol(svc))
  
  # Sample variance diagonal matrix
  D <- D.hat(svc)
  
  # S-function
  S.func(sqrt(n) * diag(D)^(-1/2) * m.avg, Omega.hat(svc))
}

#' @title Obtain the test statistic by minimizing the S-function over the
#' feasible region \Beta(r)
#' 
#' @param beta.init Starting value of minimization algorithm.
#' @param data Data frame.
#' @param par.space Matrix containing the bounds on the parameter space.
#' @param t Time point at which to evaluate beta. Also allowed to
#' be a vector of time points (used in estimating the model under assumed time-
#' independent coefficients).
#' @param hp List of hyperparameters.
#' @param c Projection vector
#' @param r hypothesised value of the projection.
#' @param inst.func.evals Matrix of precomputed instrumental function
#' evaluations for each observation in the data set. If \code{NULL}, the
#' evaluations will be computed during execution of this function. Default is
#' \code{inst.func.evals = NULL}.
#' 
#' @returns A list containing the value of the test statistic and the parameter
#' at which this value was attained.
#' 
#' @import stats nloptr
#' 
get.test.statistic <- function(beta.init, data, par.space, t, hp, c, r,
                               inst.func.evals = NULL) {
  
  # Define some useful parameters
  n.param <- length(c)
  
  # If data.init represents the full vector, check whether it satisfies the
  # constraint and transform it into the unconstrained subvector
  if (length(beta.init) == n.param) {
    if (beta.init[which(c == 1)] != r) {
      stop("Given beta vector does not satisfy the constraint")
    }
    beta.init <- beta.init[which(c == 0)]
  }
  
  # Precompute the instrumental function evaluations
  if (is.null(inst.func.evals)) {
    inst.func.evals <- t(get.instrumental.function.evals(data, hp))
  }
  
  # Estimate the test statistic and the minimizer
  use.optim <- FALSE # Use stats::optim for optimization. Testing shows that
                     # this may miss the global optimum, leading to over-
                     # rejection.
  if (use.optim) {
    out <- optim(beta.init, lf.ts, data = data, t = t, hp = hp, c = c, r = r,
                 inst.func.evals = inst.func.evals,
                 method = "L-BFGS-B", lower = par.space[which(c == 0), 1],
                 upper = par.space[which(c == 0), 2], control = list(maxit = 200))
    
    Jnrh <- out$value
    beta.hat <- rep(r, length(c))
    beta.hat[which(c == 0)] <- out$par
  } else {
    out <- nloptr(beta.init,
                  eval_f = lf.ts,
                  lb = par.space[which(c == 0), 1],
                  ub = par.space[which(c == 0), 2],
                  opts = list("algorithm" = "NLOPT_LN_NEWUOA_BOUND",
                              "xtol_rel" = 1e-4,
                              "maxeval" = 1000),
                  data = data, t = t, hp = hp, c = c, r = r,
                  inst.func.evals = inst.func.evals)
    Jnrh <- out$objective
    beta.hat <- rep(r, length(c))
    beta.hat[which(c == 0)] <- out$solution
  }
  
  # Return the results
  return(list(Jnrh, beta.hat))
}

#### Calculating the critical value of the test statistic ####

#' @title Loss function to compute Delta(beta).
#' 
#' @description This function defines the loss function used in computing the
#' penalized local linear approximation of the test statistic in order to
#' construct the bootstrap distribution of the test statistic.
#' 
#' @param Delta.sub Subvector of Delta.
#' @param vnb Bootstrapped stochastic process.
#' @param phi Moment selection functions.
#' @param Gn First-order approximation matrix.
#' @param Omegan Correlation matrix of sample moment functions.
#' @param beta Coefficient vector.
#' @param c Projection vector.
#' @param r Value of projected coefficient vector.
#' @param data Data frame.
#' @param par.space Matrix containing the bounds on the parameter space.
#' @param epsilon.n Parameter used in constructing the feasible region as in
#' Example 4.1 in Bei (2024). Not used in this function.
#' @param lambda.n Weight of penalty term.
#' 
#' @returns Loss function evaluation evaluated at the given Delta.
#' 
#' @references Bei, X. (2024). Local linearieation based subvector inference in
#' moment inequality models. Journal of Econometrics. 238:105549-
#' 
lf.delta.beta1 <- function(Delta.sub, vnb, phi, Gn, Omegan, beta, c, r, data,
                           par.space, epsilon.n, lambda.n) {
  
  # Extract some useful parameters
  n <- nrow(data)
  
  # Make the completed parameter vector
  Delta <- rep(0, length(Delta.sub) + 1)
  Delta[which(c == 0)] <- Delta.sub
  
  # Value of the loss function
  S.func(vnb + phi + Gn %*% Delta, Omegan) + lambda.n/n * sum(Delta^2)
  
}

#' @title Compute the critical value of the test statistic.
#' 
#' @description This function computes the critical value following the 
#' algorithm of Section 4.3 in Bei (2024).
#' 
#' @param BetaI.r Matrix containing in its columns the minimizers of the
#' S-function leading to the test statistic.
#' @param data Data frame.
#' @param t Time point of interest. Also allowed to
#' be a vector of time points (used in estimating the model under assumed time-
#' independent coefficients).
#' @param hp List of hyperparameters.
#' @param c Projection vector.
#' @param r Result of projection of parameter vector onto \code{c}.
#' @param par.space Bounds on the parameter space.
#' @param inst.func.evals Matrix of precomputed instrumental function
#' evaluations for each observation in the data set. If \code{NULL}, the
#' evaluations will be computed during execution of this function. Default is
#' \code{inst.func.evals = NULL}.
#' @param alpha Confidence level.
#' 
#' @returns The critical value for the test statistic.
#' 
#' @import stats
#' 
#' @references Bei, X. (2024). Local linearieation based subvector inference in
#' moment inequality models. Journal of Econometrics. 238:105549-
#' 
get.cvLLn <- function(BetaI.r, data, t, hp, c, r, par.space,
                      inst.func.evals = NULL, alpha = 0.95) {
  
  # Define variables that will be useful throughout
  n <- nrow(data)
  J <- hp[["n.inst.func"]]*2
  n.beta <- ncol(BetaI.r)
  n.param <- nrow(BetaI.r)
  B <- hp[["B"]]
  kappa.n <- hp[["kappa.n"]]
  epsilon.n <- hp[["epsilon.n"]]
  lambda.n <- hp[["lambda.n"]]
  min.var <- hp[["min.var"]]
  
  # Precompute instrumental function evaluations
  if (is.null(inst.func.evals)) {
    inst.func.evals <- t(get.instrumental.function.evals(data, hp))
  }
  
  # Precompute moment function evaluations for all parameters in BetaI.r
  mi.tens <- array(dim = c(J*length(t), n, n.beta))
  for (beta.idx in 1:n.beta) {
    beta <- BetaI.r[, beta.idx]
    mi.tens[, , beta.idx] <- get.mi.mat(data, beta, t, hp, inst.func.evals)
  }
  
  # Precompute sample averages of moment functions for all beta in BetaI.r
  m.avg.mat <- matrix(nrow = J*length(t), ncol = n.beta)
  for (beta.idx in 1:n.beta) {
    beta <- BetaI.r[, beta.idx]
    m.avg.mat[, beta.idx] <- m.bar(data, beta, t, hp,
                                   mi.mat = mi.tens[, , beta.idx])
  }
  
  # Precompute variance-covariance matrix of moment functions for all beta.
  # Ensure the invertibility of each.
  Sigma.tens <- array(dim = c(J*length(t), J*length(t), n.beta))
  for (beta.idx in 1:n.beta) {
    beta <- BetaI.r[, beta.idx]
    Sigma.tens[, , beta.idx] <- Sigma.hat(data, beta, t, hp,
                                          m.avg = m.avg.mat[, beta.idx],
                                          mi.mat = mi.tens[ , , beta.idx])
    
    Sigma.tens[, , beta.idx] <- Sigma.tens[, , beta.idx] + min.var * diag(ncol(Sigma.tens[, , beta.idx]))
  }
  
  # Precompute the variance diagonal matrices for all beta in BetaI.r (stored
  # as vectors)
  D.diag.mat <- matrix(nrow = J*length(t), ncol = n.beta)
  for (beta.idx in 1:n.beta) {
    D.diag.mat[, beta.idx] <- diag(Sigma.tens[, , beta.idx])
  }
  
  # Precompute the square root of the inverse of the diagonal variance matrices
  D.inv.sqrt.diag.mat <- D.diag.mat^(-1/2)
  
  # Precompute phi(xi(beta)) of each beta in BetaI.r
  phi.mat <- matrix(nrow = J*length(t), ncol = n.beta)
  for (beta.idx in 1:n.beta) {
    beta <- BetaI.r[, beta.idx]
    phi <- rep(0, J)
    for (j in 1:J) {
      phi[j] <- max(sqrt(n) * m.avg.mat[j, beta.idx] * D.inv.sqrt.diag.mat[j, beta.idx] / kappa.n, 0)
    }
    phi.mat[, beta.idx] <- phi
  }
  
  # Precompute Gn(theta)
  Gn.tens <- array(dim = c(J*length(t), n.param, n.beta))
  for (beta.idx in 1:n.beta) {
    beta <- BetaI.r[, beta.idx]
    Gn.tens[, , beta.idx] <- G.hat(data, beta, t, hp,
                                   mi.mat = mi.tens[ , , beta.idx],
                                   m.avg = m.avg.mat[ , beta.idx],
                                   D = diag(D.diag.mat[, beta.idx], nrow  = J*length(t)))
  }
  
  # Initialize object that will store all bootstrapped test statistics
  JnLLb.r_vct <- NULL
  
  # For each bootstrap iteration, compute the bootstrapped test statistic
  for (b in 1:B) {
    
    # Simulate i.i.d standard normal random variables.
    zeta <- rnorm(n)
    
    # Initialize matrix to store all S-function evaluations needed to determine
    # the bootstrapped test statistic.
    S.evals <- NULL
    
    # Loop over all beta vectors inside BetaI.r and compute the corresponding
    # S-function evaluation.
    for (col.idx in 1:n.beta) {
      
      # Select the beta corresponding to this iteration
      beta <- BetaI.r[, col.idx]
      
      # Matrix of moment function evaluations
      mi.mat <- mi.tens[, , col.idx]
      
      # Sample average of the moment functions
      m.avg <- m.avg.mat[, col.idx]
      
      # Sample variance-covariance matrix
      Sigman <- Sigma.tens[, , col.idx]
      
      # Sample correlation matrix
      Omegan <- Omega.hat(Sigman)
      
      # (Inverse square root) sample variance matrix
      D.inv.sqrt <- diag(D.inv.sqrt.diag.mat[, col.idx], J*length(t))
      
      # Compute vnb
      vnb <- rep(0, J*length(t))
      for (i in 1:n) {
        vnb <- vnb + sqrt(1/n) * D.inv.sqrt %*% (mi.mat[,i] - m.avg) * zeta[i]
      }
      
      # phi(xi(theta))
      phi <- phi.mat[, col.idx]
      
      # Compute \hat{G}_n
      Gn <- Gn.tens[, , col.idx]
      
      # Evaluate the loss function at the origin. Store the result in a matrix
      delta.search <- matrix(c(rep(0, n.param), S.func(vnb + phi, Omegan)), ncol = n.param + 1)
      colnames(delta.search) <- c(paste0("X", 0:(n.param - 1)), "val")
      
      # Find the minimum of the loss function, starting from each corner point of
      # the feasible region.
      for (comb.nbr in 1:2^(n.param - 1)) {
        
        # Define the parameter bounds for the optimization
        optim.lb <- sqrt(n) * (par.space[which(c == 0), 1] + epsilon.n - beta[which(c == 0)])
        optim.ub <- sqrt(n) * (par.space[which(c == 0), 2] - epsilon.n - beta[which(c == 0)])
        
        # Get the corner point corresponding to this iteration
        comb <- base(comb.nbr - 1, 2, num.digits = n.param - 1)
        corner.point <- (1 - comb) * optim.lb + comb * optim.ub
        
        # Perform the optimization
        out <- optim(corner.point, lf.delta.beta1, vnb = vnb, phi = phi, Gn = Gn,
                     Omegan = Omegan, beta = beta, c = c, r = r, data = data,
                     par.space = par.space, epsilon.n = epsilon.n,
                     lambda.n = lambda.n, method = "L-BFGS-B",
                     lower = optim.lb, upper = optim.ub)
        
        # Extract the results
        val <- out$value
        solution <- rep(0, length(c))
        solution[which(c == 0)] <- out$par
        delta.search <- rbind(delta.search, c(solution, val))
      }
      
      # Obtain the 'global' minimum
      Delta.beta <- delta.search[which.min(delta.search[, "val"]), 1:n.param]
      
      # S function evaluation
      S.evals <- c(S.evals, S.func(vnb + phi + Gn %*% Delta.beta, Omegan))
    }
    
    # Compute the bootstrapped LL test statistic
    JnLLb.r_vct <- c(JnLLb.r_vct, min(S.evals))
  }
  
  # Obtain the (1 - \alpha)-quantile of the bootstrap distribution
  cvLLn <- quantile(JnLLb.r_vct, probs = alpha)
  cvLLn
}

#### Analyzing the results ####

#' @title Obtain identified set based on results of main estimation algorithm.
#' 
#' @description
#' Takes the results of the main estimation algorithm as input and outputs the 
#' 1 dimensional identified set.
#' 
#' @param test.results Results of main algorithm.
#' 
#' @noRd
#' 
get.identified.set <- function(test.results) {
  
  # If no feasible points were found, return [-\infty, \infty]
  if (length(which(test.results[, 2] <= test.results[, 3])) == 0) {
    return(c(-Inf, Inf))
  }
  
  # This is precisely step 6 in the algorithm described in Bei, 2024
  lb <- min(test.results[which(test.results[, 2] <= test.results[, 3]), 1])
  ub <- max(test.results[which(test.results[, 2] <= test.results[, 3]), 1])
  
  # Return results
  c(lb, ub)
}

