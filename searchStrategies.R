
#### Import dependencies ####

require(SPOT)

#### Helper functions ####

#' @title Obtain next point for feasible point search.
#' 
#' @description
#' Function to obtain the next point to evaluate in the search for a feasible
#' point. This function evaluates the entire parameter space of the component of
#' theta as evenly as possible. Used in the initialization step
#' (feasible_point_search.R)
#' 
#' @param evaluations Matrix of evaluations performed so far.
#' @param lb.theta Lower bound on the parameter of interest.
#' @param ub.theta Upper bound on the parameter of interest.
#' 
#' @returns Next point in the feasible point search.
#' 
get.next.point <- function(evaluations, lb.theta, ub.theta) {
  
  # Determine whether a next evaluation point is needed. I.e., whether there
  # already is a feasible point in the given set of evaluated points.
  if (any(evaluations[, "t.stat"] <= evaluations[, "crit.val"])) {
    return(list(theta.next = NULL, idx.after = NULL, stop = TRUE))
  }
  
  # If not, get the next point to evaluate.
  
  # To start, append the bounds of the theta values to the vector of evaluated
  # theta points, if necessary. Also store flag of whether bounds were added
  lb.added <- FALSE
  ub.added <- FALSE
  eval.points.extended <- evaluations[, "theta"]
  if (!(lb.theta %in% evaluations[, "theta"])) {
    eval.points.extended <- c(lb.theta, eval.points.extended)
    lb.added <- TRUE
  }
  if (!(ub.theta %in% evaluations[, "theta"])) {
    eval.points.extended <- c(eval.points.extended, ub.theta)
    ub.added <- TRUE
  }
  
  # Obtain the lengths of in the intervals between two consecutive points
  int.len <- diff(eval.points.extended)
  
  # Obtain the indices of the unique elements
  idx.uniques <- which.unique(matrix(int.len, nrow = length(int.len)), tol = 1e-10)
  
  # If all elements in int.len are the same (and hence length(idx.uniques) = 1),
  # return the midpoint between lb.theta and the smallest evaluated point.
  if (length(idx.uniques) == 1) {
    return(list(theta.next = sum(eval.points.extended[1:2])/2,
                idx.after = ifelse(lb.added, 0, 1),
                stop = FALSE))
  }
  
  # If not, length(idx.uniques) will necessarily be 2. Return the point obtained
  # by adding the smallest of the two interval lengths to the end point of the
  # last smallest interval.
  return(list(theta.next = eval.points.extended[idx.uniques[1] + 1] + int.len[1],
              idx.after = idx.uniques[1] + ifelse(lb.added, 0, 1),
              stop = FALSE))
}

#' @title Insert row into a matrix at a given row index
#' 
#' @description
#' Used in initalization step (feasible_point_search.R).
#' 
#' @param evaluations Matrix of violation function evaluations.
#' @param row Row (evaluations) to be added to the evaluation matrix.
#' @param idx.after Index of the row of \code{evaluations} after which the given
#' row should be placed.
#' 
#' @returns Evaluation matrix.
#' 
insert.row <- function(evaluations, row, idx.after) {
  
  if (idx.after == 0) {
    evaluations <- rbind(row, evaluations)
  } else if (idx.after == nrow(evaluations)) {
    evaluations <- rbind(evaluations, row)
  } else {
    evaluations <- rbind(evaluations[1:idx.after,], row, evaluations[(idx.after + 1):nrow(evaluations), ])
  }
  
  evaluations
}

#' @title Expected improvement
#'
#' @description Used in the M-step (M_step.R). Note: predict(fit.krige, ...)
#' has weird beheviour when making predictions for a single value in terms of
#' standard error. We work around this issue in this implementation.
#' 
#' @param theta Vector of coefficients.
#' @param test.fun Test function (cf. \code{EstimationAlgorithmBei.R}).
#' @param fit.krige Fitted Kriging model.
#' @param theta.hash Tentative optimal value for theta, i.e., the largest or
#' smallest feasible value for theta (if dir = 1 or dir = -1, respectively). A
#' 'feasible value' is one that satisfies all moment restrictions.
#' @param dir Search direction. \code{dir = 1} corresponds to looking for an
#' upper bound. \code{dir = -1} corresponds to looking for a lower bound.
#' 
#' @returns The expected improvement.
#' 
#' @import SPOT stats
#' 
EI <- function(theta, test.fun, fit.krige, theta.hash, dir) {
  
  # Select another point for which to make the prediction (see note in the
  # description section of the documentation above).
  other.point <- theta + 1
  
  # Predicted (test statistic - critical value) based on kriging model
  pred.kriging <- predict(fit.krige, matrix(c(theta, other.point), nrow = 2))
  violation.theta <- pred.kriging[1, "y"]
  
  # Standard deviation of prediction
  sL.theta <- sqrt(pred.kriging[1, "s"])
  
  # Expected improvement
  dir * (theta - theta.hash) * (1 - pnorm(violation.theta/sL.theta))
}

#' @title Draw initial set of starting values for optimizing the expected
#' improvement.
#'
#' @description  Used in the M-step (get.starting.values.R). ToDo: Adapt this
#' code so as to also perform sample space contractions as in the MatLab
#' implementation of Bei (2024).
#' 
#' @param theta.hash Tentative optimal value for theta, i.e., the largest or
#' smallest feasible value for theta (if dir = 1 or dir = -1, respectively). A
#' 'feasible value' is one that satisfies all moment restrictions.
#' @param dir Search direction. \code{dir = 1} corresponds to looking for an
#' upper bound. \code{dir = -1} corresponds to looking for a lower bound.
#' @param hyperparams List of hyperparameters.
#' @param EI.fun Function used to compute the expected improvement. See also
#' \code{EI}.
#' 
#' @returns Initial set of starting values.
#' 
#' @import stats
#' 
#' @references Bei, X. (2024). Local linearieation based subvector inference in
#' moment inequality models. Journal of Econometrics. 238:105549-
#'
draw.sv.init <- function(theta.hash, dir, hyperparams, EI.fun) {
  
  #### Extract the necessary sampling hyperparameters ####
  
  # Minimum and maximum distance of sampled points from current best theta 
  min.dist <- hyperparams[["min.dist"]]
  max.dist <- hyperparams[["max.dist"]]
  
  # Bounds of the parameter space for theta
  theta.lb <- hyperparams[["theta.lb"]]
  theta.ub <- hyperparams[["theta.ub"]]
  
  # Total number of sampled points required in initial drawing process
  nbr.init.sample.points <- hyperparams[["nbr.init.sample.points"]]
  
  # Number of points sampled per iteration in the initial drawing process
  nbr.points.per.iter.init <- hyperparams[["nbr.points.per.iter.init"]]
  
  # Total number of uniformly drawn points in the initial set of starting
  # values
  nbr.init.unif <- hyperparams[["nbr.init.unif"]]
  
  #### Draw points from a box ####
  
  # Initialize matrix that will store all the drawn points and their expected
  # improvement.
  init.draws <- matrix(nrow = 0, ncol = 2)
  colnames(init.draws) <- c("init.val", "EI")
  
  # Iterate until a prespecified number of points have been drawn or the
  # minimum distance of sampled points from theta.hash exceeds the maximum
  # distance (at the end of each iteration, the maximum distance is halved).
  while ((nrow(init.draws) < nbr.init.sample.points) & (min.dist < max.dist)) {
    
    # Compute bounds in which to draw points [See ToDo].
    if (dir == 1) {
      lb.init.sample.space <- theta.hash
      ub.init.sample.space <- min(theta.hash + max.dist, theta.ub)
    } else if (dir == -1) {
      lb.init.sample.space <- max(theta.hash - max.dist, theta.lb)
      ub.init.sample.space <- theta.hash
    }
    
    # Sample points in [lb.init.sample.space, ub.init.sample.space]
    draws <- runif(nbr.points.per.iter.init, min = lb.init.sample.space,
                   max = ub.init.sample.space)
    
    # For each draw, compute the expected improvement. Store it.
    for (draw in draws) {
      init.draws <- rbind(init.draws, c(draw, EI.fun(draw)))
    }
    
    # Only keep the draws with positive expected improvement
    keep.idxs <- which(init.draws[, "EI"] > 1e-10)
    init.draws <- init.draws[keep.idxs, , drop = FALSE]
    
    # Reduce the maximum distance from theta.hash
    max.dist <- max.dist / 2
  }
  
  # Only keep the top 'nbr.init.sample.points' points with highest expected
  # improvement, taking into account that init.draws can still be empty.
  if (nrow(init.draws) > 0) {
    init.draws <- init.draws[order(init.draws[,2], decreasing = TRUE), , drop = FALSE]
    init.draws <- init.draws[1:min(nbr.init.sample.points, nrow(init.draws)), , drop = FALSE]
  }
  
  #### Draw points uniformly from the parameter space ####
  
  # Draw from uniform distribution
  unif.draws <- runif(nbr.init.unif, min = theta.lb, max = theta.ub)
  
  # Set the expected improvement of the uniformly drawn points equal to zero.
  # (This might not be true but we do not require EI to be positive for these
  # points, so there is no need to compute it)
  unif.draws.EI <- rep(0, nbr.init.unif)
  
  # Append to other starting values
  init.draws <- rbind(init.draws, cbind(unif.draws, unif.draws.EI))
  
  #### Return the results ####
  
  init.draws
}

#' @title Analogue to KMS_AUX4_MSpoints(...) in MATLAB code of Bei (2024).
#'
#' @description Create starting values for EI maximization. Used in the M-step
#' (get.starting.values.R).
#' 
#' @param draws.init Initial draws.
#' 
#' @references Bei, X. (2024). Local linearieation based subvector inference in
#' moment inequality models. Journal of Econometrics. 238:105549-
#' 
MSpoint <- function(draws.init) {
  X <- draws.init[, 1]
  n.eval <- length(X)
  
  # Initialize object that will store ms_points
  ms_points <- matrix(0, nrow = n.eval, ncol = 2)
  
  # Sort the starting values
  temp <- sort(X)
  
  # Replace the k-th coordinate with averages with nearest neighbours
  temp1 <- c((temp[1] + temp[2])/2, (temp[-1] + temp[-length(temp)])/2)
  temp2 <- c(temp1[-1], (temp[length(temp) - 1] + temp[length(temp)])/2)
  sel <- rbinom(n.eval, 1, 0.5)
  
  # Randomize between candidates
  ms_points[ , 1] <- sel * temp1 + (1 - sel) * temp2
  
  # Use the ones that were not selected for the second dimension
  ms_points[ , 2] <- (1 - sel) * temp1 + sel * temp2
  
  # Reshape and return the evaluation points
  c(ms_points[, 1], ms_points[, 2])
}

#' @title Main function for obtaining the starting values of the expected
#' improvement maximization step.
#' 
#' @description Obtain starting values used in the M-step (M_step.R).
#' 
#' @param theta.hash Tentative optimal value for theta, i.e., the largest or
#' smallest feasible value for theta (if dir = 1 or dir = -1, respectively). A
#' 'feasible value' is one that satisfies all moment restrictions.
#' @param dir Search direction. \code{dir = 1} corresponds to looking for an
#' upper bound. \code{dir = -1} corresponds to looking for a lower bound.
#' @param EI.Mstep Function to compute expected improvements.
#' @param hyperparams List of hyperparameters.
#' 
#' @returns 
get.starting.values <- function(theta.hash, dir, EI.Mstep, hyperparams) {
  
  # Extract necessary hyperparameters for readibility
  nbr.start.vals <- hyperparams[["nbr.start.vals"]]
  
  # Draw initial set of starting values
  draws.init <- draw.sv.init(theta.hash, dir, hyperparams, EI.Mstep)
  
  # Append theta.hash and modify them with the function MSpoints.R
  draws.init <- rbind(draws.init, c(theta.hash, 0))
  start.vals <- MSpoint(draws.init)
  
  # Only keep unique starting values
  start.vals <- start.vals[which.unique(matrix(start.vals, ncol = 1))]
  
  # Compute the expected improvement for each starting value
  EI.start.vals <- rep(0, length(start.vals))
  for (i in 1:length(start.vals)) {
    EI.start.vals[i] <- EI.Mstep(start.vals[i])
  }
  
  # Keep top nbr.start.vals starting values with highest expected improvement
  idxs.largest <- tail(order(EI.start.vals), nbr.start.vals)
  start.vals <- start.vals[idxs.largest]
  
  # Also include theta.hash and a slight perturbation of it
  theta.eps <- theta.hash + 1e-4
  start.vals <- c(start.vals, theta.hash, theta.eps)
  
  # Return the results
  start.vals
}

#' @title Optimize the expected improvement
#' 
#' @description This function finds the point for which the expected improvement
#' is optimal, based on a given set of starting values. (M_step.R)
#' 
#' @param start.vals Starting values for optimization.
#' @param EI.Mstep Function to compute expected improvements.
#' @param hyperparams List of hyperparameters.
#' 
#' @returns Maximum of the expected imrpovement function.
#' 
#' @import stats
#' 
do.optimization.Mstep <- function(start.vals, EI.Mstep, hyperparams) {
  
  # Bounds of the parameter space for theta
  theta.lb <- hyperparams[["theta.lb"]]
  theta.ub <- hyperparams[["theta.ub"]]
  
  # Number of optimal values to return
  nbr.opt.EI <- hyperparams[["nbr.opt.EI"]]
  
  # Initialize object that will store the results of the optimization
  opt.EI <- rep(0, length(start.vals))
  opt.theta <- rep(0, length(start.vals))
  
  # Starting from each starting value, optimize the expected improvement
  for (i in 1:length(start.vals)) {
    
    # Set the starting value of this iteration
    start.val <- start.vals[i]
    
    # Maximize the expected improvement function
    opt.out <- optim(start.val, EI.Mstep, method = "L-BFGS-B", lower = theta.lb,
                     upper = theta.ub, control = list(fnscale = -1))
    opt.EI[i] <- opt.out$value
    opt.theta[i] <- opt.out$par
  }
  
  # Keep the top nbr.top.EI unique points with the highest expected improvement
  idxs.unique <- which.unique(matrix(opt.theta, ncol = 1))
  opt.EI <- opt.EI[idxs.unique]
  opt.theta <- opt.theta[idxs.unique]
  
  idxs.keep <- tail(order(opt.EI), nbr.opt.EI)
  opt.EI <- opt.EI[idxs.keep]
  opt.theta <- opt.theta[idxs.keep]
  
  # Return the results
  cbind(opt.theta, opt.EI)
}

#' @title Get extra evaluation points for E-step
#' 
#' @description Function used to obtain extra theta values to be supplied to the
#' E-step in the next iteration (M_step.R). Note: this function should be 
#' changed when implementing the sample space contractions (see comment made in
#' documentation of \code{M_step}).
#' 
#' @param dir Search direction. \code{dir = 1} corresponds to looking for an
#' upper bound. \code{dir = -1} corresponds to looking for a lower bound.
#' @param theta.hash Tentative optimal value for theta, i.e., the largest or
#' smallest feasible value for theta (if dir = 1 or dir = -1, respectively). A
#' 'feasible value' is one that satisfies all moment restrictions.
#' @param maxviol.hash Violation curve evaluated at  \code{theta.hash}.
#' @param hyperparams List of hyperparameters.
#' 
#' @returns Points to evaluate in E-step.
#' 
#' @import stats
#' 
get.extra.Estep.points <- function(dir, theta.hash, maxviol.hash,
                                   hyperparams) {
  
  # Extract the necessary information
  theta.lb <- hyperparams[["theta.lb"]]
  theta.ub <- hyperparams[["theta.ub"]]
  nbr.extra <- hyperparams[["nbr.extra"]]
  EAM_thetadistort <- hyperparams[["min.improvement"]]
  
  # Randomly sampled points
  rsp <- runif(nbr.extra, min = theta.lb, max = theta.ub)
  
  # Initialize object that will store perturbed theta.hash values
  theta_eps <- rep(0, 3)
  
  # Slight perturbations of theta.hash
  delta1 <- abs(maxviol.hash)
  delta2 <- EAM_thetadistort
  delta3 <- 10*EAM_thetadistort
  if (dir == 1) {
    theta_eps[1] <- min(theta.hash + delta1, theta.ub)
    theta_eps[2] <- min(theta.hash + delta2, theta.ub)  
    theta_eps[3] <- min(theta.hash + delta3, theta.ub)  
  } else if (dir == -1) {
    theta_eps[1] <- max(theta.hash - delta1, theta.lb)
    theta_eps[2] <- max(theta.hash - delta2, theta.lb)  
    theta_eps[3] <- max(theta.hash - delta3, theta.lb)  
  }
  
  # Return results
  c(rsp, theta_eps)
}

#### Steps in EAM algorithm ####

#' @title Method for finding initial points of the EAM algorithm
#' 
#' @description Also called the 'initialization' step in KMS19, this method
#' tries to find at least one initial feasible point, which is required to run
#' the EAM algorithm.
#' ToDo: Investigate whether the feasible point search of Bei (2024) is better.
#' If so, implement it.
#' 
#' @param test.fun Function that takes a parameter vector as a first argument
#' and returns the test statistic, as well as the critical value.
#' @param lb.theta Lower bound on theta.
#' @param ub.theta Upper bound on theta
#' @param min.eval Minimum number of theta values at which to evaluate the test
#' statistic and critical value.
#' @param max.eval Maximum number of theta values at which to evaluate the test
#' statistic and critical value.
#' @param picturose Picturosity flag. If \code{TRUE}, a plot illustrating the
#' workings of the algorithm will updated during runtime. Default is
#' \code{picturose = FALSE}.
#' @param parallel Flag for whether or not parallel computing should be used.
#' Default is \code{parallel = FALSE}.
#' 
#' @returns Results of the initial feasible point search.
#' 
#' @refernces Kaido et al. (2019). Confidence intervals for projections of
#' partially identified parameters. Econometrica. 87(4):1397-1432.
#' 
feasible_point_search <- function(test.fun, hyperparams, verbose,
                                  picturose = FALSE, parallel = FALSE) {
  
  #### Precondition checks ####
  
  # For ease of readability
  lb.theta <- hyperparams[["theta.lb"]]
  ub.theta <- hyperparams[["theta.ub"]]
  min.eval <- hyperparams[["min.eval"]]
  max.eval <- hyperparams[["max.eval"]]
  
  # Bounds on theta should be finite
  if ((abs(lb.theta) == Inf) | (abs(ub.theta) == Inf)) {
    stop("Bound on theta must be finite.")
  }
  
  # Minimum number of points to evaluate should be strictly positive
  if (min.eval <= 0) {
    stop("Minimum number of points to evaluate must be strictly positive.")
  }
  
  # Minimum number of evaluations should be smaller than maximum number of
  # evaluations
  if (min.eval > max.eval) {
    stop("min.eval should be smaller than max.eval.")
  }
  
  #### Initial points to evaluate. If feasible point is found, return ####
  
  # Define initial set of points to evaluate
  pte <- seq(lb.theta, ub.theta, length.out = min.eval)
  
  # If plot of estimation procedure should be drawn, do so
  if (picturose) {
    plot.addpte(pte, col = "white")
  }
  
  # Initialize object that will store the test statistics and critical value for
  # each point in this initial set of points to evaluate.
  evaluations <- matrix(nrow = min.eval, ncol = 3)
  colnames(evaluations) <- c("theta", "t.stat", "crit.val")
  
  # For each point, obtain the test statistic and critical value.
  if (parallel) { # Using parallel computing
    
    # Update the user. Add to plot.
    if (verbose >= 3) {
      message("Evaluating initial points in parallel")
    }
    
    # Compute batch size and batch indices
    batch.size <- length(clust)
    n.batchs <- ceiling(length(pte)/batch.size)
    for (batch.idx in 1:n.batchs) {
      
      # Compute indices of pte vector to handle in this batch
      i.start <- (batch.idx - 1)*batch.size + 1
      i.end <- min(batch.idx * batch.size, length(pte))
      
      # If required, add points to plot
      if (picturose) {
        plot.addpte(pte[i.start:i.end])
      }
      
      # Evaluate batch of points.
      suppressWarnings({evaluations[i.start:i.end, ] <-
        foreach(i = i.start:i.end, .combine = 'rbind',
                .export = c("par.space", "hp", "c", "inst.func.evals", "t",
                            "data", "test.fun", "options")
                ) %dopar% {
        
        # Load all necessary packages
        source("simulationFunctions.R")
        
        # Select theta of this iteration
        theta <- pte[i]
        
        # Run the test
        test.out <- test.fun(theta)
        
        # Return the results
        c(theta, test.out[["t.stat"]], test.out[["crit.val"]])
      }})
      
      # If required, add points to plot
      if (picturose) {
        plot.addpte.eval(evaluations[i.start:i.end, , drop = FALSE])
      }
    }
    
  } else { # Using sequential computing
    
    for (i in 1:min.eval) {
      if (verbose >= 3) {
        message(sprintf("Checking initial points (%s / %s)", i, min.eval))
      }
      
      # Select theta of this iteration
      theta <- pte[i]
      
      # If necessary, plot it
      if (picturose) {
        plot.addpte(theta)
      }
      
      # Run the test
      test.out <- test.fun(theta)
      
      # Store the results
      evaluations[i,] <- c(theta, test.out[["t.stat"]], test.out[["crit.val"]])
      
      # If plot of estimation procedure should be drawn, do so.
      if (picturose) {
        plot.addpte.eval(evaluations[i, , drop = FALSE])
      }
    }
  }
  
  # If there is at least one feasible point in the set of evaluated points,
  # return the set of evaluated points.
  if (any(evaluations[, "t.stat"] <= evaluations[, "crit.val"])) {
    return(list(evaluations = evaluations))
  }
  
  #### If no feasible point was found, continue search ####
  
  if (parallel) { # Using parallel computing
    
    # Set some parameters
    batch.size <- length(clust)
    eval.nbr <- min.eval
    stop <- FALSE
    
    # Evaluate (in batch) additional points until stopping criterion is reached
    while(eval.nbr < max.eval & !stop) {
      
      # Obtain next batch of points to evaluate
      evaluations.dummy <- evaluations
      pte <- c()
      pte.idxafter <- c()
      for (i in 1:batch.size) {
        gnp.out <- get.next.point(evaluations.dummy, lb.theta, ub.theta)
        theta.next <- gnp.out[["theta.next"]]
        idx.after <- gnp.out[["idx.after"]]
        row <- c(theta.next, 1, 0)
        evaluations.dummy <- insert.row(evaluations.dummy, row, idx.after)
        pte <- c(pte, theta.next)
        pte.idxafter <- c(pte.idxafter, idx.after)
      }
      
      # Evaluate batch
      evaluations.add <-
        foreach(i = 1:batch.size, .combine = 'rbind') %dopar% {
          
          # Load all necessary packages
          source("simulationFunctions.R")
          
          # Select theta of this iteration
          theta <- pte[i]
          
          # Run the test
          test.out <- test.fun(theta)
          
          # Return the results
          c(theta, test.out[["t.stat"]], test.out[["crit.val"]])
        }
      
      # Store the results
      for (i in 1:batch.size) {
        evaluations <- insert.row(evaluations, evaluations.add[i, ], pte.idxafter[i])
      }
      
      # Update stopping criteria
      stop <- any(evaluations[, "t.stat"] < evaluations[, "crit.val"])
      eval.nbr <- eval.nbr + batch.size
    }
    
  } else { # Using sequential computing
    
    # Initialization for while-loop
    gnp.out <- get.next.point(evaluations, lb.theta, ub.theta)
    theta.next <- gnp.out[["theta.next"]]
    idx.after <- gnp.out[["idx.after"]]
    stop <- gnp.out[["stop"]]
    eval.nbr <- min.eval + 1
    
    # If no feasible point has been found yet, test if the midpoints
    # between evaluated points are feasible. Continue this procedure until either
    # a feasible point is found, or the the stopping criterion is reached.
    while (!stop & (eval.nbr <= max.eval)) {
      
      # Update user
      if (verbose >= 3) {
        message(sprintf("Checking additional points (%s / %s)", eval.nbr, max.eval))
      }
      if (picturose) {
        plot.addpte(theta.next)
      }
      
      # Run the test for this point
      test.out <- test.fun(theta.next)
      
      # Add the results to the set of evaluated points
      row <- c(theta.next, test.out[["t.stat"]], test.out[["crit.val"]])
      evaluations <- insert.row(evaluations, row, idx.after)
      
      # Get next point
      gnp.out <- get.next.point(evaluations, lb.theta, ub.theta)
      theta.next <- gnp.out[["theta.next"]]
      idx.after <- gnp.out[["idx.after"]]
      stop <- gnp.out[["stop"]]
      
      # Increment number of evaluations
      eval.nbr <- eval.nbr + 1
      
      # Update plot
      if (picturose) {
        plot.addpte.eval(theta, test.out)
      }
    }
    
  }
  
  #### Return the results ####
  
  list(evaluations = evaluations)
}

#' @title E-step in the EAM algorithm as described in KMS19.
#' 
#' @description This function performs the estimation step in the EAM algorithm.
#' 
#' @param thetas Points at which to perform the E-step. Usually the result of
#' the M-step.
#' @param test.fun Function returning the test statistic, as well as the critical
#' value.
#' @param dir Direction in which to optimize. For finding upper bounds, set
#' \code{dir = 1}, for finding lower bounds, set \code{dir = -1}.
#' @param evaluations Matrix containing each point that was already evaluated,
#' alongside the corresponding test statistic and critical value, as its rows.
#' @param Verbose Verbosity parameter.
#' 
#' @returns Results of the E-step.
#' 
E_step <- function(thetas, test.fun, dir, evaluations, verbose) {
  
  # Loop over all values of theta to be checked
  for (i in 1:length(thetas)) {
    
    # Get the theta value of this iteration
    theta <- thetas[i]
    
    # Check whether the evaluation of the current theta value was already
    # carried out and can hence be skipped
    if (any(abs(theta - evaluations[, "theta"]) < 1e-10)) {
      if (verbose >= 3) {
        message(sprintf(
          "\t Evaluating point %s out of %s... (skipped)", i, length(thetas)
          ))
      }
      next
    }
    
    if (verbose >= 3) {
      message(sprintf("\t Evaluating point %s out of %s...", i, length(thetas)))
    }
    
    # Obtain the test statistic and critical value of the given point
    test.out <- test.fun(theta)
    t.stat <- test.out[["t.stat"]]
    crit.val <- test.out[["crit.val"]]
    
    # Append the point to the set of evaluated points
    evaluations <- rbind(evaluations, c(theta, t.stat, crit.val))
    evaluations <- evaluations[order(evaluations[, 1]),]
  }
  
  # Indices of feasible points
  feas.idxs <- which(evaluations[, "t.stat"] <= evaluations[, "crit.val"])
  
  # Tentative optimal value
  if (dir == 1) {
    theta.hash <- max(evaluations[feas.idxs, "theta"])
  } else if (dir == -1) {
    theta.hash <- min(evaluations[feas.idxs, "theta"])
  }
  
  # Return the results
  list(evaluations = evaluations, theta.hash = theta.hash)
}

#' @title A-step in the EAM algorithm described in KMS19
#' 
#' @description This function performs the approximation step in the EAM
#' algorithm. More specifically, it fits a Gaussian-process regression model
#' (Kriging) to the evaluated data points (\theta, c(\theta)).
#' 
#' @param evaluations Matrix containing each point that was already evaluated,
#' alongside the corresponding test statistic and critical value, as its rows.
#' @param verbose Verosity parameter.
#' 
#' @returns Results of the A-step.
#' 
#' @import SPOT
#' 
A_step <- function(evaluations, verbose = 0) {
  
  # Matrix of theta values
  theta.mat <- matrix(evaluations[, "theta"], nrow = nrow(evaluations))
  
  # Make sure design points are not too close together, as it will cause 
  # unwanted behaviour of the kriging model.
  idxs.unique <- which.unique(theta.mat, tol = 1/nrow(theta.mat))
  theta.mat <- theta.mat[idxs.unique, , drop = FALSE]
  
  # Vector of violations (test statistic - critical value)
  violations.vct <- evaluations[idxs.unique, "t.stat"] - evaluations[idxs.unique, "crit.val"]
  
  # Control parameters
  control = list(regr = regpoly0, corr = corrkriging, target = c("y", "s"))
  
  # Fit the Kriging model
  fit.krige <- buildKrigingDACE(theta.mat, violations.vct, control = control)
  
  # If asked, plot the Kriging model
  if (verbose >= 3) {
    x.vals <- seq(min(evaluations[, "theta"]), max(evaluations[, "theta"]),
                  length.out = 500)
    predictions <- predict(fit.krige, matrix(x.vals, nrow = length(x.vals)))
    y.vals <- predictions$y
    sd.vals <- predictions$s
    
    plot(x.vals, y.vals, type = 'l', xlab = "theta", ylab = "predicted violation",
         main = "Kriging model")
    lines(x.vals, y.vals + 2 * sd.vals, type = 'l', lty = 2, col = "red")
    lines(x.vals, y.vals - 2 * sd.vals, type = 'l', lty = 2, col = "red")
    abline(h = 0, col = "blue", lty = 2)
  }
  
  # Return the Kriging model
  fit.krige
}

#' @title M-step in the EAM algorithm described in KMS19.
#' 
#' @description  This function performs the maximization step in the EAM
#' algorithm. More specifically, it maximizes the expected improvement.
#' ToDo: implement sample space contractions (see comment made in documentation
#' of \code{draw.sv.init}).
#' 
#' @param dir Direction to search in. \code{dir = 1} corresponds to finding the
#' upper bound of the confidence interval. \code{dir = -1} corresponds to
#' finding the lower bound.
#' @param evaluations Matrix containing each point that was already evaluated,
#' alongside the corresponding test statistic and critical value, as its rows.
#' @param theta.hash Tentative best value of theta. Obtained from the E-step.
#' @param fit.krige Kriging model obtained from the A-step.
#' @param c Projection vector.
#' @param par.space Bounds of the parameter space.
#' @param hyperparams Parameters used in obtaining initial values
#' for the maximization algorithm. If \code{NULL}, default values are used.
#' Default is \code{hyperparams = NULL}.
#' @param verbose Verbosity parameter.
#' 
M_step <- function(dir, evaluations, theta.hash, fit.krige, c, par.space, 
                   hyperparams, verbose) {
  
  #### Set some useful variables ####
  
  # Compute tentative optimal value for theta if not supplied
  if (is.null(theta.hash)) {
    
    # Indices of feasible points
    feas.idxs <- which(evaluations[, "t.stat"] <= evaluations[, "crit.val"])
    
    # Find current best value
    if (dir == 1) {
      theta.hash <- max(evaluations[feas.idxs, "theta"])
    } else if (dir == -1) {
      theta.hash <- min(evaluations[feas.idxs, "theta"])
    }
  }
  
  # Current largest value of feasible theta
  theta.hash.idx <- which(evaluations[, "theta"] == theta.hash)[1]
  maxviol.hash <- evaluations[theta.hash.idx, "t.stat"] - evaluations[theta.hash.idx, "crit.val"]
  
  # Define expected improvement function wrt theta.hash
  EI.Mstep <- function(theta) {EI(theta, test.fun, fit.krige, theta.hash, dir)}
  
  #### Main code for M-step ####
  
  # Obtain starting values for the maximization of the EI
  start.vals <- get.starting.values(theta.hash, dir, EI.Mstep, hyperparams)
  
  # Compute the optimal expected improvement and corresponding theta-value.
  opt.res <- do.optimization.Mstep(start.vals, EI.Mstep, hyperparams)
  
  # Add a randomly drawn theta s.t. dir * theta > dir * theta.hash to the set
  # of points to evaluate in the next E-step.
  extra.points <- get.extra.Estep.points(dir, theta.hash, maxviol.hash,
                                         hyperparams)
  
  # If asked, plot the expected improvement function, theta.hash and the most
  # promising theta value.
  if (verbose >= 3) {
    x.vals <- seq(min(evaluations[, "theta"]), max(evaluations[, "theta"]),
                  length.out = 500)
    y.vals <- Vectorize(EI.Mstep)(x.vals)
    
    plot(x.vals, y.vals, type = 'l', xlab = "theta", ylab = "EI",
         main = "Expected improvement function M-step")
    abline(v = theta.hash, col = "red")
    abline(v = opt.res[1], col = "green")
  }
  
  # Return the results
  list(opt.res = opt.res, extra.points = extra.points)
}

#' @title Check convergence of the EAM algorithm.
#'
#' @description This function checks the convergence of the EAM algorithm.
#' ToDo: Get rid of hard coding stop of algorithm when no more improvement of
#' theta values (maybe related to parameter space contraction, since the problem
#' is that the given points to check in the E-step of the following iteration
#' can always be the same and always be rejected (except of course for the
#' randomly chosen one), while the most promising theta value continues to be
#' the same, infeasible value. In this way, it is possible that
#' theta.hash - mp.theta.next at some point will never decrease).
#' 
#' @param opt.val.prev Previous optimal theta value.
#' @param evaluations Matrix of violation curve evaluations.
#' @param mp.theta.next Most promising value of theta for which to run the
#' E-step in the following iteration
#' @param iter.nbr Number of iterations of the EAM algorithm run so far.
#' @param dir Search direction.
#' @param hyperparams List of hyperparameters used in the EAM algorithm.
#' @param verbose Verbosity parameter.
#' 
#' @returns Boolean value whether or not algorithm has converged.
#' 
EAM.converged <- function(opt.val.prev, evaluations, mp.theta.next, iter.nbr,
                          dir, hyperparams, verbose) {
  
  #### Extract necessary information from set of hyperparameters ####
  
  # Minimum amount that this step should have improved upon the previous step
  min.improvement <- hyperparams[["min.improvement"]]
  
  # Minimum amount of improvement that can be the result of running the next
  # step.
  min.possible.improvement <- hyperparams[["min.possible.improvement"]]
  
  # Minimum number of EAM iterations that should be run.
  EAM.min.iter <- hyperparams[["EAM.min.iter"]]
  
  # Maximum number of EAM iterations that should be run.
  max.iter <- hyperparams[["max.iter"]]
  
  #### Check all convergence criteria ####
  
  # Indices of feasible points
  feas.idxs <- which(evaluations[, "t.stat"] <= evaluations[, "crit.val"])
  
  # Best and second best feasible point
  if (dir == 1) {
    opt.val <- sort(evaluations[feas.idxs, "theta"], decreasing = TRUE)[1]
  } else if (dir == -1) {
    opt.val <- sort(evaluations[feas.idxs, "theta"], decreasing = FALSE)[1]
  }
  
  # If this is run in the first iteration, opt.val.prev = NULL. In this case,
  # set opt.val.prev = -Inf
  if (is.null(opt.val.prev)) {
    opt.val.prev <- -Inf
  }
  
  conv1 <- (dir * (opt.val - opt.val.prev) < min.improvement)
  conv2 <- (dir * (mp.theta.next - opt.val) < min.possible.improvement)
  conv3 <- (iter.nbr > EAM.min.iter)
  
  # Determine whether or not to stop the algorithm (stopping criterion based on
  # implementation of Bei, 2024).
  stop <- (conv1 & conv2 & conv3) | (iter.nbr > max.iter)
  
  # If there is no more improvement in theta value, stop the algorithm. This is
  # valid, since each iteration tries to improve the current best theta value by
  # 'min.improvement'.
  # 
  # !!! See also [ToDo] !!!
  if (abs(opt.val - opt.val.prev) < 1e-12) {
    stop <- TRUE
  }
  
  # If necessary, print convergence information to the console
  if (verbose >= 3) {
    message(
      "_______________________________________________________________________________")
    message(sprintf("Iteration %s:\n", iter.nbr))
    message(
      "Impr. wrt previous | Possible improvement | min.iter reached | max.iter reached"
      )
    message(
      "-------------------|----------------------|------------------|-----------------"
    )
    message(sprintf(
      "%18f | %20f | %16s | %15s", dir * (opt.val - opt.val.prev),
      dir * (mp.theta.next - opt.val), ifelse(conv3, "TRUE", "FALSE"),
      ifelse(iter.nbr > max.iter, "TRUE", "FALSE")
    ))
    message(
      "_______________________________________________________________________________")
  } else if (verbose %in% c(1, 2)) {
    if (iter.nbr == 1) {
      message(
        "______________________________________________________________________________________"
        )
      message(
        "iter | Impr. wrt previous | Possible improvement | min.iter reached | max.iter reached"
      )
      message(
        "-----|--------------------|----------------------|------------------|-----------------"
      )
    }
    message(sprintf(
      "%4d | %18f | %20f | %16s | %15s", iter.nbr, dir * (opt.val - opt.val.prev),
      dir * (mp.theta.next - opt.val), ifelse(conv3, "TRUE", "FALSE"),
      ifelse(iter.nbr > max.iter, "TRUE", "FALSE")
    ))
  }
  
  # Return the result
  stop
}

#### Steps in grid search ####

#' @title Rudimentary, bidirectional 1D grid search algorithm.
#' 
#' @description
#' This function implements a rudimentary, bidirectional search algorithm. It 
#' works by expanding a grid with given step.size in both directions, starting
#' from an initial feasible point.
#' 
#' @param test.results Matrix containing the evaluations of the test statistic
#' and critical value.
#' @param max.iter Maximum number of iterations.
#' @param step.size Step size based on which the grid is constructed.
#' 
#' @returns The next point to evaluate in the grid search.
#' 
gs.algo.bidir <- function(test.results, max.iter, step.size) {
  
  # Define some useful variables
  n.test <- nrow(test.results)
  
  # Point around which the grid is centered
  center.point <- test.results[1, 1]
  
  # Initialize some variables
  stop <- FALSE
  
  #### Some edge cases ####
  
  # If test.results only contains the initial point, confirm that the test
  # did not lead to a rejection and if so, return a new r value to the right of
  # the original one. Else stop the algorithm
  if (n.test == 1) {
    if (test.results[n.test, 2] <= test.results[n.test, 3]) {
      r <- center.point + step.size
      return(list(r = r, stop = stop))
    } else {
      return(list(r = NULL, stop = TRUE))
    }
  }
  
  # If test.results contains the initial point and one initial point, return an
  # r evaluation point in the other direction of the initial point.
  if (n.test == 2) {
    r <- center.point - sign(test.results[n.test, 1] - center.point) * step.size
    return(list(r = r, stop = stop))
  }
  
  # If maximum number of grid evaluations is reached, stop the algorithm
  if (n.test >= max.iter) {
    return(list(r = NULL, stop = TRUE))
  }
  
  #### Main logic ####
  
  # If the test evaluations for the last two r values where at opposite sides of
  # the center point...
  if ((test.results[n.test, 1] - center.point) * (test.results[n.test - 1, 1] - center.point) < 0) {
    
    # If the second to last r value led to a non-rejection...
    if (test.results[n.test - 1, 2] <= test.results[n.test - 1, 3]) {
      
      # Return an r value one step further than this second to last r value
      r <- test.results[n.test - 1, 1] + sign(test.results[n.test - 1, 1]) * step.size
    }
    
    # If the second to last r value led to a rejection...
    else {
      
      # Check whether the last r value led to a rejection. If it didn't lead to
      # a rejection... 
      if (test.results[n.test, 2] <= test.results[n.test, 3]) {
        
        # Return an r value one step further than the last r value
        r <- test.results[n.test, 1] + sign(test.results[n.test, 1]) * step.size
      }
      
      # If it did lead to a rejection...
      else {
        
        # Stop the algorithm
        r <- NULL
        stop <- TRUE
        
      }
    }
  }
  
  # If the test evaluation for the last two r values was at the same side of the
  # center point...
  else {
    
    # Test whether the last r value evaluation led to a rejection. If it didn't...
    if (test.results[n.test, 2] <= test.results[n.test, 3]) {
      
      # Return an r value one step further than the last r value
      r <- test.results[n.test, 1] + sign(test.results[n.test, 1]) * step.size
    }
    
    # If if did...
    else {
      
      # Stop the algorithm
      r <- NULL
      stop <- TRUE
    }
  }
  
  # Return results
  return(list(r = r, stop = stop))
}

#' @title Return the next point to evaluate when doing regular grid search
#' 
#' @description This function implements a unidirectional grid search, that
#' works by expanding a grid starting from a given feasible point in the
#' given direction.
#' 
#' @param evaluations Matrix of evaluated test statistics and critical values.
#' @param dir Search direction.
#' @param iter.nbr Iteration number.
#' @param hp List of hyperparameters.
#' 
#' @returns Next point to evaluate in the search algorithm.
gs.regular <- function(evaluations, dir, iter.nbr, hp) {
  
  # Extract the necessary hyperparameters
  step.size <- hp[["step.size"]]
  max.iter <- hp[["max.iter"]]
  
  # Find current largest (smallest) feasible value, called theta.hash
  feas.idx <- which(evaluations[, "t.stat"] <= evaluations[, "crit.val"])
  theta.hash <- dir * max(dir * evaluations[feas.idx, "theta"])
  
  # If theta.hash is larger (smaller) than the upper (lower) bound, stop the
  # algorithm.
  if (((dir == 1) & (theta.hash > hp$theta.ub)) | 
      ((dir == -1) & (theta.hash < hp$theta.lb))) {
    theta.to.eval <- theta.hash + dir
    return(list(theta.to.eval = theta.to.eval, stop = TRUE))
  }
  
  # Get next theta value
  theta.to.eval <- theta.hash + dir * step.size
  
  # Check whether theta.to.eval has already been evaluated and if so, whether or
  # not it was feasible.
  stop <- FALSE
  if (theta.to.eval %in% evaluations[, "theta"]) {
    idx.tte <- which(evaluations[, "theta"] == theta.to.eval)
    stop <- evaluations[idx.tte, "t.stat"] > evaluations[idx.tte, "crit.val"]
  }
  if (iter.nbr > max.iter) {
    stop <- TRUE
  }
  
  # Return the results
  list(theta.to.eval = theta.to.eval, stop = stop)
}

#' @title Return the next point to evaluate when doing binary search
#' 
#' @description This function implements the binary search algorithm, that
#' starts from a given feasible point and looks in the given direction for the
#' root of the violation curve.
#' 
#' @param evaluations Matrix of evaluated test statistics and critical values.
#' @param dir Search direction.
#' @param iter.nbr Iteration number.
#' @param hp List of hyperparameters.
#' 
#' @returns The next point to evaluate.
#' 
gs.binary <- function(evaluations, dir, iter.nbr, hp) {
  
  # Extract the necessary hyperparameters
  bin.search.tol <- hp[["bin.search.tol"]]
  max.iter <- hp[["max.iter"]]
  
  # Check if maximum number of iterations has been reached. If so, stop the
  # algorithm
  if (iter.nbr > max.iter) {
    return(list(theta.to.eval = NULL, stop = TRUE))
  }
  
  # Find current largest (smallest) feasible value, called theta.hash
  feas.idx <- which(evaluations[, "t.stat"] <= evaluations[, "crit.val"])
  theta.hash <- dir * max(dir * evaluations[feas.idx, "theta"])
  
  # If theta.hash is larger (smaller) than or equal to the upper (lower) bound,
  # stop the algorithm.
  if (((dir == 1) & (theta.hash >= hp$theta.ub)) |
      ((dir == -1) & (theta.hash <= hp$theta.lb))) {
    return(list(theta.to.eval = NULL, stop = TRUE))
  }
  
  # Matrix of all infeasible points
  infeas.idxs <- which(evaluations[, "t.stat"] > evaluations[, "crit.val"])
  evals.infeas <- matrix(nrow = 0, ncol = 3)
  colnames(evals.infeas) <- colnames(evaluations)
  for (idx in infeas.idxs) {
    evals.infeas <- rbind(evals.infeas, evaluations[idx, ])
  }
  
  # Indices of points in infeasible point set that are larger (smaller) than
  # theta.hash
  idxs <- which(dir * evals.infeas[, "theta"] > dir * theta.hash)
  
  # If no infeasible point which is larger (smaller) than theta.hash is found
  # yet, search for one.
  #
  # NOTE: when both theta.lb and theta.ub are checked in an initial stage, this
  #       code block will never be ran. Indeed, either theta.lb (theta.ub) is
  #       feasible, in which case the root finding algorithm would already 
  #       have been stopped (checked above). If it is infeasible, then an
  #       infeasible point that is smaller (larger) than theta.hash is known.
  iter.nbr <- 1
  if (length(idxs) == 0) {
    
    # Obtain largest (smallest) and second largest (smallest) theta values
    out.sort <- dir * sort(dir * evaluations[, "theta"], decreasing = TRUE)[1:2]
    theta.largest <- out.sort[1]
    theta.second <- out.sort[2]
    
    # Theta value to evaluate
    dist <- 2 * max(abs(theta.largest - theta.second), 1)
    theta.to.eval <- dir * (dir * theta.largest + dist)
    
    # Return theta
    return(list(theta.to.eval = theta.to.eval, stop = FALSE))
  }
  
  # Else, determine theta.tilde
  theta.tilde <- dir * min(dir * evals.infeas[idxs, "theta"])
  
  # Midpoint between theta.tilde and theta.hash
  theta.to.eval <- (theta.hash + theta.tilde)/2
  
  # Stopping criterion
  stop <- FALSE
  if (abs(theta.to.eval - theta.hash) < bin.search.tol) {
    stop <- TRUE
  }
  
  # Return the results
  list(theta.to.eval = theta.to.eval, stop = stop)
}

#' @title Return the next point to evaluate when doing interpolation search
#' 
#' @description This function implements the interpolation search algorithm,
#' that starts from a given feasible point and looks in the given direction for
#' the root of the violation curve.
#' 
#' @param evaluations Matrix of evaluated test statistics and critical values.
#' @param dir Search direction.
#' @param iter.nbr Iteration number.
#' @param hp List of hyperparameters.
#' 
#' @returns The next point to evaluate.
#' 
gs.interpolation <- function(evaluations, dir, iter.nbr, hp) {
  
  # Extract the necessary hyperparameters
  int.search.tol <- hp[["bin.search.tol"]]
  max.iter <- hp[["max.iter"]]
  
  # Check if maximum number of iterations has been reached. If so, stop the
  # algorithm
  if (iter.nbr > max.iter) {
    return(list(theta.to.eval = NULL, stop = TRUE))
  }
  
  # Find current largest (smallest) feasible value, called theta.hash
  feas.idx <- which(evaluations[, "t.stat"] <= evaluations[, "crit.val"])
  theta.hash <- dir * max(dir * evaluations[feas.idx, "theta"])
  
  # If theta.hash is larger (smaller) than or equal to the upper (lower) bound,
  # stop the algorithm.
  if (((dir == 1) & (theta.hash >= hp$theta.ub)) |
      ((dir == -1) & (theta.hash <= hp$theta.lb))) {
    return(list(theta.to.eval = NULL, stop = TRUE))
  }
  
  # Matrix of all infeasible points
  infeas.idxs <- which(evaluations[, "t.stat"] > evaluations[, "crit.val"])
  evals.infeas <- matrix(nrow = 0, ncol = 3)
  colnames(evals.infeas) <- colnames(evaluations)
  for (idx in infeas.idxs) {
    evals.infeas <- rbind(evals.infeas, evaluations[idx, ])
  }
  
  # Indices of points in infeasible point set that are larger (smaller) than
  # theta.hash
  idxs <- which(dir * evals.infeas[, "theta"] > dir * theta.hash)
  
  # If no infeasible point which is larger (smaller) than theta.hash is found
  # yet, search for one.
  #
  # NOTE: when both theta.lb and theta.ub are checked in an initial stage, this
  #       code block will never be ran. Indeed, either theta.lb (theta.ub) is
  #       feasible, in which case the root finding algorithm would already 
  #       have been stopped (checked above). If it is infeasible, then an
  #       infeasible point that is smaller (larger) than theta.hash is known.
  iter.nbr <- 1
  if (length(idxs) == 0) {
    
    # Obtain largest (smallest) and second largest (smallest) theta values
    out.sort <- dir * sort(dir * evaluations[, "theta"], decreasing = TRUE)[1:2]
    theta.largest <- out.sort[1]
    theta.second <- out.sort[2]
    
    # Theta value to evaluate
    dist <- 2 * max(abs(theta.largest - theta.second), 1)
    theta.to.eval <- dir * (dir * theta.largest + dist)
    
    # Return theta
    return(list(theta.to.eval = theta.to.eval, stop = FALSE))
  }
  
  # Else, determine theta.tilde
  theta.tilde <- dir * min(dir * evals.infeas[idxs, "theta"])
  idx.theta.tilde <- which(evals.infeas[, "theta"] == theta.tilde)
  
  # Obtain the values of the violation curve at both points
  theta.tilde.viol <- evals.infeas[idx.theta.tilde, "t.stat"] - evals.infeas[idx.theta.tilde, "crit.val"]
  idx.theta.hash <- which(evaluations[, "theta"] == theta.hash)
  theta.hash.viol <- evaluations[idx.theta.hash, "t.stat"] - evaluations[idx.theta.hash, "crit.val"]
  theta.to.eval <- theta.hash - theta.hash.viol * (theta.tilde - theta.hash)/(theta.tilde.viol - theta.hash.viol)
  
  # Stopping criterion
  stop <- FALSE
  if (abs(theta.to.eval - theta.hash) < int.search.tol) {
    stop <- TRUE
  }
  
  # Return the results
  list(theta.to.eval = theta.to.eval, stop = stop)
}

#### Main functions EAM ####

#' @title Set default hyperparameters for EAM algorithm
#' 
#' @description This function returns a list with the (default) hyperparameters
#' used in the EAM algorithm
#' 
#' @param options A list of user-specified values for (some of) the
#' hyperparameters. These hyperparameters can include:
#' \begin{itemize}
#'  \item 'min.dist'/'max.dist': The minimum/maximum distance of sampled points
#'  from the current best value for the coefficient of interest.
#'  \item 'min.eval'/'max.eval': The minimum/maximum number of points evaluated
#'  in the initial feasible point search.
#'  \item 'nbr.init.sample.points': The total number of drawn points required in
#'  the initial drawing process.
#'  \item 'nbr.init.unif': The total number of uniformly drawn points in the
#'  initial set of starting values.
#'  \item 'nbr.points.per.iter.init': Number of points sampled per iteration in
#'  the initial drawing process.
#'  \item 'nbr.start.vals': Number of starting values for which to run the
#'  optimization algorithm for the expected improvement.
#'  \item 'nbr.opt.EI': Number of optimal theta values found by the optimization
#'  algorithm to return.
#'  \item 'nbr.extra': Number of extra randomly drawn points to add to the set
#'  of optimal theta values (to be supplied to the next E-step).
#'  \item 'min.improvement': Minimum amount that the current best root of the
#'  violation curve should improve by wrt. the its previous value.
#'  \item 'min.possible.improvement': Minimum amount that the next iteration
#'  should be able to improve upon the current best value of the root.
#'  \item 'EAM.min.iter': Minimum amount of EAM iterations to run.
#'  \item 'max.iter': Maximum amount of EAM iterations to run.
#' \end{itemize}
#' 
#' @returns List of hyperparameters for the EAM algotithm.

set.EAM.hyperparameters <- function(options) {
  
  # Define the list of hyperparameters
  hyperparams <- list(
    
    # Minimum and maximum distance of sampled points from current best theta 
    min.dist = 1e-4,
    max.dist = 1,
    
    # Minimum and maximum number of points to evaluate in initial feasible
    # search.
    min.eval = 10,
    max.eval = 100,
    
    # Total number of drawn points required in initial drawing process
    nbr.init.sample.points = 10,
    
    # Total number of uniformly drawn points in the initial set of starting
    # values
    nbr.init.unif = 20,
    
    # Number of points sampled per iteration in the initial drawing process
    nbr.points.per.iter.init = 4,
    
    # Number of starting values with which to run the optimization algorithm
    # for the expected improvement.
    nbr.start.vals = 80,
    
    # Number of optimal theta values found by the optimization algorithm to
    # return
    nbr.opt.EI = 1,
    
    # Number of extra randomly drawn points to add to the set of optimal
    # theta values (to be supplied to the next E-step).
    nbr.extra = 1,
    
    # [NOT USED] Distortion to be applied to theta.hash in order to obtain extra
    # evaluation points for the next E-step. In the code, 'min.improvement' is
    # used as the distortion value for theta.
    EAM_thetadistort = 0.005,
    
    # Minimum amount that this step should have improved upon the previous step
    min.improvement = 0.0001,
    
    # Minimum amount of improvement that can be the result of running the next
    # step.
    min.possible.improvement = 0.005,
    
    # Minimum number of EAM iterations that should be run
    EAM.min.iter = 4,
    
    # Maximum number of EAM iterations that should be run
    max.iter = 100
  )
  
  # Overwrite default values with user specified values if needed
  hyperparams[names(options)] <- options
  
  # Return the result
  hyperparams
} 

#' @title Main function to run the EAM algorithm
#' 
#' @description This function implements the EAM search strategy as described in
#' Kaido, Molinari and Stoye (2019). This is a generic function, requiring the
#' specification of a test function (\code{test.fun}), as well as the
#' specification of the parameter space (\code{hyperparams}).
#' 
#' @param dir The direction of the test. \code{dir = 1} corresponds to finding
#' the upper bound of the identified set, \code{dir = -1} corresponds to finding
#' the lower bound.
#' @param test.fun The test function to be inverted in order to obtain the
#' identified set. It should take a scalar parameter as argument (i.e. the
#' specified value of a component of the full parameter vector) and return a
#' list with named elements \code{list(theta, t.stat, crit.val)}, where
#' \code{theta} is the scalar value that was tested, \code{t.stat} is the value
#' of the test statistic and \code{crit.val} is the critical value to be used in
#' determining whether to reject or not reject.
#' @param hyperparams A list of hyperparameters needed in order for this method
#' to run (see \code{set.EAM.hyperparameters.R}).
#' @param evaluations Matrix of already evaluated points, of which at least one
#' is feasible. When \code{evaluations = NULL} (default), the initial feasible
#' point search will be executed first. 
#' @param time.run.duration Boolean value indicating whether to time each step
#' in the EAM algorithm. Requires \code{chronometer.R}. Default is
#' \code{time.run.duration = FALSE}.
#' @param verbose Boolean value indicating whether or not to print run time
#' updates to the console. Default is \code{verbose = FALSE}.
#'
#' @returns List containing the evaluations of the test statistic and critical
#' values, convergence information, and run times.
#'
#' @refernces Kaido et al. (2019). Confidence intervals for projections of
#' partially identified parameters. Econometrica. 87(4):1397-1432.
#' 
EAM <- function(dir, test.fun, hyperparams, evaluations = NULL,
                time.run.duration = FALSE, verbose = 0, picturose = FALSE) {
  
  # If the run times should be recorded, initialize a chronometer object
  if (time.run.duration) {
    source("chronometer.R")
    chronometer <- Chronometer$new()
    chronometer$start()
  } else {
    chronometer <- NULL
  }
  
  # Returns a set of points, of which at least one will be feasible
  if (!is.null(evaluations)) {
    fps.fail_flag <- FALSE
    if (all(evaluations[, "t.stat"] > evaluations[, "crit.val"])) {
      fps.fail_flag <- TRUE
    }
  } else {
    fps.out <- feasible_point_search(test.fun, hyperparams, verbose)
    evaluations <- fps.out[["evaluations"]]
    fps.fail_flag <- fps.out[["stop"]]
    if (time.run.duration) {
      chronometer$record.leg("finished.fps")
    }
  }
  
  # Get the tentative best value for theta
  feas.idxs <- which(evaluations[, "t.stat"] <= evaluations[, "crit.val"])
  theta.hash <- dir * max(dir * evaluations[feas.idxs, "theta"])
  
  # Initialize variables used in the main loop
  iter.nbr <- 1               # The index of the current iteration
  converged <- fps.fail_flag  # Convergence flag
  ptc.Estep <- c()            # Vector of points to check in the E-step
  opt.theta.prev <- NULL      # Optimal theta value of previous iteration
  
  # Perform main loop
  if (verbose >= 3) {
    message("Entering main loop...")
  }
  while (!converged) {
    
    #### Evaluation step ####
    
    # Run the test on the theta values returned in the previous iteration.
    # (This step is skipped during the first iteration)
    if (iter.nbr > 1) {
      if (verbose >= 3) {
        message("Doing E-step...")
      }
      E_step.out <- E_step(ptc.Estep, test.fun, dir, evaluations, verbose)
      evaluations <- E_step.out[["evaluations"]]
      theta.hash <- E_step.out[["theta.hash"]]
    }
    
    if (time.run.duration) {
      chronometer$record.leg("finished.Estep")
    }
    
    #### Approximation step ####
    
    if (verbose >= 3) {
      message("Doing A-step...")
    }
    
    # Fit a Kriging model using the points (theta, t.stat(theta) - c(theta)).
    fit.krige <- A_step(evaluations, verbose)
    
    if (time.run.duration) {
      chronometer$record.leg("finished.Astep")
    }
    
    #### Maximization step ####
    
    if (verbose >= 3) {
      message("Doing M-step...")
    }
    
    # Find promising value(s) of theta to be checked in the next iteration of
    # this algorithm. To counter greediness of the search, also include some
    # randomly selected points.
    M_step.out <- M_step(dir, evaluations, theta.hash, fit.krige, c, par.space,
                         hyperparams, verbose)
    opt.res <- M_step.out[["opt.res"]]
    extra.points <- M_step.out[["extra.points"]]
    
    if (time.run.duration) {
      chronometer$record.leg("finished.Mstep")
    }
    
    #### Check convergence of the algorithm and prepare for next iteration ####
    
    # Get the most promising theta value as given by the M-step
    mp.theta.next <- opt.res[which.max(opt.res[, "opt.EI"]), 1]
    
    # Check convergence
    converged <- EAM.converged(opt.theta.prev, evaluations, mp.theta.next,
                               iter.nbr, dir, hyperparams, verbose)
    
    # Increment iteration number
    iter.nbr <- iter.nbr + 1
    
    # Prepare for the next iteration: define vector of points to be checked in
    # the E-step of the next iteration.
    ptc.Estep <- c(opt.res[1], extra.points)
    opt.theta.prev <- theta.hash
    
  } # End main loop
  
  # Convergence of algorithm
  if (fps.fail_flag) {
    converge <- 2 # No initial feasible point was found
  } else if (iter.nbr <= hyperparams[["max.iter"]]) {
    converge <- 1 # Algorithm converged within max number of iterations
  } else {
    converge <- 0 # Algorithm didn't converge within max number of iterations
  }
  
  # Stop the chronometer
  if (time.run.duration) {
    chronometer$stop("algorithm finished")
  }
  
  # Return the results
  list(evaluations = evaluations, converge = converge, 
       chronometer = chronometer)
}

#### Main functions grid search ####

#' @title Set default hyperparameters for grid search algorithm
#' 
#' @description This function returns a list with the (default) hyperparameters
#' used in the grid search algorithm
#' 
#' @param options A list of user-specified values for (some of) the
#' hyperparameters. These hyperparameters could include:
#' \begin{itemize}
#'  \item 'min.eval'/'max.eval':
#'  \item 'next.gs.point':
#'  \item 'step.size':
#'  \item 'bin.search.tol':
#'  \item 'max.iter':
#' \end{itemize}
#' 
#' @returns List of hyperparameters for the gridsearch and binary search
#' algorithms.
#' 
set.GS.hyperparameters <- function(options) {
  
  # Define the list of hyperparameters
  hyperparams <- list(
    
    # Minimum and maximum number of points to evaluate in initial feasible
    # search. min.eval must be at least 2 in order for the binary search
    # algorithm to work properly.
    min.eval = 10,
    max.eval = 100,
    
    # Type of grid search to be carried out
    next.gs.point = gs.binary,
    
    # Step size to be used in rudimentary grid search
    step.size = 1,
    
    # Convergence tolerance to be used in binary search grid search
    bin.search.tol = 1e-3,
    
    # Maximum number of iterations to run in the grid search algorithm
    max.iter = 100
  )
  
  # Precondition check
  if (identical(hyperparams[["next.gs.point"]], gs.binary) & (hyperparams[["min.eval"]] < 2)) {
    stop("When binary search is used, min.eval should be at least 2.")
  }
  
  # Overwrite default values with user specified values if needed
  hyperparams[names(options)] <- options
  
  # Return the result
  hyperparams
}

#' @title Grid search algorithm for finding the identified set
#' 
#' @description This function implements the gridsearch and binary search
#' algorithms used to compute the roots of the violation curve and hence in 
#' estimating the identified intervals. 
#' 
#' @param dir Search direction.
#' @param test.fun The test function to be inverted in order to obtain the
#' identified set. It should take a scalar parameter as argument (i.e. the
#' specified value of a component of the full parameter vector) and return a
#' list with named elements \code{list(theta, t.stat, crit.val)}, where
#' \code{theta} is the scalar value that was tested, \code{t.stat} is the value
#' of the test statistic and \code{crit.val} is the critical value to be used in
#' determining whether to reject or not reject.
#' @param hyperparams List of hyperparameters.
#' @param evaluations Matrix of already evaluated points, of which at least one
#' is feasible. When \code{evaluations = NULL} (default), the initial feasible
#' point search will be executed first.
#' @param time.run.duration Boolean value indicating whether to time each step
#' in the EAM algorithm. Requires \code{chronometer.R}. Default is
#' \code{time.run.duration = FALSE}.
#' @param verbose Boolean value indicating whether or not to print run time
#' updates to the console. Default is \code{verbose = FALSE}.
#' 
#' @returns List containing the evaluations of the test statistic and critical
#' values, convergence information, and run times.
#' 
gridSearch <- function(dir, test.fun, hyperparams, evaluations = NULL,
                       time.run.duration = FALSE, verbose = 0,
                       picturose = FALSE) {
  
  # If the run times should be recorded, initialize a chronometer object
  if (time.run.duration) {
    source("chronometer.R")
    chronometer <- Chronometer$new()
    chronometer$start()
  } else {
    chronometer <- NULL
  }
  
  # Returns a set of points, of which at least one will be feasible
  if (!is.null(evaluations)) {
    fps.fail_flag <- FALSE
    if (all(evaluations[, "t.stat"] > evaluations[, "crit.val"])) {
      fps.fail_flag <- TRUE
    }
  } else {
    fps.out <- feasible_point_search(test.fun, hyperparams, verbose)
    evaluations <- fps.out[["evaluations"]]
    fps.fail_flag <- fps.out[["stop"]]
    if (time.run.duration) {
      chronometer$record.leg("finished.fps")
    }
  }
  
  # Initialize variables and function used in the main loop
  iter.nbr <- 0
  stop <- fps.fail_flag
  max.iter <- hyperparams[["max.iter"]]
  next.gs.point <- hyperparams[["next.gs.point"]]
  
  # Main loop
  if (verbose >= 3) {
    message("Entering main loop...")
  }
  while (!stop) {
    
    # Increment iteration number
    iter.nbr <- iter.nbr + 1
    
    # Get new r value and stopping criterion
    gnp.out <- next.gs.point(evaluations, dir, iter.nbr, hyperparams)
    r <- gnp.out[["theta.to.eval"]]
    stop <- gnp.out[["stop"]]
    
    # If a next test should be performed...
    if (!is.null(r)) {
      
      # Message to user. Update plot
      if (verbose >= 2) {
        message(sprintf("Iteration %s. Testing for r = %.3f", iter.nbr, r))
      }
      if (picturose) {
        plot.addpte(r)
      }
      
      # Perform the test
      res.Bei <- test.fun(r)
      
      # Store the results
      evaluations <- rbind(evaluations, res.Bei)
      
      # Record the run time of the iteration
      if (time.run.duration) {
        chronometer$record.leg(paste("Iter", iter.nbr))
      }
      
      # Update plot
      if (picturose) {
        plot.addpte.eval(evaluations[nrow(evaluations), , drop = FALSE])
      }
      
    } else if (!stop) {
      stop(paste("Next point is NULL while stop = FALSE. This should not",
                 "happen. Please contact the devs."))
    }
  } # End main loop
  
  # Convergence of algorithm
  if (fps.fail_flag) {
    converge <- 2 # No initial feasible point was found
  } else if (iter.nbr <= hyperparams[["max.iter"]]) {
    converge <- 1 # Algorithm converged within max number of iterations
  } else {
    converge <- 0 # Algorithm didn't converge within max number of iterations
  }
  
  # Stop the chronometer
  if (time.run.duration) {
    chronometer$stop("algorithm finished")
  }
  
  # Return the results
  list(evaluations = evaluations, chronometer = chronometer, 
       converge = converge)
}



