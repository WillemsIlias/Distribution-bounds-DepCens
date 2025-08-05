
# Load dependencies
if (!exists("simulationFunctionsSources")) {
  source("simulationFunctions.R")
}

#' @title Check argument consistency.
#' 
#' @description This function checks whether the arguments supplied to the
#' main estimation function \code{pi.surv} are valid. When arguments are
#' invalid, the an exception is thrown.
#' 
#' @inheritParams pi.surv
#' 
check.args.pisurv <- function(data, idx.param.of.interest, idxs.c, t, par.space,
                              search.method, add.options) {
  
  #### Checks for the data ####
  
  # The provided data should be given as a data frame
  if (class(data) != "data.frame") {
    stop("The provided data should be given as a data frame.")
  }
  
  # Check proper naming of columns
  colnames.unchecked <- colnames(data)
  if (!("Y" %in% colnames(data))) {
    stop(paste("The column containing the survival times should be named 'Y'",
                "(case sensitive)."))
  } else {
    colnames.unchecked <- setdiff(colnames.unchecked, "Y")
  }
  if (!("Delta" %in% colnames(data))) {
    stop(paste("The column containing the censoring indicator should be named",
                "Delta (case\n       sensitive)."))
  } else {
    colnames.unchecked <- setdiff(colnames.unchecked, "Delta")
  }
  if (!("X0" %in% colnames(data))) {
    stop(paste("The given data should contain an intercept column, named 'X0'",
                "(case sensitive)."))
  } else {
    colnames.unchecked <- setdiff(colnames.unchecked, "X0")
  }
  if (!("X1" %in% colnames(data))) {
    stop(paste("The given data should contain at least one covariate.",
                "covariates should be named\n       'X1', 'X2', ..."))
  }
  if (sum(grepl("[X][[:digit:]]+$", colnames(data))) !=
      max(as.numeric(gsub("X", "", colnames(data))[grepl("[X][[:digit:]]+$", colnames(data))])) + 1) {
    stop(paste("Invalid naming of the covariates detected. Covariates should",
                "be named 'X1', 'X2',\n       ... E.g. the case of two covariates",
               "named 'X1' and 'X3' is not allowed."))
  } else {
    colnames.unchecked <- setdiff(colnames.unchecked, colnames(data)[grepl("[X][[:digit:]]+$", colnames(data))])
  }
  if (length(colnames.unchecked) > 0) {
    stop(paste("Columns detected which do not have one of the",
               "required/allowed column names\n       (see documentation). Please",
               "remove these before excecution to ensure proper\n      ",
               "functioning of this function.",
               sprintf("\n\n       (Problematic columns: %s)", 
                       paste(colnames.unchecked, collapse = ", "))))
  }
  
  # Check valid variable types
  if (any(apply(data, 2, class) != "numeric")) {
    stop("All variables should be given as numeric.")
  }
  
  # Check valid values for censoring indicator
  if (!all(data$Delta %in% c(0, 1))) {
    stop("Censoring indicator values can only be 0 or 1.")
  }
  
  # Check correct specification of idx.param.of.interest
  n.cov <- sum(grepl("X[[:digit:]]+$", colnames(data))) - 1
  n.param <- n.cov + 1
  if (class(idx.param.of.interest) == "character") {
    if (idx.param.of.interest != "all") {
      stop("Invalid specification for idx.param.of.interest.")
    }
  } else {
    if (!(class(idx.param.of.interest)) %in% c("numeric", "integer")) {
      stop("When idx.param.of.interest is not 'all', it should be an integer")
    }
    if (idx.param.of.interest != round(idx.param.of.interest)) {
      stop("When idx.param.of.interest is not 'all', it should be an integer")
    }
    if (idx.param.of.interest > n.param | idx.param.of.interest < 1) {
      stop("Invalid index of parameter of interest")
    }
  }
  
  # Check correct specification of idxs.c
  if (!all(class(idxs.c) %in% c("numeric", "integer"))) {
    stop("idxs.c should be a vector of numerics.")
  }
  if (!all(idxs.c == round(idxs.c))) {
    stop("idxs.c should be integer valued.")
  }
  if (any(idxs.c > n.cov) | any (idxs.c < 1)) {
    stop("Invalid indices detected in idxs.c")
  }
  
  # Check valid specification of time point of interest
  if (class(t) != "numeric") {
    stop("Time point of interest t should be numeric.")
  }
  
  # Check correct specification of parameter space.
  if (!("matrix" %in% class(par.space))) {
    stop("Class of par.space should be 'matrix'.")
  }
  if (any(apply(par.space, 1:2, class) != "numeric")) {
    stop(paste("Elements of the matrix containing the bounds on the parameter",
               "space should be numeric."))
  }
  if (ncol(par.space) != 2) {
    stop("par.space must have 2 columns.")
  }
  if (any(apply(par.space, 1, function(row) {row[1] >= row[2]}))) {
    stop("Invalid bounds on parameter space detected.")
  }
  if (nrow(par.space) != n.cov + 1) {
    stop("Invalid number of rows in parameter space given the covariates",
         "provided in the data set.")
  }
  
  # Check correct specification of search.method.
  if (class(search.method) != "character") {
    stop("Invalid specification of search method.")
  }
  if (!(search.method %in% c("EAM", "GS"))) {
    stop("Invalid specification of search.method")
  }
}

#' @title Partially identify the coefficients in the model \Lambda(x^\top \beta(t))
#' for the given data set. This methodology implements the one described in
#' Willems et al. (2024+).
#' 
#' @description This function estimates bounds on the coefficients the single-
#' index model \Lambda(x^\top \beta(t)) for the conditional CDF of the event time.
#' 
#' @param data Data frame containing the data on which to fit the model. The
#' columns should be named as follows: 'Y' = observed timed, 'Delta' = censoring
#' indicators, 'X0' = intercept column, 'X1' - 'Xp' = covariates.
#' @param idx.param.of.interest Index of element in the covariate vector for
#' which the identified interval should be estimated. It can also be specified
#' as \code{idx.param.of.interest = "all"}, in which case identified intervals
#' will be computed for all elements in the parameter vector.
#' @param idxs.c Vector of indices of the continuous covariates. Suppose the
#' given data contains 5 covariates, of which 'X2' and 'X5' are continuous, this
#' argument should be specified as \code{idxs.c = c(2, 5)}.
#' @param t Time point for which to estimate the identified set of \beta(t).
#' @param par.space Matrix containing bounds on the space of the parameters. The
#' first column corresponds to lower bounds, the second to upper bounds. The i'th
#' row corresponds to the bounds on the i'th element in the parameter vector.
#' @param search.method The search method to be used to find the identified
#' interval. Default is \code{search.method = "GS"}.
#' @param add.options List of additional options to be specified to the method.
#' These options can range from 'standard' hyperparameters such as the
#' confidence level of the test and number of instrumental functions to be used,
#' to technical hyperparameters regarding the search method and test
#' implementation. For the latter, we refer to the documentations of
#' \code{set.hyperparameters}, \code{set.EAM.hyperparameters} and
#' \code{set.GS.hyperparameters}. We recommend to use the default parameters,
#' unless you really know what you are doing.
#' @param picturose Picturosity flag. If \code{TRUE}, a plot illustrating the
#' workings of the algorithm will updated during runtime. Default is
#' \code{picturose = FALSE}.
#' @param parallel Flag for whether or not parallel computing should be used.
#' Default is \code{parallel = FALSE}. When \code{parallel = TRUE}, this
#' implementation will use \code{min(detectCores() - 1, 10)} cores to construct
#' the parallel back-end.
#' 
#' @returns Matrix containing the identified intervals of the specified
#' coefficients, as well as corresponding convergence information of the
#' estimation algorithm.
#' 
#' @examples
#' \donttest{
#'   ## Example 1:
#'
#'   #     - Link function: AFT link function (default setting)
#'   #     - Number of IF: 5 IF per continuous covariate (default setting)
#'   #     - Search method: Binary search
#'   #     - Type of IF: Cubic spline functions for continuous covariate, indicator
#'   #       function for discrete covariate (default setting).
#'   
#'   # Load 'survival' package in R.
#'   library("survival")
#'   
#'   # Load and preprocess data
#'   data <- survival::lung
#'   data[, "intercept"] <- rep(1, nrow(data))
#'   data[, "status"] <- data[, "status"] - 1
#'   data <- data[, c("time", "status", "intercept", "age", "sex")]
#'   colnames(data) <- c("Y", "Delta", "X0", "X1", "X2")
#'   
#'   # Settings for main estimation function
#'   idx.param.of.interest <- 1 # Interest in effect of age
#'   idxs.c <- 1                # X1 (age) is continuous
#'   t <- 200                   # Model imposed at t = 200
#'   search.method <- "GS"      # Use binary search
#'   par.space <- matrix(rep(c(-10, 10), 3), nrow = 3, byrow = TRUE)
#'   
#'   # Estimate the identified intervals
#'   pi.surv(data, idx.param.of.interest, idxs.c, t, par.space, search.method)
#'   
#'   # Example 2:
#'   
#'   #     - Link function: Cox link function
#'   #     - Number of IF: 7 IF per continuous covariate
#'   #     - Search method: EAM algorithm
#'   #     - Type of IF: Cubic spline functions for continuous covariate, indicator
#'   #       function for discrete covariate (default setting).
#'   #     - Confidence level: 1 - 0.05/3 (e.g. for Bonferroni correction)
#'   
#'   # Load 'survival' package in R.
#'   library("survival")
#'   
#'   # Load and preprocess data
#'   data <- survival::lung
#'   data[, "intercept"] <- rep(1, nrow(data))
#'   data[, "status"] <- data[, "status"] - 1
#'   data <- data[, c("time", "status", "intercept", "age", "sex")]
#'   colnames(data) <- c("Y", "Delta", "X0", "X1", "X2")
#'   
#'   # Settings for main estimation function
#'   idx.param.of.interest <- 1                  # Interest in effect of age
#'   idxs.c <- 1                                 # X1 (age) is continuous
#'   t <- 200                                    # Model imposed at t = 200
#'   search.method <- "EAM"                      # Use EAM algorithm
#'   add.options <- list(alpha = 1 - 0.05/3,     # Specify confidence level of test 
#'                       n.if.per.cov = 5)       # Specify nbr. of inst. func.
#'   par.space <- matrix(rep(c(-10, 10), 3), nrow = 3, byrow = TRUE)
#'   
#'   # Estimate the identified intervals
#'   pi.surv(data, idx.param.of.interest, idxs.c, t, par.space, search.method,
#'           add.options)
#' }
#' 
#'
#' @references Willems, I., Beyhum, J. and Van Keilegom, I. (2024+). Partial
#' identification for a class of survival models under dependent censoring.
#' (In preparation).
#' 
pi.surv <- function(data, idx.param.of.interest, idxs.c, t, par.space,
                    search.method = "GS", add.options = list(),
                    picturose = FALSE, parallel = FALSE) {
  
  #### Consistency checks ####
  
  check.args.pisurv(data, idx.param.of.interest, idxs.c, t, par.space,
                    search.method, add.options)
  
  # If required, set-up parallel back-end. The variable 'clust' is defined
  # globally to facilitate usage of parallel cluster throughout the
  # implementation
  if (parallel) {
    n.cores <- min(detectCores() - 1, 10)
    clust <<- makeCluster(n.cores)
    registerDoParallel(clust)
  }
  
  #### Set the hyperparameters of the algorithm ####
  
  # Number of covariates
  n.cov <- sum(grepl("X[1-9][[:digit:]]*", colnames(data)))
  
  # Transform argument 'idx.param.of.interest' to list of unit vectors with a
  # 1 placed at the respective index.
  c.to.check <- list()
  if (class(idx.param.of.interest) == "character") {
    for (i in 1:(n.cov + 1)) {
      c.vec <- rep(0, n.cov + 1)
      c.vec[i] <- 1
      c.to.check <- c(c.to.check, list(c.vec))
    }
  } else {
    c.vec <- rep(0, n.cov + 1)
    c.vec[idx.param.of.interest] <- 1
    c.to.check <- list(c.vec)
  }
  
  # Pre-set algorithm hyperparameters
  options <- list(n.if.per.cov = 5,
                  K.bar = 3,
                  B = 600,
                  next.gs.point = gs.binary,
                  alpha = 0.95,
                  link.function = "AFT_ll",
                  inst.func.family = "cd",
                  degree = 3,
                  G.c = G.spline,
                  idxs.c = idxs.c)
  
  # Overwrite pre-sets with user-specified options
  options[names(add.options)] <- add.options
  
  #### Find the identified interval(s) #### 
  
  # Initialize object that will store the results
  results <- matrix(nrow = length(c.to.check), ncol = 4)
  rownames(results) <- paste0("X", matrix(unlist(c.to.check), nrow = length(c.to.check)) %*% 0:n.cov)
  colnames(results) <- c("lower", "upper", "conv.l", "conv.u")
  
  # Run the algorithm for the selected elements of the parameter vector
  for (c in c.to.check) {
    
    # Get the identified set of the covariate of interest
    fis.out <- find.identified.set(c, t, par.space, data, search.method, options,
                                   verbose = 3, picturose = picturose,
                                   parallel = parallel)
    
    # Store the results
    row.name <- paste0("X", which(c == 1) - 1)
    results[row.name, c("lower", "upper")] <- fis.out$ident.set
    results[row.name, c("conv.l", "conv.u")] <- c(fis.out$converge1, fis.out$converge2)
  }
  
  # Remove parallel back-end
  if (parallel) {
    stopCluster(clust)
  }
  
  # Return the results
  results
}

#' @title Subset a given data set based on the provided criteria.
#' 
#' @description This function subsets a data set based on the given criteria.
#' 
#' @param data Data frame (see documentation of 'pi.surv.R')
#' @param criteria A data frame wherein each row contains the variable for which
#' the criterion should apply, a comparing operator, and the reference value.
#' @param keep.covs Names of covariates to keep. Default \code{keep.covs = NULL}.
#' @param max.row Maximum rows in the resulting data set. Default is
#' \code{max.row = NULL}
#' 
#' @noRd
#' 
subset.data <- function(data, criteria, keep.covs = NULL, max.row = NULL) {
  data.sub <- data
  for (row.idx in 1:nrow(criteria)) {
    
    # Get variable name or index for which to apply the criterion
    var.name <- criteria[row.idx, 1]
    
    # Apply criterion over data set
    expr <- sprintf("data.sub[which(data.sub[, var.name] %s criteria[row.idx, 3]),]",
                    criteria[row.idx, 2])
    data.sub <- eval(parse(text = expr))
  }
  
  # Retain only the specified covariates
  name.dict <- NULL
  if (!is.null(keep.covs)) {
    data.sub <- data.sub[, c("Y", "Delta", "X0", keep.covs)]
    name.dict <- data.frame(
      old = c("Y", "Delta", "X0", keep.covs),
      new = c("Y", "Delta", paste0("X", 0:length(keep.covs)))
    )
    colnames(data.sub) <- name.dict[match(colnames(data.sub), name.dict[, 1]), 2]
  }
  
  # Retain only the first 'max.row' rows of the data frame
  if (!is.null(max.row)) {
    data.sub <- data.sub[1:min(nrow(data.sub), max.row), ]
  }
  
  # Return subsetted data and additional information
  list(data.sub, name.dict)
}

#' @title Combine bounds based on majority vote.
#' 
#' @description This function combines a list of individual identified intervals
#' to a single one based on majority vote. Note that the intersection of all
#' intervals can be viewed as a majority vote as well, so that it is included as
#' a special case.
#' 
#' @param results.list List object containing the individual identified
#' intervals.
#' @param threshold Threshold proportion of identified intervals a given value
#' should be contained in in order for it to be included in the combined
#' identified interval. For intersection bounds, set this value to \code{1}.
#' 
#' @returns The combined identified interval.
#' 
cbMV <- function(results.list, threshold) {
  
  # Obtain some information about the results
  cov.names <- rownames(results.list[[1]])
  
  # Obtain vector of lower and upper bounds
  lbs <- matrix(ncol = length(cov.names), nrow = length(results.list))
  ubs <- matrix(ncol = length(cov.names), nrow = length(results.list))
  missspec <- FALSE
  for (results.idx in 1:length(results.list)) {
    
    # Get results of this iteration
    results <- results.list[[results.idx]]
    
    # Extract the bounds on the parameters
    for (cov.name.idx in 1:length(cov.names)) {
      cov.name <- cov.names[cov.name.idx]
      lbs[results.idx, cov.name.idx] <- results[cov.name, "lower"]
      ubs[results.idx, cov.name.idx] <- results[cov.name, "upper"]
    }
    
    # Check if the model was determined to be misspecified.
    if (any(results[, c("conv.l", "conv.u")] == 2)) {
      missspec <- TRUE
    }
  }
  colnames(lbs) <- colnames(ubs) <- cov.names
  
  # Initialize object that will store the combined bounds
  bounds.MV <- matrix(rep(c(-Inf, Inf, 1, 1), length(cov.names)),
                          nrow = length(cov.names),
                          byrow = TRUE)
  colnames(bounds.MV) <- colnames(results.list[[1]])
  rownames(bounds.MV) <- cov.names
  
  # If model was misspecified at any of the tested points, also the combined
  # model is misspecified.
  if (missspec) {
    bounds.MV[, 3:4] <- 2
    return(bounds.MV)
  }
  
  # Else, determine the combined bounds for each parameter of interest.
  for (cov.name in cov.names) {
    
    # Get bounds of each part
    ths <- sort(unique(c(lbs[, cov.name], ubs[, cov.name])))
    parts <- matrix(rep(ths, c(1, rep(2, length(ths) - 2), 1)), ncol = 2, byrow = TRUE)
    
    # For each part, obtain the number of votes
    parts <- cbind(parts, rep(0, nrow(parts)))
    for (part.idx in 1:nrow(parts)) {
      parts[part.idx, 3] <- sum((lbs[, cov.name] < mean(parts[part.idx, 1:2])) &
                                  (ubs[, cov.name] > mean(parts[part.idx, 1:2])))
    }
    
    # Subset to parts getting majority vote
    MV.parts <- parts[parts[, 3]/length(lbs[, cov.name]) >= threshold, 1:2, drop = FALSE]
    
    # Recombine neighbouring bounds
    row.idx <- 1
    while (row.idx < nrow(MV.parts)) {
      if (MV.parts[row.idx, 2] == MV.parts[row.idx + 1, 1]) {
        MV.parts[row.idx, 2] <- MV.parts[row.idx + 1, 2]
        MV.parts <- MV.parts[setdiff(1:nrow(MV.parts), row.idx + 1), , drop = FALSE]
      } else {
        row.idx <- row.idx + 1
      }
    }
    
    # If the result is not an interval, take lower and upper bounds, but warn
    # the user.
    if (nrow(MV.parts) > 1) {
      MV.parts <- matrix(c(min(MV.parts), max(MV.parts)), nrow = 1)
      warning("Non-connected identified set. Returning outer identified interval.")
    }
    
    # Store the result
    bounds.MV[cov.name, ] <- cbind(MV.parts, 1, 1)
  }
  
  # Return the result
  bounds.MV
}

#' @title Get the file name of the results from Pancreas data applications
#' 
#' @description This function obtains the names of the files containing the
#' results of the Pancreas data application.
#' 
#' @param args The vector of arguments supplied to the data application wrapper
#' function. 
#' 
#' @noRd
#' 
get.file.name.pancreas <- function(args, idx.param.of.interest, master.dir) {
  
  # Create directory if necessary
  check_create.dir(master.dir)
  
  # Get name of variable
  if (idx.param.of.interest == 1) {
    var.name <- "intercept"
  } else if (idx.param.of.interest == 2) {
    var.name <- "age"
  } else if (idx.param.of.interest == 3) {
    var.name <- "tumorsize"
  }
  
  # Create file name
  file.name <- sprintf("crit-%d__nifpc-%d_var-%s.csv", args["criteria.idx"],
                       args["n.if.per.cov"], var.name)
  
  # Return path name
  paste0(c(master.dir, file.name), collapse = "/")
}

#### Plotting functions (picturose) ####

#' @title Clear plotting window
#' 
#' @description This function clears the plotting window
#' 
clear.plt.wdw <- function() {
  tryCatch(invisible(dev.off()), error = function(e) {invisible(e)})
}

#' @title Draw base plot
#' 
#' @description This functon draws the base plot, used when
#' \code{picturose = TRUE}.
#' 
#' @param c Projection vector
#' @param hp List of hyperparameters
#' 
plot.base <- function(c, hp) {
  
  # Clear plotting window
  clear.plt.wdw()
  
  # Make base plot
  c.idx <- which(c == 1) - 1
  plot(x = c(hp$theta.lb, hp$theta.ub), y = c(0, 0), type = 'l',
       xlab = bquote("parameter"~"space"~"for"~theta[.(c.idx)]),
       ylab = "",
       yaxt = 'n')
}

#' @title Draw points to be evaluated
#' 
#' @description This function draws the points to be evaluated.
#' 
#' @param pte Vector of points to be evaluated.
#' @param col Color of the points.
#' 
plot.addpte <- function(pte, col = "orange") {
  points(x = pte, y = rep(0, length(pte)), pch = 16, col = col)
  points(x = pte, y = rep(0, length(pte)), pch = 1, col = "black")
}

#' @title Draw evaluated points.
#' 
#' @description This function draws evaluated points. Feasible points are
#' indicated in green, red points correspond to infeasible points.
#' 
#' @param evaluations Matrix of evaluations to be drawn.
#' 
plot.addpte.eval <- function(evaluations) {
  feas <- evaluations[, 2] <= evaluations[, 3]
  col <- ifelse(feas, "green", "red")
  points(x = evaluations[, 1], y = rep(0, nrow(evaluations)), col = col, pch = 16)
  points(x = evaluations[, 1], y = rep(0, nrow(evaluations)), pch = 1, col = "black")
}




