
# Load dependencies of this script
if ((sys.nframe() > 0) & exists("main.func.dir.name")) {
  source(paste0(main.func.dir.name, "/searchStrategies.R"))
  source(paste0(main.func.dir.name, "/EstimationAlgorithmBei.R"))
  source(paste0(main.func.dir.name, "/lowLevelFunctions.R"))
  source(paste0(main.func.dir.name, "/chronometer.R"))
} else {
  source("searchStrategies.R")
  source("EstimationAlgorithmBei.R")
  source("lowLevelFunctions.R")
  source("chronometer.R")
  require(doParallel)
}

# Set global flag that the functions in this file have been sourced
simulationFunctionsSources <<- TRUE


#' @title Define the hyperparameters used for finding the identified interval
#' 
#' @description This function defines all the necessary hyperparameters used to
#' run the methodology.
#' 
#' @param data Data frame.
#' @param par.space Bounds on the parameter space.
#' @param c Projection vector.
#' @param search.method Search method to use (\code{"EAM"} or \code{"GS"})
#' @param options List of user specified hyperparameters that will substitute
#' the corresponding default values. This list can contain the entries:
#' \begin{itemize}
#'  \item 'cov.ranges': known bounds on each of the covariates in the data set.
#'  \item 'norm.func.name': Name of the normalization function to be used. Can
#'  be either "normalize.covariates1" or "normalize.covariates2" (default).
#'  The former is a simple elementwise rescaling. The latter uses the PCA
#'  approach as discussed in Willems et al. (2024+).
#'  \item 'inst.func.family': Family of instrumental functions to be used for
#'  all covariates. Options are "box", "spline" and "cd". The former two are
#'  only applicable for continuous covariates. The latter can also handle
#'  discrete covariates. Default is "cd".
#'  \item 'G.c': The class of instrumental functions used for the continuous
#'  covariates in the model, in case "cd" is selected as
#'  \code{inst.func.family}. Options are "box" and "spline". Default is "spline".
#'  \item 'degree': The degree of the B-spline functions, should they be used as
#'  instrumental functions for the continuous covariates. Default is 3.
#'  \item 'link.function': Name of the link function to be used. Options are
#'  "AFT_ll" for the AFT model with log-logistic baseline, or "Cox_wb" for the
#'  Cox PH model (originally with Weibull baseline, but now for a general)
#'  baseline hazard).
#'  \item 'K.bar': Number of refinement steps when obtaining the critical value.
#'  See Bei (2024).
#'  \item 'B': Number of bootstrap samples to be used when obtaining the
#'  bootstrap distribution of the test statistic.
#'  \item 'ignore.empty.IF': Boolean value indicating whether instrumental
#'  functions with empty support should be ignored (cf. Willems et al., 2024).
#'  Default is FALSE. The feature \code{ignore.empty.IF = TRUE} is experimental,
#'  so there might exist edge cases for which the implementation will fail to
#'  run.
#' \end{itemize}
#' Other (hidden) options can also be overwritten, though we highly discourage
#' this. If necessary, you can consult the source code of this functions to
#' find the names of the desired parameters and add their name alongside their
#' desired value as an entry in \code{options} (e.g.
#' \code{options$min.var <- 1e-4}. Again, not recommended!).
#' 
#' @returns The list of hyperparameters.
#' 
set.hyperparameters <- function(data, par.space, c, search.method, options) {
  
  #### General tuning hyperparameters ####
  
  # Number of covariates (excluding the intercept)
  n.cov <- sum(grepl("X[[:digit:]]+", colnames(data))) - 1
  n <- nrow(data)
  
  # Type of covariates. DGP option is used when running the simulations of
  # Willems et al. (2024+).
  cov.idxs <- 1:n.cov
  if (is.null(options[["idxs.c"]])) {
    if (options[["DGP"]] <= 20) {
      idxs.c <- cov.idxs
    } else if (20 < options[["DGP"]] & options[["DGP"]] <= 40) {
      idxs.c <- cov.idxs[1:ceiling(n.cov/2)]
    } else if (40 < options[["DGP"]] & options[["DGP"]] <= 60) {
      idxs.c <- integer(0)
    }
  } else {
    idxs.c <- options[["idxs.c"]]
  }
  if (!all(idxs.c %in% 1:n.cov)) {
    stop("Invalid value for 'idxs.c'")
  }
  
  # Tuning parameters
  kappa.n <- sqrt(log(n))
  lambda.n <- sqrt(n) * kappa.n^2
  epsilon.n <- sqrt(log(kappa.n^2)/n)
  delta.n <- min(1/n, 10^(-4))
  
  # Artificial minimum variance (ensure non-nullity of variance).
  min.var <- 1e-6
  
  # Instrumental function hyperparameters
  n.if.per.cov <- options[["n.if.per.cov"]]
  if (is.null(n.if.per.cov)) {
    stop("Number of instumental functions per covariate should be specified")
  }
  inst.func.family <- options[["inst.func.family"]]
  if (is.null(inst.func.family)) {
    inst.func.family <- "cd"
  }
  ignore.empty.IF <- FALSE
  if ("ignore.empty.IF" %in% names(options)) {
    if (class(options[["ignore.empty.IF"]]) == "logical") {
      ignore.empty.IF <- options[["ignore.empty.IF"]]
    } else {
      stop("Specified value for 'ignore.empty.IF' should be TRUE/FALSE.")
    }
  }
  cov.ranges <- options[["cov.ranges"]]
  
  # Normalization function (used in instrumental function)
  norm.func.name <- options[["norm.func.name"]]
  if (is.null(norm.func.name)) {
    norm.func.name <- "normalize.covariates2"
  }
  if (norm.func.name == "normalize.covariates1") {
    norm.func <- normalize.covariates
  } else if (norm.func.name == "normalize.covariates2") {
    norm.func <- normalize.covariates2
  }
  
  # Precompute normalized covariates. When inst.func.family == cd, this should
  # only be done on the continuous covariates.
  if (inst.func.family == "cd") {
    cols.to.include <- which(!(colnames(data) %in% paste0("X", setdiff(1:n.cov, idxs.c))))
    data.c <- data[, cols.to.include]
    cov.ranges.c <- cov.ranges[, cols.to.include]
    norm.cov.out <- norm.func(data = data.c, cov.ranges = cov.ranges.c,
                              idxs.c = "all")
  } else {
    norm.cov.out <- norm.func(data = data, cov.ranges = cov.ranges)
  }
  
  # Number of instrumental functions (computed differently depending on whether
  # or not G.cd or G.cd.mc is selected).
  if (inst.func.family == "cd") {
    
    # Number of instrumental functions pertaining to continuous covariates
    names.cov.d <- setdiff(paste("X", setdiff(1:n.cov, idxs.c), sep = ""), "X")
    n.inst.func.c <- n.if.per.cov^(n.cov - length(names.cov.d))
    
    # Number of instrumental functions pertaining to discrete covariates, taking
    # into account that some combinations of discrete covariate levels might
    # be empty.
    covariates.d <- data[, names.cov.d, drop = FALSE]
    n.inst.func.d <- max(nrow(unique(covariates.d)), 1)
    
    # Total number of covariates
    n.inst.func <- n.inst.func.c * n.inst.func.d
    
  } else if (inst.func.family == "cd.manycov") {
    
    # Precompute some necessary information
    cov.names <- colnames(data)[grep("X[1-9][[:digit:]]*$", colnames(data))]
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
    
    # Get number of instrumental functions
    n.inst.func <- max(info.manycov$cumsum)
    
  } else {
    n.inst.func <- n.if.per.cov^n.cov
  }
  
  # If "G.cd" is selected, obtain all levels of the 'combined' discrete
  # covariates present in the data.
  discrete.covariate.levels <- NULL
  if (inst.func.family == "cd") {
    names.cov.d <- setdiff(paste("X", setdiff(1:n.cov, idxs.c), sep = ""), "X")
    discrete.covariate.levels <- unique(data[, names.cov.d, drop = FALSE])
  }
  
  # Instrumental function
  if (n.cov > 3 & inst.func.family != "cd.manycov") {
    warning(paste("It is recommended to specify 'inst.func.family = 'cd.manycov'",
                  "when the number of covariates exceeds 3."))
  }
  if (inst.func.family == "box") {
    degree <- 0
    G <- function(x, g.idx) {G.box(x, g.idx, data, n.if.per.cov,
                                   cov.ranges = cov.ranges,
                                   norm.func = norm.func,
                                   norm.cov.out = norm.cov.out)}
    
  } else if (inst.func.family == "spline") {
    degree <- options[["degree"]]
    if (is.null(degree)) {degree <- 3}
    G <- function(x, g.idx) {G.spline(x, g.idx, data, n.if.per.cov,
                                      degree = degree, cov.ranges = cov.ranges,
                                      norm.func = norm.func,
                                      norm.cov.out = norm.cov.out)}
    
  } else if (inst.func.family == "cd") {
    G.c <- options[["G.c"]]
    if (is.null(G.c)) {stop("G.c should be specified when G.cd is used")}
    degree <- options[["degree"]]
    if (is.null(degree)) {degree <- 3}
    G <- function(x, g.idx) {G.cd(x = x, g.idx = g.idx, data = data,
                                  n.if.per.cov = n.if.per.cov, idxs.c = idxs.c,
                                  G.c = G.c, norm.func = norm.func,
                                  discrete.covariate.levels = discrete.covariate.levels,
                                  cov.ranges = cov.ranges,
                                  norm.cov.out = norm.cov.out, degree = degree)}
  
  } else if (inst.func.family == "cd.manycov") {
    G.c <- options[["G.c"]]
    if (is.null(G.c)) {stop("G.c should be specified when G.cd is used")}
    degree <- options[["degree"]]
    if (is.null(degree)) {degree <- 3}
    G <- function(x, g.idx) {
      G.cd.mc(x = x, g.idx = g.idx, data = data, n.if.per.cov = n.if.per.cov,
              idxs.c = idxs.c, G.c = G.c, norm.func = norm.func,
              info.manycov = info.manycov, cov.ranges = cov.ranges,
              degree = degree)
      }
    
  } else {
    stop("Unknown instrumental function family specified.")
  }
  
  # Link function to be used
  link.function <- "AFT_ll"
  if ("link.function" %in% names(options)) {
    link.function <- options[["link.function"]]
  }
  if (link.function == "AFT_ll") {
    Lambda <- Lambda_AFT_ll
    dLambda <- dLambda_AFT_ll
    inv.Lambda <- Lambda_inverse_AFT_ll
  } else if (link.function == "Cox_wb") {
    Lambda <- Lambda_Cox_wb
    dLambda <- dLambda_Cox_wb
    inv.Lambda <- Lambda_inverse_Cox_wb
  }
  
  # Number of initial values in step 4 of the algorithm of Bei (2024)
  K.bar <- options[["K.bar"]]
  if (is.null(K.bar)) {
    stop(paste0("Number of initial values in step 4 of the algorithm of Bei ", 
                "(2024) should be specified."))
  } 
  
  # Number of bootstrap samples
  B <- options[["B"]]
  if (is.null(B)) {
    stop("Number of bootstrap samples should be specified")
  } 
  
  hp <- list(
    
    # Link function
    "Lambda" = Lambda,
    "dLambda" = dLambda,
    "inv.Lambda" = inv.Lambda,
    
    # Instrumental functions
    "n.if.per.cov" = n.if.per.cov,
    "G" = G,
    "n.inst.func" = n.inst.func,
    "discrete.covariate.levels" = discrete.covariate.levels,
    "norm.func.name" = norm.func.name,
    
    # Tuning parameters 
    "kappa.n" = kappa.n,
    "lambda.n" = lambda.n,
    "epsilon.n" = epsilon.n,
    "delta.n" = delta.n,
    "min.var" = min.var,
    
    # Bounds of the parameter space for theta
    "theta.lb" = par.space[which(c == 1), 1],
    "theta.ub" = par.space[which(c == 1), 2],
    
    # Hyperparameters for computing the test statistic and critical value
    "K.bar" = K.bar,    # Number of initial values in step 4 of algorithm
                        # described in Bei (2024).
    "B" = B             # Number of bootstrap samples
  )
  
  # Overwrite default values with user specified values if needed
  hp[names(options)] <- options
  
  #### Search strategy-specific hyperparameters ####
  
  if (search.method == "EAM") {
    ss.hp <- set.EAM.hyperparameters(options)
  } else if (search.method == "GS") {
    ss.hp <- set.GS.hyperparameters(options)
  }
  
  #### Return all hyperparamters ####
  
  hp[names(ss.hp)] <- ss.hp
  return(hp)
}

#' @title Estimate the identified set of the parameters in the model.
#' 
#' @description This function estimates the identified set of the parameters in
#' the model. It does so by running the selected search algorithm twice (once
#' for minimization of the nonlinear program and once for maximization). Not to
#' be confused with the low level function 'get.identified.set.R'.
#' 
#' @param c Projection vector
#' @param t Time point of interest.
#' @param par.space Bounds on the parameter space.
#' @param data Data frame.
#' @param search.method String value indicating the search method to be used
#' for finding the identified set. Can be 'EAM' or 'GS' (grid search).
#' @param options List of user specified hyperparameters that will substitute
#' the corresponding default values.
#' @param verbose Verbosity parameter. Higher values indicate larger degrees of
#' verbosity. Default is \code{verbose = 0} (no verbosity).
#' @param time.run.duration Boolean value indicating whether run durations
#' should be timed.
#' @param picturose Picturosity flag. If \code{TRUE}, a plot illustrating the
#' workings of the algorithm will updated during runtime. Default is
#' \code{picturose = FALSE}.
#' @param parallel Flag for whether or not parallel computing should be used.
#' Default is \code{parallel = FALSE}.
#' 
#' @returns List containing the estimated identified set, convergence
#' information and run times.
#'
#' @noRd
#' 
find.identified.set <- function(c, t, par.space, data, search.method, options,
                                verbose = 0, time.run.duration = FALSE,
                                picturose = FALSE, parallel = FALSE) {
  
  #### Precondition checks ####
  
  if (verbose >= 2) {
    message("  Performing precondition checks...")
  }
  
  # Check if the projection vector is valid
  if (!all(c %in% c(0, 1)) | sum(c^2) != 1) {
    stop("Invalid argument for c")
  }
  
  # Check if the bounds of the parameter space are valid
  if (!is.matrix(par.space)) {
    stop("par.space should be a matrix")
  } else {
    if (ncol(par.space) != 2) {
      stop("par.space should contain 2 columns")
    }
  }
  
  # Check if the data is specified as required
  if (!("X0" %in% colnames(data))) {
    stop("Data frame should contain intercept column named 'X0'.")
  }
  if (!any(grepl("X[1-9]+", colnames(data)))) {
    warning("No covariates detected in the data")
  }
  
  #### Initialize hyperparameters ####
  
  if (verbose >= 2) {
    message("  Initializing hyperparameters...")
  }
  
  # Hyperparameters
  hp <- set.hyperparameters(data, par.space, c, search.method, options)
  
  # Set search algorithm
  if (search.method == "EAM") {
    search.algorithm <- EAM
  } else if (search.method == "GS") {
    search.algorithm <- gridSearch
  }
  
  # Precompute instrumental function evaluations
  inst.func.evals <- t(get.instrumental.function.evals(data, hp))
  
  # Find possible instrumental functions with empty support. If indicated,
  # remove these from the analysis.
  idx.empty.IF <- which(rowSums(inst.func.evals) == 0)
  if (!is.null(hp$ignore.empty.IF)) {
    if (hp$ignore.empty.IF & length(idx.empty.IF) > 0) {
      
      # Remove empty IF
      inst.func.evals <- inst.func.evals[setdiff(1:nrow(inst.func.evals), idx.empty.IF), ]
      
      # Update the hyperparameter information
      hp$n.inst.func <- nrow(inst.func.evals)
      
      # Warn the user
      warning(sprintf("%d empty instrumental functions removed!", length(idx.empty.IF)))
    }
  }
  
  # Test whether each instrumental function contains at least one observation.
  # If not, return with appropriate failure flag
  if (any(rowSums(inst.func.evals) == 0)) {
    return(list(ident.set = c(-Inf, Inf),
                converge1 = 3,
                converge2 = 3,
                chronometer1 = Chronometer$new(),
                chronometer2 = Chronometer$new()))
  }
  
  # Set test function
  test.fun <- function(theta) {
    test.point_Bei(theta, c, t, par.space, data, hp, verbose = FALSE,
                   inst.func.evals = inst.func.evals, alpha = options$alpha)
  }
  
  #### Find identified set ####
  
  # If required, update user and initialize plot
  if (verbose >= 2) {
    message("  Starting search for initial feasible points...")
  }
  if (picturose) {
    plot.base(c, hp)
  }
  
  # Pre-search for feasible points
  fps.out <- feasible_point_search(test.fun, hp, verbose, picturose = picturose,
                                   parallel = parallel)
  evaluations <- fps.out[["evaluations"]]
  
  # Re-set test function, this time using parallel computing, if necessary.
  test.fun <- function(theta) {
    test.point_Bei(theta, c, t, par.space, data, hp, verbose = FALSE,
                   inst.func.evals = inst.func.evals, alpha = options$alpha,
                   parallel = parallel)
  }
  
  # Run search algorithm in dir = 1
  if (verbose >= 2) {
    message("  Starting search in dir = 1")
  }
  dir <- 1
  sa.out <- search.algorithm(dir, test.fun, hp, evaluations, time.run.duration,
                             verbose, picturose = picturose)
  evaluations1 <- sa.out[["evaluations"]]
  converge1 <- sa.out[["converge"]]
  chronometer1 <- sa.out[["chronometer"]]
  
  # Run search algorithm in dir = -1
  if (verbose >= 2) {
    message("  Starting search in dir = -1")
  }
  dir <- -1
  sa.out <- search.algorithm(dir, test.fun, hp, evaluations, time.run.duration,
                             verbose, picturose = picturose)
  evaluations2 <- sa.out[["evaluations"]]
  converge2 <- sa.out[["converge"]]
  chronometer2 <- sa.out[["chronometer"]]
  
  # Combine results and find identified interval
  evaluations <- rbind(evaluations1, evaluations2)
  ident.set <- get.identified.set(evaluations)
  
  #### Return the results ####
  
  list(ident.set = ident.set, converge1 = converge1, converge2 = converge2,
       chronometer1 = chronometer1, chronometer2 = chronometer2)
}

#' @title Obtain the file name of a simulation
#' 
#' @description This function returns the filename, prefixed with relevant
#' folder paths, for the results of a simulation iteration. (This function is
#' used in simulate1D.R)
#' 
#' @param comb Either a matrix or list containing the information of the
#' simulation design.
#' @param sim.iter.nbr Number of the simulation ran (value in 1:n.sim)
#' @param seed The random seed for which the simulation is ran.
#' @param master.dir Master directory.
#' @param shortened Boolean value indicating whether a shortened version of the
#' file name should be returned.
#' @param defaults Default values of parameters that will be omitted from file
#' name if \code{shortened = TRUE}.
#' 
#' @note This function creates the relevant directories if necessary by calling
#' the function 'lowLevelFunctions::check_create.dir.R'.
#' 
#' @noRd
#' 
get.file.name <- function(comb, sim.iter.nbr, seed, master.dir,
                          shortened = FALSE, defaults = NULL) {
  
  if (class(comb)[1] == "matrix") {
    
    # Shorten inst.func.family name, if applicable.
    if (comb[, "inst.func.family"] == "cd.manycov") {
      comb[, "inst.func.family"] <- "cdmc"
    }
    
    # Directory: (link function,) search method, inst.func.family, K.bar
    sub.dir <- sprintf("Search-%s__IF-%s__Kbar-%s__Alpha-%s", comb[, "search.method"],
                       comb[, "inst.func.family"], comb[, "K.bar"], comb[, "alpha"])
    if ("link.function" %in% colnames(comb)) {
      sub.dir <- sprintf("LF-%s__%s", comb[, "link.function"], sub.dir)
    }
    if ("parametric" %in% colnames(comb)) { # Option used in misspecification simulations
      sub.dir <- sprintf("%s__Param-%s", sub.dir, comb[, "parametric"])
    }
    
    # If necessary, remove defaults
    if (shortened) {
      components <- str_split(sub.dir, "__")[[1]]
      components.shortened <- NULL
      for (component in components) {
        if (!(str_split(component, "-")[[1]][1] %in% names(defaults))) {
          components.shortened <- c(components.shortened, component)
        }
      }
      
      sub.dir <- paste(components.shortened, collapse = "__")
    }
    
    # Create directory if necessary
    check_create.dir(master.dir)
    check_create.dir(sub.dir, master.dir)
    
    # File name: (gs.method,) (n.if.per.cov,) n.sim, n, n.cov, DGP, B, t.eval
    file.name <- if (comb[, "search.method"] == "GS") {sprintf("method-%s", comb[, "gs.method"])} else NULL
    file.name <- if (comb[, "inst.func.family"] == "box") {
      paste(c(file.name, sprintf("nbpc-%s", comb[, "n.if.per.cov"])), collapse = "__")
    } else if (comb[, "inst.func.family"] == "spline") {
      paste(c(file.name, sprintf("nifpc-%s__degree-%s", comb[, "n.if.per.cov"], comb[, "degree"])), collapse = "__")
    } else if (comb[, "inst.func.family"] %in% c("cd", "cdmc")) {
      paste(c(file.name, sprintf("nifpc-%s__Gc-%s__degree-%s", comb[, "n.if.per.cov"], comb[, "inst.func.family.c"], comb[, "degree"])), collapse = "__")
    } else file.name
    gen.info <- sprintf("geninfo-%s", paste(c(comb[, "n"], comb[, "n.cov"], 
                                              comb[, "DGP"], comb[, "B"],
                                              comb[, "t.eval"], seed, sim.iter.nbr),
                                            collapse = "__"))
    file.name <- ifelse(is.null(file.name),
                        gen.info,
                        paste(c(file.name, gen.info), collapse = "__"))
    file.name <- paste0(file.name, ".csv")
    
    # If necessary, remove defaults
    if (shortened) {
      components <- str_split(file.name, "__")[[1]]
      components.shortened <- NULL
      for (component in components) {
        if (!(str_split(component, "-")[[1]][1] %in% names(defaults))) {
          components.shortened <- c(components.shortened, component)
        }
      }
      
      file.name <- paste(components.shortened, collapse = "__")
    }
    
  } else if (class(comb) == "list") {
    
    # Shorten inst.func.family name, if applicable.
    if (comb[["inst.func.family"]] == "cd.manycov") {
      comb[["inst.func.family"]] <- "cdmc"
    }
    
    # Directory: search method, inst.func.family, K.bar
    sub.dir <- sprintf("Search-%s__IF-%s__Kbar-%s__Alpha-%s", comb[["search.method"]],
                       comb[["inst.func.family"]], comb[["K.bar"]], comb[["alpha"]])
    if ("link.function" %in% names(comb)) {
      sub.dir <- sprintf("LF-%s__%s", comb[["link.function"]], sub.dir)
    }
    if ("parametric" %in% names(comb)) { # Option used in misspecification simulations
      sub.dir <- sprintf("%s__Param-%s", sub.dir, comb[["parametric"]])
    }
    
    # If necessary, remove defaults
    if (shortened) {
      components <- str_split(sub.dir, "__")[[1]]
      components.shortened <- NULL
      for (component in components) {
        if (!(str_split(component, "-")[[1]][1] %in% names(defaults))) {
          components.shortened <- c(components.shortened, component)
        }
      }
      
      sub.dir <- paste(components.shortened, collapse = "__")
    }
    
    # Create directory if necessary
    check_create.dir(master.dir)
    check_create.dir(sub.dir, master.dir)
    
    # File name: (gs.method,) (n.if.per.cov,) n.sim, n, n.cov, DGP, B, t.eval
    file.name <- if (comb[["search.method"]] == "GS") {sprintf("method-%s", comb[["gs.method"]])} else NULL
    file.name <- if (comb[["inst.func.family"]] == "box") {
      paste(c(file.name, sprintf("nbpc-%s", comb[["n.if.per.cov"]])), collapse = "__")
    } else if (comb[["inst.func.family"]] == "spline") {
      paste(c(file.name, sprintf("nifpc-%s__degree-%s", comb[["n.if.per.cov"]], comb[["degree"]])), collapse = "__")
    } else if (comb[["inst.func.family"]] %in% c("cd", "cdmc")) {
      paste(c(file.name, sprintf("nifpc-%s__Gc-%s__degree-%s", comb[["n.if.per.cov"]], comb[["inst.func.family.c"]], comb[["degree"]])), collapse = "__")
    } else file.name
    gen.info <- sprintf("geninfo-%s", paste(c(comb[["n"]], comb[["n.cov"]], 
                                              comb[["DGP"]], comb[["B"]],
                                              comb[["t.eval"]], seed, sim.iter.nbr),
                                            collapse = "__"))
    file.name <- ifelse(is.null(file.name),
                        gen.info,
                        paste(c(file.name, gen.info), collapse = "__"))
    file.name <- paste0(file.name, ".csv")
    
    # If necessary, remove defaults
    if (shortened) {
      components <- str_split(file.name, "__")[[1]]
      components.shortened <- NULL
      for (component in components) {
        if (!(str_split(component, "-")[[1]][1] %in% names(defaults))) {
          components.shortened <- c(components.shortened, component)
        }
      }
      
      file.name <- paste(components.shortened, collapse = "__")
    }
  }
  
  # Obtain full name
  path <- paste(master.dir, sub.dir, file.name, sep = "/")
  
  # Return the result
  path
}

#' @title Get matrix of simulation settings
#' 
#' @description This function takes a list of simulation parameters as argument
#' and returns a matrix where each row in the matrix corresponds to a setting
#' that should be simulated. This is especially handy if we want to work with
#' the 'worker' module on the VSC. This function is used in 'run.simulations.R'.
#' 
#' @param sim.params List of simulation paramerers.
#' @param for.VSC Boolean value indicating whether the matrix returned by
#' this function is meant to be used as input to the worker module in the super-
#' computer. If \code{for.VSC = TRUE}, additional columns will be added to the
#' the matrix, the matrix will be saved and will not be returned. Default is
#' \code{for.VSC = FALSE}.
#' @param for.cCDF Boolean value indicating whether the matrix returned by this
#' function is meant to be used in the conditional CDF simulations. If
#' \code{for.cCDF = TRUE}, the matrix will be adapted so that the same seed is
#' used among desings that differ only in the value for \code{t.eval}. Default
#' is \code{for.cCDF = FALSE}.
#' 
#' @returns Simulation setting matrix
#' 
#' @noRd
#' 
get.simulation.variable.matrix <- function(sim.params, for.VSC = FALSE,
                                           for.cCDF = FALSE) {
  
  # Extract variables used in simulation for readability
  n.sim.vct <- sim.params[["n.sim"]]
  n.vct <- sim.params[["n"]]
  n.cov.vct <- sim.params[["n.cov"]]
  DGP.vct <- sim.params[["DGP"]]
  search.method.vct <- sim.params[["search.method"]]
  alpha.vct <- sim.params[["alpha"]]
  t.eval.vct <- sim.params[["t.eval"]]
  K.bar.vct <- sim.params[["K.bar"]]
  n.if.per.cov.vct <- sim.params[["n.if.per.cov"]]
  B.vct <- sim.params[["B"]]
  inst.func.family.vct <- sim.params[["inst.func.family"]]
  inst.func.family.c.vct <- sim.params[["inst.func.family.c"]]
  degree.vct <- sim.params[["degree"]]
  gs.method.vct <- sim.params[["gs.method"]]
  iseed <- sim.params[["iseed"]]
  idx.param.of.interest <- sim.params[["idx.param.of.interest"]]
  
  # Extract variables related to the data generation
  beta.true <- sim.params[["beta.true"]]
  par.space <- sim.params[["par.space"]]
  
  # Master directory for storing simulation results
  master.dir <- sim.params[["master.dir"]]
  
  # The grid search method might not have been specified when simulations will
  # only take place using the EAM algorithm.
  if (!("GS" %in% search.method.vct)) {
    gs.method.vct <- 1
  }
  
  # The degree of the spline instrumental functions might not have been specified
  # when the family of instrumental functions to be used are box functions.
  if ((length(inst.func.family.vct) == 1) & ("box" %in% inst.func.family.vct)) {
    degree.vct <- 0
  }
  
  # The family of instrumental functions to be used when G.cd is selected might
  # not be specified when G.cd is not selected.
  if (!("cd" %in% inst.func.family.vct)) {
    inst.func.family.c.vct <- "ignore"
  }
  
  # All simulation parameters must be specified
  if (any(is.null(c(n.sim.vct, n.vct, n.cov.vct, DGP.vct, search.method.vct,
                    alpha.vct, gs.method.vct, K.bar.vct, n.if.per.cov.vct,
                    inst.func.family.vct, degree.vct, inst.func.family.c.vct,
                    B.vct, t.eval.vct, iseed, dir, c)))) {
    stop("Some simulation parameters are not specified")
  }
  
  # If the matrix of combinations is meant to be supplied to the supercomputer...
  if (for.VSC) {
    
    # Re-code the string-valued options into integers (turns out to be not
    # necessary but I'll just leave it like this).
    search.method2int <- function(sm.str) {
      if (sm.str == "GS") {return(1)}
      if (sm.str == "EAM") {return(2)}
    }
    search.method.vct <- unname(Vectorize(search.method2int)(search.method.vct))
    inst.func.fam2int <- function(iff.str) {
      if (iff.str == "box") {return(1)}
      if (iff.str == "spline") {return(2)}
      if (iff.str == "cd") {return(3)}
      else {return(0)}
    }
    inst.func.family.vct <- unname(Vectorize(inst.func.fam2int)(inst.func.family.vct))
    inst.func.family.c.vct <- unname(Vectorize(inst.func.fam2int)(inst.func.family.c.vct))
    
    # For simplicity, only allow for one option for the number of simulations to
    # be run.
    if (length(n.sim.vct) != 1) {
      stop("n.sim.vct must have length 1.")
    }
    n.sim <- n.sim.vct
    
    # Obtain matrix of all parameter combinations to run
    combinations <- expand.grid(n.vct, n.cov.vct, DGP.vct, search.method.vct,
                                alpha.vct, gs.method.vct, K.bar.vct, n.if.per.cov.vct,
                                inst.func.family.vct, degree.vct, inst.func.family.c.vct,
                                B.vct, t.eval.vct)
    colnames(combinations) <- c("n", "n.cov", "DGP", "search.method", "alpha",
                                "gs.method", "K.bar", "n.if.per.cov",
                                "inst.func.family", "degree", "inst.func.family.c",
                                "B", "t.eval")
    
    # Initialize vector of seeds to be used
    seeds.to.use <- iseed + 1:n.sim
    
    # Each simulation setting should be run n.sim amount of times, all with a
    # different initial seed.
    comb.extended <- matrix(nrow = 0, ncol = ncol(combinations) + 1)
    colnames(comb.extended) <- c("seed", colnames(combinations))
    comb.extended <- as.data.frame((comb.extended))
    for (i in 1:nrow(combinations)) {
      for (seed in seeds.to.use) {
        row <- c(seed = seed, combinations[i, ])
        comb.extended <- rbind(comb.extended, row)
      }
      seeds.to.use <- seeds.to.use + n.sim
    }
    
    # When the search method is not 'grid search', we should not simulate over
    # different values for 'grid search method'.
    comb.extended[which(comb.extended$search.method != "GS"), "gs.method"] <- 1
    
    # When the family of instrumental functions is the box family, we should not
    # simulate over different degrees.
    comb.extended[which(comb.extended$inst.func.family == 1), "degree"] <- 0
    
    # When the family of instrumental functions is not the continuous/discrete
    # family, we should not iterate over different families for the continuous
    # part.
    comb.extended[which(comb.extended$inst.func.family != 3), "inst.func.family.c"] <- 0
    
    # When the family of instrumental functions is the continuous/discrete
    # family and the family for the continuous part are box functions, we should
    # not iterate over different values of degree.
    comb.extended[which((comb.extended$inst.func.family == 3) & (comb.extended$inst.func.family.c == 1)), "degree"] <- 0
    
    # When the degree of the instrumental functions is not smaller than the
    # number of instrumental functions to be used per covariate, the method will
    # throw an error. Remove such cases.
    comb.extended <- comb.extended[!(comb.extended$n.if.per.cov <= comb.extended$degree), ]
    
    # Given the previous modifications, obtain all unique simulation settings.
    comb.extended <- unique(comb.extended)
    
    # Rename columns
    colnames(comb.extended) <- c("seed", "n", "n_cov", "DGP", "search_method", "alpha",
                                 "gs_method", "K_bar", "n_if_per_cov",
                                 "inst_func_family", "degree", "inst_func_family_c",
                                 "B", "t_eval")
    
    # If necessary, adapt the matrix for use in cCDF simulations
    if (for.cCDF) {
      
      # Get all simulation designs: remove seed and t_eval information
      designs <- comb.extended[, !(colnames(comb.extended) %in% c("seed", "t_eval"))]
      
      # Get all unique simulation designs
      unique.designs <- unique(designs)
      
      # For each unique simulation design, set the same seed to be used over all
      # values of t.eval.
      for (design.idx in 1:nrow(unique.designs)) {
        
        # Select unique desing of this iteration
        design <- unique.designs[design.idx,]
        
        # Get row indices of all simulations corresponding to this design
        equal.design.idxs <- which(apply(designs, 1, function(row) {all(row == design)}))
        
        # For each, make the seed to be used equal
        seed.to.use <- comb.extended[min(equal.design.idxs), "seed"]
        comb.extended[equal.design.idxs, "seed"] <- seed.to.use
      }
    }
    
    # Save the results
    write.csv(comb.extended, "combinations.csv", row.names = FALSE)
  }
  
  # If the matrix of combinations is not meant to be supplied to the super-
  # computer...
  if (!for.VSC) {
    
    # Obtain matrix of all parameter combinations to run
    combinations <- expand.grid(n.sim.vct, n.vct, n.cov.vct, DGP.vct, search.method.vct,
                                alpha.vct, gs.method.vct, K.bar.vct, n.if.per.cov.vct,
                                inst.func.family.vct, degree.vct, inst.func.family.c.vct,
                                B.vct, t.eval.vct, stringsAsFactors = FALSE)
    colnames(combinations) <- c("n.sim", "n", "n.cov", "DGP", "search.method",
                                "alpha", "gs.method", "K.bar", "n.if.per.cov",
                                "inst.func.family", "degree", "inst.func.family.c",
                                "B", "t.eval")
    
    # When the search method is not 'grid search', we should not simulate over
    # different values for 'grid search method'.
    combinations[which(combinations$search.method != "GS"), "gs.method"] <- 1
    
    # When the family of instrumental functions is the box family, we should not
    # simulate over different degrees.
    combinations[which(combinations$inst.func.family == "box"), "degree"] <- 0
    
    # When the family of instrumental functions is not the continuous/discrete
    # family, we should not iterate over different families for the continuous
    # part.
    combinations[which(combinations$inst.func.family != "cd"), "inst.func.family.c"] <- "NA"
    
    # When the family of instrumental functions is the continuous/discrete
    # family and the family for the continuous part are box functions, we should
    # not iterate over different values of degree.
    combinations[which((combinations$inst.func.family == "cd") & (combinations$inst.func.family.c == "box")), "degree"] <- 0
    
    # Given the previous modifications, obtain all unique simulation settings.
    combinations <- unique(combinations)
    
    # Add the initial seeds to be used
    iseed.col <- rep(iseed, nrow(combinations))
    n.sim <- 0
    for (i in 1:nrow(combinations)) {
      iseed <- iseed + n.sim
      iseed.col[i] <- iseed
      n.sim <- combinations[i, "n.sim"]
    }
    iseed <- iseed.col
    combinations <- cbind(combinations, iseed)
    
    # If necessary, adapt the matrix for use in cCDF simulations
    if (for.cCDF) {
      
      # Get all simulation designs: remove seed and t_eval information
      designs <- combinations[, !(colnames(combinations) %in% c("iseed", "t.eval"))]
      
      # Get all unique simulation designs
      unique.designs <- unique(designs)
      
      # For each unique simulation design, set the same seed to be used over all
      # values of t.eval.
      for (design.idx in 1:nrow(unique.designs)) {
        
        # Select unique desing of this iteration
        design <- unique.designs[design.idx,]
        
        # Get row indices of all simulations corresponding to this design
        equal.design.idxs <- which(apply(designs, 1, function(row) {all(row == design)}))
        
        # For each, make the seed to be used equal
        seed.to.use <- combinations[min(equal.design.idxs), "iseed"]
        combinations[equal.design.idxs, "iseed"] <- seed.to.use
      }
    }
    
    # Return the results
    return(combinations)
  }
}

#' @title Run a simulation for a specific set of variables
#' 
#' @description This function is used in 'run.simulations.R'
#' 
#' @param comb Data frame of one single row containing the parameter settings to
#' be used in the simulation.
#' @param beta.true True values of the parameters in the model; used to generate
#' the data.
#' @param idx.param.of.interest Index of the element of the parameter vector one
#' wants to create the identified interval for. I.e. if c is the projection
#' vector, then c[idx.param.of.interest] = 1, and c = 0 elsewhere.
#' @param par.space Bounds on the parameter space.
#' @param starting.seed Initial random seed to use
#' @param master.dir Name of the directory in which the simulation results
#' should be stored.
#' @param verbose Verbosity parameter.
#' 
#' @note [ToDo] Also allow for simulation using other link functions
#' 
#' @noRd
#' 
simulate1D <- function(comb, beta.true, idx.param.of.interest, par.space,
                       starting.seed, master.dir, verbose) {
  
  #### Extract hyperparameters ####
  
  # Extract the parameters (currently, cov.ranges will never be in 'comb')
  n.sim <- comb[, "n.sim"]
  n <- comb[, "n"]
  n.cov <- comb[, "n.cov"]
  DGP <- comb[, "DGP"]
  search.method <- comb[, "search.method"]
  alpha <- comb[, "alpha"]
  t.eval <- comb[, "t.eval"]
  K.bar <- comb[, "K.bar"]
  n.if.per.cov <- comb[, "n.if.per.cov"]
  B <- comb[, "B"]
  inst.func.family <- comb[, "inst.func.family"]
  degree <- comb[, "degree"]
  inst.func.family.c <- comb[, "inst.func.family.c"]
  gs.method <- comb[, "gs.method"]
  cov.ranges <- if ("cov.ranges" %in% colnames(comb)) {comb[, "cov.ranges"]} else NULL
  
  # Precondition checks. If at some point this would be adapted, also adapt the
  # code in 'conditionalCDFFunctions.R -> get.Peterson.bounds.R'
  beta0.checks <- simplify2array(lapply((-10):10, function(x) {beta.true(x)[1]}))
  if (any(beta0.checks != (-10):10)) {
    stop("Currently only able to handle beta0 = t. (*)")
    # Technically, beta0 = t + b for any b would also work.
  }
  
  # Projection vector that projects the parameter vector onto the element of
  # interest
  c <- rep(0, n.cov + 1)
  c[idx.param.of.interest] <- 1
  
  # Family of instrumental functions to be used when G.cd is selected
  G.c <- NULL
  if (inst.func.family == "cd") {
    if (inst.func.family.c == "box") {
      G.c <- G.box
    } else if (inst.func.family.c == "spline") {
      G.c <- G.spline
    } else {
      stop("Specified family of instrumental functions not implemented.")
    }
  }
  
  # If applicable, select the type of grid search to be carried out
  if (gs.method == 1) {
    next.gs.point = gs.binary  # Binary search
  } else if (gs.method == 2) {
    next.gs.point = gs.regular # Regular grid search
  }
  
  # List of hyperparameters to be supplied to various functions later on
  options <- list(n.if.per.cov = n.if.per.cov,
                  K.bar = K.bar,
                  B = B,
                  gs.method = gs.method,
                  next.gs.point = next.gs.point,
                  alpha = alpha,
                  link.function = "AFT_ll",
                  DGP = DGP,
                  inst.func.family = inst.func.family,
                  degree = degree,
                  G.c = G.c,
                  cov.ranges = cov.ranges)
  
  # Flag for recording the time that each step of the algorithm takes
  time.run.duration <- TRUE
  
  #### Perform the simulations ####
  
  # Initialize object that will store the results
  ncol.sim.results <- 7 + 3 * as.numeric(search.method == "EAM") +
                      2 * as.numeric(DGP %% 20 %in% c(10, 11, 12))
  sim.results <- matrix(nrow = 0, ncol = ncol.sim.results)
  cols <- c("ident.set.l", "ident.set.u", "conv.l", "conv.u", "total.run.time")
  if (search.method == "EAM") {
    cols <- c(cols, "E-step", "A-step", "M-step")
  }
  cols <- c(cols, "seed", "per.cens")
  if ((DGP %% 20) %in% c(10, 11, 12)) {
    cols <- c(cols, "per.cens.reg1", "per.cens.reg2")
  }
  colnames(sim.results) <- cols
  
  # Run the simulation
  for (sim.iter.nbr in 1:n.sim) {
    if (verbose >= 1) {
      message(sprintf("Starting simulation %d / %d", sim.iter.nbr, n.sim))
    }
    
    # Set random seed of this iteration
    seed.to.use <- starting.seed + sim.iter.nbr
    set.seed(seed.to.use)
    
    # Generate data 
    data <- generateData(beta.true, n, n.cov, options, plot.data = (verbose >= 3))
    
    # Find the identified set based on Bei (2024) 
    fis.out <- find.identified.set(c, t.eval, par.space, data, search.method,
                                   options, verbose, time.run.duration)
    ident.set <- fis.out$ident.set
    
    # Analyse the running times
    chrono1 <- fis.out$chronometer1
    chrono2 <- fis.out$chronometer2
    total.run.time <- chrono1$get.total.time(force = TRUE) + chrono2$get.total.time(force = TRUE)
    if (search.method == "EAM") {
      total.leg.times <- chrono1$accumulate.legs(force = TRUE)[1:3] +
                         chrono2$accumulate.legs(force = TRUE)[1:3]
    } else {
      total.leg.times <- NULL
    }
    
    # Convergence of the algorithm
    conv.l <- fis.out$converge2
    conv.u <- fis.out$converge1
    
    # Censoring percentage in the data
    per.cens <- 1 - sum(data$Delta)/n
    per.cens.region1 <- NULL
    per.cens.region2 <- NULL
    if ((DGP %% 20) %in% c(10, 11, 12)) {
      per.cens.region1 <- 1 - sum(data$Delta[inRegion1(data)])/length(inRegion1(data))
      per.cens.region2 <- 1 - sum(data$Delta[inRegion2(data)])/length(inRegion2(data))
    }
    
    # Store the results
    sim.result <- c(ident.set, conv.l, conv.u, total.run.time, total.leg.times,
                    seed.to.use, per.cens, per.cens.region1, per.cens.region2)
    sim.results <- rbind(sim.results, sim.result)
    
    # Store (also intermediate) results
    path <- get.file.name(comb, sim.iter.nbr, seed.to.use, master.dir)
    write.csv(sim.results, path, row.names = FALSE)
  }
  
  # Return the number of simulations ran (used in updating the initial seed for
  # the next simulation)
  n.sim
}

#' @title Run simulation for all combinations of given parameters.
#' 
#' @description This function runs the specified simulations.
#' 
#' @param sim.params List of parameters to be used in the simulation. For each
#' parameter, a vector of values can be provided. The function will run the
#' simulation for each combination of values for the parameters that can be
#' created.
#' @param verbose Integer value indicating the degree to which the function
#' should present run time updates to the user.
#' 
#' @noRd
#' 
run.simulations <- function(sim.params, verbose) {
  
  #### Obtain matrix for variables of the simulation ####
  
  combinations <- get.simulation.variable.matrix(sim.params)
  
  # Update the user
  if (verbose >= 1) {
    message("\n \n")
    message(sprintf("--- Starting simulation over the %d combinations ---",
                    nrow(combinations)))
    message("\n \n")
  }
  
  #### Run the simulations ####
  
  # Initialize some variables
  iseed <- sim.params[["iseed"]]
  seed.to.use <- iseed
  
  # Run the simulations
  foreach(row.idx = 1:nrow(combinations)) %dopar% {
    
    # Load dependencies
    source("simulationFunctions.R")
    
    # Update the user
    if (verbose >= 1) {
      message(sprintf("Simulation configuration %d / %d", row.idx,
                      nrow(combinations)))
    }
    
    # Increment the seed to be used
    if (row.idx > 1) {
      seed.to.use <- iseed + sum(combinations[1:(row.idx - 1), "n.sim"])
    }
    
    # Get the parameter values corresponding to this iteration
    comb <- combinations[row.idx, ]
    
    # Get true values for coefficients, alongside their parameter space
    beta.true <- sim.params[["beta.true"]]
    par.space <- sim.params[["par.space"]]
    
    # Get index of the parameter of interest
    idx.param.of.interest <- sim.params[["idx.param.of.interest"]]
    
    # Get name of master directory
    master.dir <- sim.params[["master.dir"]]
    
    # Run the simulation.
    simulate1D(comb, beta.true, idx.param.of.interest, par.space, iseed,
               master.dir, verbose)
  }
}

#' @title Estimate an identified model assuming independence.
#' 
#' @description This function estimates an identified model assuming that the
#' event and censoring time independent. Specifically, it estimates a regular
#' AFT model or Cox PH model. This is used as a reference in the simulations. 
#' 
#' @param c Projection vector
#' @param t.eval Time point of interest.
#' @param data Data frame.
#' @param link.function Name of the link function to use ("AFT_ll" or "Cox_wb").
#' 
#' @returns The coefficients of the identified model.
#' 
#' @import survival
#' 
#' @noRd
#' 
get.coef.identified.model <- function(c, t.eval, data, link.function) {
  
  # Get number of covariates
  n.cov <- sum(grepl("X[1-9][[:digit:]]*$", colnames(data)))
  
  # Estimate the model
  if (link.function == "AFT_ll") {
    fml <- as.formula(paste0("Surv(Y, Delta) ~ ", paste(paste0("X", 1:n.cov), collapse = " + ")))
    fit <- survreg(fml, dist = "loglogistic", data = data)
    coef <- (-fit$coefficients) %*% c
    
  } else if (link.function == "Cox_wb") {
    fml <- as.formula(paste0("Surv(Y, Delta) ~ ", paste(paste0("X", 1:n.cov), collapse = " + ")))
    fit <- coxph(fml, data = data)
    coef <- fit$coefficients %*% c[-1]
  }
  
  # Return the coefficient of interest
  coef
}



