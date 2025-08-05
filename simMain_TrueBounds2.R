
#### Preamble ####

on.VSC <- grepl("vsc36272", getwd())

# Load dependencies
if (on.VSC) {
  
  # Load dependencies
  Sys.setenv(TZ = "Europe/Amsterdam")
  lib.loc <- paste0(dirname(getwd()), "/R/")
  pkg.list <- c("MASS", "lubridate", "SPOT", "EnvStats", "splines2", "copula",
                "gplm", "ks", "PLRModels", "np", "tidyr", "tzdb", "readr",
                "forcats", "tidyverse", "survival", "doParallel", "nloptr")
  for (pkg in pkg.list) {
    eval(bquote(require(.(pkg), lib.loc = lib.loc)))
  }
  source("simulationFunctions.R")
  
  # Set master directory
  master.dir <- "Simulations_Main_trueBounds"
  check_create.dir(master.dir)
  
  # Arguments supplied to the function
  args <- commandArgs(TRUE)
  
} else {
  
  # Clear workspace
  rm(list = setdiff(ls(), "on.VSC"))
  
  require("MASS")
  require("lubridate")
  require("SPOT")
  require("EnvStats")
  require("splines2")
  require("copula")
  source("simulationFunctions.R")
  
  # Set master directory
  master.dir <- "Simulations_Main_trueBounds"
  check_create.dir(master.dir)
  
  # Dummy argument vector
  args <- c(26)
}
# IMPORTANT NOTE:
# This script computes the true identified interval FOR A GIVEN DATA SET.
# Indeed, note that the instrumental functions are data-driven, and that the 
# true identified set depends on these instrumental functions.

# Store starting time
time.start <- Sys.time()

#### Create matrix of DGPs to check ####

if (!on.VSC) {
  n.reps <- 4
  seed.inits <- 1000*(1:n.reps)
  DGP.options <- as.vector(outer(seed.inits, 21:32, FUN = "+"))
  designs <- matrix(DGP.options, ncol = 1)
  colnames(designs) <- c("DGP")
  write.csv(designs, "simMain_trueBoundsDesigns.csv", row.names = FALSE)
}

#### Extract arguments ####

seed <- as.numeric(args[1])
DGP <- seed %% 1000

# Use a different seed for each DGP
set.seed(seed)

#### Define the settings of the main simulation function ####

# Time point of interest
t <- 1

# Sample size
n.on.VSC <- 50000
n <- n.on.VSC 

# Number of covariates
n.cov <- 2

# True parameter vectors
H0 <- function(t) {t}
H0.inv <- function(t) {t}
beta1 <- 1
beta2 <- -1
beta3 <- 2
beta4 <- -1
beta.true <- function(t) {c(H0(t), beta1, beta2, beta3, beta4)}
par.space.lb <- rep(-10, n.cov + 1)
par.space.ub <- rep(10, n.cov + 1)
par.space <- matrix(c(par.space.lb, par.space.ub), ncol = 2)

# Parameter of interest
idx.param.of.interest <- 2
c <- rep(0, n.cov + 1)
c[idx.param.of.interest] <- 1

# Other hyperparameters
n.if.per.cov <- 6
K.bar <- 3
B <- 600
search.method <- "GS"
next.gs.point <- gs.binary
alpha <- 0.95
link.function <- ifelse(DGP %in% c(21, 22, 23, 24, 25, 26), "AFT_ll",
                 ifelse(DGP %in% c(27, 28, 29, 30, 31, 32), "Cox_wb",
                 stop("Unknown DGP")))
inst.func.family <- "cd"
degree <- 3
G.c <- G.spline
cov.ranges <- NULL
options <- list(n.if.per.cov = n.if.per.cov,
                K.bar = K.bar,
                B = B,
                search.method = search.method,
                next.gs.point = next.gs.point,
                alpha = alpha,
                link.function = link.function,
                DGP = DGP,
                inst.func.family = inst.func.family,
                degree = degree,
                G.c = G.c,
                cov.ranges = cov.ranges)

# Generate data
data <- generateData_simMain(beta.true, n, n.cov, options, H0.inv,
                             plot.data = FALSE)

# Set hyperparameters
hp <- set.hyperparameters(data, par.space, c, search.method, options)

#### Define all relevant functions ####

# Number of instrumental functions
J <- hp$n.inst.func

# Instrumental functions
g <- function(x_noint, j) {
  hp$G(x_noint, j)
}

# Link function
Lambda <- function(lin.comb) {
  hp$Lambda(lin.comb)
}

# Unconditional moment functions
get.g.evals.mat <- function(x) {
  x_noint <- x[,-1]
  g.evals.mat <- matrix(nrow = nrow(x), ncol = J)
  for (j in 1:J) {
    g.evals.mat[, j] <- apply(x_noint, 1, g, j = j)
  }
  g.evals.mat
}
m1_j <- function(y, delta, x, t, j, beta, g.evals) {
  (as.numeric(y <= t) - Lambda(x %*% beta)) * g.evals
}
m2_j <- function(y, delta, x, t, j, beta, g.evals) {
  (Lambda(x %*% beta) - as.numeric(y <= t & delta == 1)) * g.evals
}
m1 <- function(y, delta, x, t, beta, g.evals.mat) {
  res <- matrix(nrow = length(y), ncol = J)
  for (j in 1:J) {
    g.evals <- g.evals.mat[, j]
    res[, j] <- m1_j(y, delta, x, t, j, beta, g.evals)
  }
  res
}
m2 <- function(y, delta, x, t, beta, g.evals.mat) {
  res <- matrix(nrow = length(y), ncol = J)
  for (j in 1:J) {
    g.evals <- g.evals.mat[, j]
    res[,j] <- m2_j(y, delta, x, t, j, beta, g.evals)
  }
  res
}
m <- function(y, delta, x, t, beta, g.evals.mat) {
  cbind(m1(y, delta, x, t, beta, g.evals.mat),
        m2(y, delta, x, t, beta, g.evals.mat))
}
m.datawrapper <- function(MC.data, t, beta, g.evals.mat) {
  y <- MC.data[, "Y"]
  delta <- MC.data[, "Delta"]
  x <- as.matrix(MC.data[, grepl("X[[:digit:]]+", colnames(MC.data))])
  m(y, delta, x, t, beta, g.evals.mat)
}

# Monte Carlo integration function
integrate.MC <- function(t, beta, g.evals.mat = NULL) {
  
  # Generate data according to true DGP
  # MC.data <- generateData_simMain(beta.true, sam.size, n.cov, options, H0.inv,
  #                                 plot.data = FALSE)
  MC.data <- data
  
  # Compute instrumental function evaluations
  if (is.null(g.evals.mat)) {
    x <- as.matrix(MC.data[, grepl("X[[:digit:]]+", colnames(MC.data))])
    g.evals.mat <- get.g.evals.mat(x)
  }
  
  # Compute MC sample moments
  evals <- m.datawrapper(MC.data, t, beta, g.evals.mat)
  colMeans(evals)
}

# Remove duplicate rows and order matrix
order.matrix <- function(m) {
  m <- round(m, 3)
  m <- unique(m)
  call <- sprintf("order(%s)", paste(paste("m[, ", 1:ncol(m), "]", sep = ""), collapse = ", "))
  m[eval(parse(text = call)), , drop = FALSE]
}

# Check if two vectors are equal up to tolerance
compare.vecs <- function(vec1, vec2, err.tol) {
  diffvec <- vec1 - vec2
  diffvec.disc <- ifelse(diffvec < -err.tol, -1, ifelse(diffvec > err.tol, 1, 0))
  if (all(diffvec.disc == 0)) {
    return(0)  # Equal
  } else if (diffvec.disc[diffvec.disc != 0][1] == 1) {
    return(1)  # vec1 larger
  } else {
    return(-1) # vec1 smaller
  }
}

# Compute set difference based on ordered matrices up to tolerance
set.diff <- function(m1, m2) {
  dupe.idxs <- NULL
  m1.row.idx <- 1
  m2.row.idx <- 1
  while ((m1.row.idx <= nrow(m1)) & (m2.row.idx <= nrow(m2))) {
    vec1 <- m1[m1.row.idx, ]
    vec2 <- m2[m2.row.idx, ]
    comp <- compare.vecs(vec1, vec2, err.tol = 1e-8)
    if (comp == 0) {
      dupe.idxs <- c(dupe.idxs, m1.row.idx)
      m1.row.idx <- m1.row.idx + 1
      m2.row.idx <- m2.row.idx + 1
    } else if (comp == 1) {
      m2.row.idx <- m2.row.idx + 1
    } else if (comp == -1) {
      m1.row.idx <- m1.row.idx + 1
    }
  }
  
  if (!is.null(dupe.idxs)) {
    m1 <- m1[-dupe.idxs, , drop = FALSE]
  }
  
  order.matrix(m1)
}

# Compute set union based on matrices, of which m2 is ordered
set.union <- function(m1, m2) {
  
  # If m1 and m2 are empty, return empty matrix
  if (nrow(m1) == 0 & nrow(m2) == 0) {
    return(m1)
  }
  
  # If m1 is empty, return sorted m2
  if (nrow(m1) == 0) {
    return(order.matrix(m2))
  }
  
  # If m2 is empty, return sorted m1
  if (nrow(m2) == 0) {
    return(order.matrix(m1))
  }
  
  # Order m1
  m1 <- order.matrix(m1)
  
  # Iterate over rows in m1, inserting them into m2 when appropriate.
  m1.row.idx <- 1
  m2.row.idx <- 1
  while (m1.row.idx <= nrow(m1)) {
    
    # Get vector in m1 of this iteration
    vec1 <- m1[m1.row.idx, ]
    
    if (nrow(m2) > 1) {
      
      # Find vectors vec2 and vec3 in m2 such that vec2 < vec1 < vec3
      vec1.inserted <- FALSE
      while (m2.row.idx < nrow(m2) & !vec1.inserted) {
        
        # Get vec2 and vec3. Compare to vec1.
        vec2 <- m2[m2.row.idx, ]
        vec3 <- m2[m2.row.idx + 1, ]
        comp1 <- compare.vecs(vec1, vec2, err.tol = 1e-8)
        comp2 <- compare.vecs(vec1, vec3, err.tol = 1e-8)
        
        # Decide on next action
        if (comp1 == 1 & comp2 == -1) { # If vec2 < vec1 < vec3
          
          # Insert vec1 into m2
          dummy <- matrix(nrow = nrow(m2) + 1, ncol = ncol(m2))
          dummy[-(m2.row.idx + 1), ] <- m2
          dummy[m2.row.idx + 1, ] <- vec1
          m2 <- dummy
          
          # Set-up for next iteration of outer while
          m2.row.idx <- m2.row.idx + 1
          vec1.inserted <- TRUE
          
        } else if (comp1 == 0 | comp2 == 0) { # If vec1 already in m2
          
          # Break loop
          vec1.inserted <- TRUE
          
        } else if (comp1 == -1) { # If vec1 < vec2
          
          # Insert vec1 into m2
          dummy <- matrix(nrow = nrow(m2) + 1, ncol = ncol(m2))
          dummy[-m2.row.idx, ] <- m2
          dummy[m2.row.idx, ] <- vec1
          m2 <- dummy
          
          # Set-up for next iteration of outer while
          vec1.inserted <- TRUE
          
        } else {
          
          # Set up for next iteration if inner while
          m2.row.idx <- m2.row.idx + 1
        }
      }
      
      if (m2.row.idx == nrow(m2) & !vec1.inserted) {
        m2 <- rbind(m2, vec1)
        m2.row.idx <- nrow(m2)
      }
      
    } else {
      vec2 <- m2[m2.row.idx, ]
      comp1 <- compare.vecs(vec1, vec2, err.tol = 1e-8)
      if (comp1 == 1) {
        m2 <- rbind(m2, vec1)
        m2.row.idx <- nrow(m2)
      } else if (comp1 == -1) {
        m2 <- rbind(vec1, m2)
        m2.row.idx <- 1
      }
    }
    
    # Increment m1.idx
    m1.row.idx <- m1.row.idx + 1
  }
  
  # Return the result
  order.matrix(m2)
}

# Obtain points to check in refinement stage
get.beta.to.test <- function(feas.beta, infeas.beta, by.val) {
  
  # If feas.beta is empty, return no new points to test
  if (nrow(feas.beta) == 0) {
    return(feas.beta)
  }
  
  # Intialize matrix of beta to check
  beta.to.check <- matrix(nrow = 0, ncol = ncol(feas.beta))
  
  # Matrix of directions in which to expand the grid at a given point
  call <- sprintf("expand.grid(%s)", paste(rep("c(-1, 0, 1)", ncol(feas.beta)), collapse = ", "))
  combinations <- as.matrix(eval(parse(text = call)))
  combinations[, 1] <- 0
  combinations <- unique(combinations)
  combinations <- combinations[!apply(combinations, 1, function (row) {all(row == 0)}), ]
  
  # Obtain beta values in a grid around feas.beta
  for (row.idx in 1:nrow(feas.beta)) {
    
    # Obtain beta of this iteration
    fb <- feas.beta[row.idx, , drop = FALSE]
    
    # Obtain beta.to.check of this iteration
    beta.to.check.fb <- matrix(nrow = 0, ncol = ncol(feas.beta))
    
    # Get all points in a grid around this value
    for (comb.idx in 1:nrow(combinations)) {
      comb <- as.numeric(combinations[comb.idx, , drop = FALSE])
      beta.to.check.fb <- set.union(beta.to.check.fb, round(fb + by.val * comb, 3))
    }
    
    # Remove already checked points
    beta.to.check.fb <- set.diff(beta.to.check.fb, feas.beta)
    
    # Add points to beta.to.check
    beta.to.check <- set.union(beta.to.check, beta.to.check.fb)
  }
  
  # Remove infeasible points
  beta.to.check <- set.diff(beta.to.check, infeas.beta)
  
  # Return results
  order.matrix(beta.to.check)
}

# Obtain time difference
get.time.diff <- function(t.start, t.end) {
   h <- floor(difftime(t.end, t.start, units = 'hours'))
   m <- as.numeric(floor(difftime(t.end, t.start, units = 'mins'))) %% 60
   s <- as.numeric(floor(difftime(t.end, t.start, units = 'secs'))) %% 60
   return(sprintf("%02d:%02d:%02d", h, m ,s))
}

#### Find identified set ####

# Search hyperparameters
by.val <- 0.1
beta.grid <- seq(-10, 10, by = by.val)
beta0.grid <- beta.grid
beta1.grid <- seq(beta1 - 3*by.val, beta1 + 3*by.val, by = by.val)
beta2.grid <- seq(beta2 - 3*by.val, beta2 + 3*by.val, by = by.val)

# Precompute instrumental function evaluations
g.evals.mat.file <- sprintf("%s/g.evals.mat%d.Rdata", master.dir, seed)
x <- as.matrix(data[, grepl("X[[:digit:]]+", colnames(data))])
if (file.exists(g.evals.mat.file)) {
  message(sprintf("Seed %d: Precomputing instrumental functions skipped...", seed))
  load(g.evals.mat.file)
} else {
  message(sprintf("Seed %d: Precomputing instrumental function evaluations...", seed))
  g.evals.mat <- get.g.evals.mat(x)
}
end.time.precompute <- Sys.time()
message(sprintf("Seed %d: Precomputing finished (runtime = %s)",
                seed, get.time.diff(time.start, end.time.precompute)))

# Save intermediate results
save(g.evals.mat, file = sprintf("%s/g.evals.mat%d.Rdata", master.dir, seed))

# File name of resulting eval.array
eval.array.file <- sprintf("%s/eval_array%d_intermediate.Rdata", master.dir, seed)

# Check if initial grid search has already been performed. If so, load the
# results and skip to the next step.
if (file.exists(eval.array.file)) {
  load(eval.array.file)
  message(sprintf("Seed %d: Initial grid search skipped", seed))
  
# If not, do the initial grid search.
} else {
  eval.array <- array(Inf, dim = c(length(beta0.grid), length(beta1.grid), length(beta2.grid)))
  message(sprintf("Seed %d: Starting identified set search.", seed))
  for (beta0.idx in 1:length(beta0.grid)) {
    for (beta1.idx in 1:length(beta1.grid)) {
      for (beta2.idx in 1:length(beta2.grid)) {
        beta0 <- beta0.grid[beta0.idx]
        beta1 <- beta1.grid[beta1.idx]
        beta2 <- beta2.grid[beta2.idx]
        
        beta <- c(beta0, beta1, beta2)
        eval.array[beta0.idx, beta1.idx, beta2.idx] <- sum(pmin(0, integrate.MC(t, beta, g.evals.mat))^2)
      }
    }
    
    # Store intermediate results
    save(eval.array, file = eval.array.file)
    
    # Update user
    if ((beta0.idx %% (length(beta0.grid) %/% 10)) == 0) {
      message(sprintf("Seed %d: %.2f%% completion. ",
                      seed, 100 * (beta0.idx - 1)/length(beta.grid)),
              appendLF = FALSE)
      message(sprintf("runtime so far = %s.",
                      get.time.diff(end.time.precompute, Sys.time())))
    }
  }
}

# Plot the results, if appropriate
if (!on.VSC & (length(beta1.grid) == 1) & (length(beta2.grid) == 1)) {
  plot(beta0.grid, eval.array[, 1, 1])
}

# Refinement stage
message(sprintf("Seed %d: Starting refinement stage.", seed))
min.by.val <- 0.05
tol <- 0
if (DGP == 22) {
  tol <- 1e-8
}
feas.beta.idxs <- arrayInd(which(eval.array <= tol), .dim = dim(eval.array))
if (nrow(feas.beta.idxs) == 0) {
  stop(sprintf("No initial feasible point found (Seed = %d).", seed))
}
feas.beta <- t(matrix(apply(feas.beta.idxs, 1, function(row) {c(beta0.grid[row[1]], beta1.grid[row[2]], beta2.grid[row[3]])}),
                    ncol = nrow(feas.beta.idxs)))
feas.beta <- order.matrix(feas.beta)
infeas.beta.idxs <- arrayInd(which(eval.array > tol), .dim = dim(eval.array))
infeas.beta <- t(matrix(apply(infeas.beta.idxs, 1, function(row) {c(beta0.grid[row[1]], beta1.grid[row[2]], beta2.grid[row[3]])}),
                        ncol = nrow(infeas.beta.idxs)))
infeas.beta <- order.matrix(infeas.beta)

# Initialize matrix that will store all betas that were tested in the refinement
# stage which turned out to be infeasible.
beta.to.test.infeas <- matrix(nrow = 0, ncol = ncol(feas.beta))

# Load latest checkpoint
if (file.exists(sprintf("%s/eval_array%d_refined.Rdata", master.dir, seed))) {
  message(sprintf("Seed %d: Loading latest checkpoint of refinement stage...",
                  seed))
  load(sprintf("%s/eval_array%d_refined.Rdata", master.dir, seed))
}

while (by.val >= min.by.val) {
  
  # Update the user
  message(sprintf("Seed %d: Using by.val = %.3f.", seed, by.val))
  
  # Update infeas.beta
  infeas.beta <- set.union(infeas.beta, beta.to.test.infeas)
  
  # Find beta to test
  message(sprintf("Seed %d: Finding beta.to.test.with by.val = %.2f.",
                  seed, by.val))
  beta.to.test <- get.beta.to.test(feas.beta, infeas.beta, by.val)
  message(sprintf("Seed %d: Found beta.to.test (%d rows, by.val = %.2f).",
                  seed, nrow(beta.to.test), by.val))
  
  # While there are still beta left to test...
  while (nrow(beta.to.test) > 0) {
    
    # Plot results
    if (!on.VSC) {
      plot(feas.beta[, 2], feas.beta[, 3], xlim = c(0, 2), ylim = c(-3, 1), xlab = "beta1", ylab = "beta2")
      Sys.sleep(0)
      Sys.sleep(2)
      points(beta.to.test[, 2], beta.to.test[, 3], col = "orange")
    }
    
    # Initialize object that will store new beta to test
    beta.to.test.feas <- matrix(nrow = 0, ncol = ncol(beta.to.test))
    
    # Test each beta.
    for (beta.to.test.idx in 1:nrow(beta.to.test)) {
      
      # Extract beta of this iteration
      b <- beta.to.test[beta.to.test.idx, , drop = FALSE]
      b.vec <- beta.to.test[beta.to.test.idx, ]

      # Evaluate beta
      b.eval <- sum(pmin(0, integrate.MC(t, b.vec, g.evals.mat))^2)
      
      # Depending on whether it is a feasible point, store it as such. If feasible,
      # add grid of points around it to beta.to.test.
      if (b.eval <= tol) {
        beta.to.test.feas <- set.union(beta.to.test.feas, b)
        feas.beta <- set.union(feas.beta, b)
      } else {
        beta.to.test.infeas <- set.union(beta.to.test.infeas, b)
      }
    }
    
    # Plot results
    if (!on.VSC) {
      points(beta.to.test.feas[, 2], beta.to.test.feas[, 3], col = "green")
      points(beta.to.test.infeas[, 2], beta.to.test.infeas[, 3], col = "red")
    }
    
    # Once all beta to be tested are tested, obtain new batch of beta to test
    message(sprintf("Seed %d: Finding beta.to.test.with by.val = %.2f.",
                    seed, by.val))
    beta.to.test.new <- get.beta.to.test(beta.to.test.feas, infeas.beta, by.val)
    message(sprintf("Seed %d: Found beta.to.test (%d rows, by.val = %.2f).",
                    seed, nrow(beta.to.test), by.val))
    beta.to.test.new <- set.diff(beta.to.test.new, feas.beta)
    beta.to.test.new <- set.diff(beta.to.test.new, beta.to.test.infeas)
    beta.to.test <- beta.to.test.new
    
    # Plot results
    if (!on.VSC) {
      points(beta.to.test.new[, 2], beta.to.test.new[, 3], col = "blue")
    }
    
    # Update user
    message(sprintf("Seed %d: Saving intermediate refined results. Length feas.beta = %d. Length beta.to.test = %d.",
                    seed, nrow(feas.beta), nrow(beta.to.test)))
    
    # Save intermediate refined values
    save.image(file = sprintf("%s/eval_array%d_refined.Rdata", master.dir, seed))
  }
  
  # Decrease by.val
  by.val <- by.val/2
}

# Save refined values
save.image(file = sprintf("%s/eval_array%d_refined.Rdata", master.dir, seed))




