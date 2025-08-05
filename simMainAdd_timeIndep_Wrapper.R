
#### Preamble ####

on.VSC <- grepl("vsc36272", getwd())
makeDesigns <- FALSE

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
  master.dir <- "Simulations_Add_timeIndep"
  check_create.dir(master.dir)
  
  # Arguments supplied to the function
  args <- commandArgs(TRUE)
  
} else {
  
  # Clear workspace
  rm(list = ls())
  
  require("MASS")
  require("lubridate")
  require("SPOT")
  require("EnvStats")
  require("splines2")
  require("copula")
  source("simulationFunctions.R")
  
  # Set master directory
  master.dir <- "Simulations add timeIndep"
  check_create.dir(master.dir)
  
  # Arguments supplied to the function (for testing purposes)
  args <- c(6000141,1000,8,27,2,2)
}

if (makeDesigns) {
  
  # Set simulation settings
  n.to.check <- c(1000, 2000)
  n.if.per.cov.to.check <- c(8)
  DGP.to.check <- 27:32
  link.function.idx.to.check <- c(2)
  t.vec.idx.to.check <- 1:3
  
  # Make simulation design matrix
  designs <- expand.grid(n.to.check,
                         n.if.per.cov.to.check,
                         DGP.to.check,
                         link.function.idx.to.check,
                         t.vec.idx.to.check)
  designs <- cbind(6000123 + 1:nrow(designs), designs)
  colnames(designs) <- c("seed", "n", "n_if_per_cov", "DGP", "link_function_idx",
                         "t_vec_idx")
  
  # Save the result
  write.csv(designs, "simMainDesign1.csv", row.names = FALSE)
}

#### Extract and set parameters ####

# Extract arguments for simulation design
seed <- as.numeric(args[1])
n <- as.numeric(args[2])
n.if.per.cov <- as.numeric(args[3])
DGP <- as.numeric(args[4])
link.function.idx <- as.numeric(args[5])
t.vec.idx <- as.numeric(args[6])

# Select the correct vector of time points to take up in the model
t.vec.list <- list(
  c(1, 3, 5),
  c(0.3333, 0.6666, 1),
  c(2, 4, 6, 8, 10)
)
t <- t.vec.list[[t.vec.idx]]

# Set fixed settings
n.cov <- 2                       # Number of covariates
K.bar <- 3                       # Numer of critical value refinements.
# Recommendation in Bei (2024) (somewhere between 2 and 10)
alpha <- 0.95                    # 1 - Level of the test
search.method.idx <- 1           # Do not use EAM algorithm.
B <- 600                         # Number of bootstraps. Recommendation in Bei (2024)
inst.func.family.idx <- 3        # Instrumental function family to use
inst.func.family.c.idx <- 2      # Instrumental function family for continuous covariates
degree <- 3                      # Degree of spline function
gs.method <- 1                   # Use binary search
cov.ranges <- NULL    

# Set fixed parameters
idx.param.of.interest <- 2
H0 <- function(t) {t}
beta1 <- 1
beta2 <- -1
beta3 <- 2
beta4 <- -1
H0.inv <- function(t) {t}

# Verbose
verbose <- 3

#### Set and define variables ####

# Re-code the link.function vector
int2link.function <- function(int) {
  if (int == 1) {return("AFT_ll")}
  if (int == 2) {return("Cox_wb")}
}
link.function <- unname(int2link.function(link.function.idx))

# Re-code the search.method vector
int2search.method <- function(int) {
  if (int == 1) {return("GS")}
  if (int == 2) {return("EAM")}
}
search.method <- unname(int2search.method(search.method.idx))

# Re-code the instrumental function vector
int2inst.func.fam <- function(int) {
  if (int == 0) {return("ignore")}
  if (int == 1) {return("box")}
  if (int == 2) {return("spline")}
  if (int == 3) {return("cd")}
  if (int == 4) {return("cd.manycov")}
}
inst.func.family <- unname(int2inst.func.fam(inst.func.family.idx))
inst.func.family.c <- unname(int2inst.func.fam(inst.func.family.c.idx))

# Projection vector that projects the parameter vector onto the element of
# interest
c <- rep(0, n.cov + length(t))
c[idx.param.of.interest + length(t) - 1] <- 1

# Family of instrumental functions to be used when G.cd is selected
G.c <- NULL
if (inst.func.family %in% c("cd", "cd.manycov")) {
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
  next.gs.point = gs.binary
} else if (gs.method == 2) {
  next.gs.point = gs.regular
}

# Flag for recording the time that each step of the algorithm takes
time.run.duration <- TRUE

# Parameters related to data
beta.true <- function(t) {c(H0(t), beta1, beta2, beta3, beta4)}
par.space.lb <- rep(-10, n.cov + 1)
par.space.ub <- rep(10, n.cov + 1)
par.space <- matrix(c(par.space.lb, par.space.ub), ncol = 2)


#### Do adapted test of Bei (2024) ####

# Augment the parameter space, if necessary
if (length(t) > 1) {
  par.space.aug <- matrix(0, nrow = nrow(par.space) + length(t) - 1, ncol = ncol(par.space))
  par.space.aug[-(2:length(t)), ] <- par.space
  par.space.aug[2:length(t), ] <- matrix(rep(c(0, 10), length(t) - 1), ncol = 2, byrow = TRUE)
}
par.space <- par.space.aug

# Define hyperparameter list
options <- list(n.if.per.cov = n.if.per.cov,
                K.bar = K.bar,
                B = B,
                gs.method = gs.method,
                next.gs.point = next.gs.point,
                alpha = alpha,
                link.function = link.function,
                inst.func.family = inst.func.family,
                degree = degree,
                DGP = DGP,
                G.c = G.c,
                cov.ranges = cov.ranges)

# Generate the data
set.seed(seed)
data <- generateData_simMain(beta.true, n, n.cov, options, H0.inv,
                             plot.data = (verbose >= 4))

# Define the list of hyperparameters
hp <- set.hyperparameters(data, par.space, c, search.method, options)

# Run the test on a grid of values
inst.func.evals <- NULL
parallel <- FALSE
alpha <- 0.95
verbose <- TRUE
r.grid <- seq(par.space[which(c == 1), 1], par.space[which(c == 1), 2],
              length.out = 21)
results <- matrix(nrow = 0, ncol = 3)
for (r in r.grid) {
  
  # Run the test
  res <- test.point_Bei_MT(r = r,
                           c = c,
                           t = t,
                           par.space = par.space,
                           data = data,
                           hp = hp,
                           verbose = verbose,
                           inst.func.evals = inst.func.evals,
                           alpha = alpha,
                           parallel = parallel)
  
  # Store the results
  results <- rbind(results, res)
}

#### Store the results ####

file.name <- paste0(paste(args, collapse = "__"), ".csv")
path.name <- sprintf("%s/%s", master.dir, file.name)
write.csv(results, path.name, row.names = FALSE)



