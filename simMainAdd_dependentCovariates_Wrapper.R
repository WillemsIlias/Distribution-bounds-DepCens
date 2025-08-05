
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
  master.dir <- "Simulations_Main_depCov"
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
  source("Data application functions.R")
  
  # Set master directory
  master.dir <- "Simulations main"
  check_create.dir(master.dir)
  
  # Arguments supplied to the function (for testing purposes)
  args <- c(3038663,2000,8,27,2,1)
}

#### Extract and set parameters ####

# Extract arguments for simulation design
seed <- as.numeric(args[1])
n <- as.numeric(args[2])
n.if.per.cov <- as.numeric(args[3])
DGP <- as.numeric(args[4])
link.function.idx <- as.numeric(args[5])
t.eval <- as.numeric(args[6])

# Set fixed settings
K.bar <- 3                       # Numer of critical value refinements.
# Recommendation in Bei (2024) (somewhere between 2 and 10)
n.cov <- 2                       # Number of covariates
alpha <- 0.95                    # 1 - Level of the test
search.method.idx <- 1           # Do not use EAM algorithm.
B <- 600                         # Number of bootstraps. Recommendation in Bei (2024)
inst.func.family.idx <- 3        # Overall instrumental function family
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

#### Some preliminary steps ####

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
}
inst.func.family <- unname(int2inst.func.fam(inst.func.family.idx))
inst.func.family.c <- unname(int2inst.func.fam(inst.func.family.c.idx))

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
  next.gs.point = gs.binary
} else if (gs.method == 2) {
  next.gs.point = gs.regular
}

# List of hyperparameters to be supplied to various functions later on
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
                cov.ranges = cov.ranges,
                ignore.empty.IF = TRUE)

# Flag for recording the time that each step of the algorithm takes
time.run.duration <- TRUE

# Parameters related to data
beta.true <- function(t) {c(H0(t), beta1, beta2, beta3, beta4)}
par.space.lb <- rep(-10, n.cov + 1)
par.space.ub <- rep(10, n.cov + 1)
par.space <- matrix(c(par.space.lb, par.space.ub), ncol = 2)

#### Perform the simulations ####

# Initialize object that will store the results
ncol.sim.results <- 8 + 3 * as.numeric(search.method == "EAM")
sim.results <- matrix(nrow = 0, ncol = ncol.sim.results)
cols <- c("ident.set.l", "ident.set.u", "conv.l", "conv.u", "total.run.time")
if (search.method == "EAM") {
  cols <- c(cols, "E-step", "A-step", "M-step")
}
cols <- c(cols, "seed", "per.cens", "coef.ident.model")
colnames(sim.results) <- cols

# Set random seed of this iteration
set.seed(seed)

# Generate data 
data <- generateData_simAdd(beta.true, n, n.cov, options, H0.inv,
                            plot.data = (verbose >= 4))

# Find the identified set
fis.out <- find.identified.set(c, t.eval, par.space, data, search.method,
                               options, verbose = verbose,
                               time.run.duration = time.run.duration)
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

#### Estimate an identified model assuming independence (for reference) ####

coef.identified.model <- get.coef.identified.model(c, t.eval, data, link.function)

#### Store all results ####

# Censoring percentage in the data
per.cens <- 1 - sum(data$Delta)/n

# Store the results
sim.result <- c(ident.set, conv.l, conv.u, total.run.time, total.leg.times,
                seed, per.cens, coef.identified.model)
sim.results <- rbind(sim.results, sim.result)

#### Save results ####

comb <- list("n" = n, "n.cov" = n.cov, "link.function" = link.function, "DGP" = DGP,
             "search.method" = search.method, "alpha" = alpha, "gs.method" = gs.method,
             "K.bar" = K.bar, "n.if.per.cov" = n.if.per.cov,
             "inst.func.family" = inst.func.family, "degree" = degree,
             "inst.func.family.c" = inst.func.family.c, "B" = B, "t.eval" = t.eval)
defaults <- list(Kbar = 3,
                 Alpha = 0.95,
                 IF = "cd",
                 method = 1,
                 Search = "GS",
                 Gc = "spline",
                 degree = 3)
results.idx <- unique(((seed - 3000123) %/% 5) %% 100)
results.dir <- sprintf("%s/res%s", master.dir, results.idx)
check_create.dir(results.dir)
path <- get.file.name(comb, 1, seed, results.dir, shortened = TRUE,
                      defaults = defaults)
write.csv(sim.results, path, row.names = FALSE)