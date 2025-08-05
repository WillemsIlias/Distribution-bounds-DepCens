
#### Preamble ####

on.VSC <- grepl("vsc36272", getwd())

# Load dependencies
if (on.VSC) {
  
  # Load dependencies
  Sys.setenv(TZ = "Europe/Amsterdam")
  lib.loc <- paste0(dirname(getwd()), "/R/")
  pkg.list <- c("MASS", "lubridate", "SPOT", "EnvStats", "splines2", "copula",
                "gplm", "ks", "PLRModels", "np", "tidyr", "tzdb", "readr",
                "forcats", "tidyverse", "survival", "stringr", "doParallel",
                "nloptr")
  for (pkg in pkg.list) {
    eval(bquote(require(.(pkg), lib.loc = lib.loc)))
  }
  source("simulationFunctions.R")
  
  # Set master directory
  master.dir <- "Simulations_Main_MT"
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
  master.dir <- "Simulations main MT"
  check_create.dir(master.dir)
  
  # Arguments supplied to the function (for testing purposes)
  args <- c(773,1000,100,31,2,1)
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
                link.function = link.function,
                DGP = DGP,
                inst.func.family = inst.func.family,
                degree = degree,
                G.c = G.c,
                cov.ranges = cov.ranges)

# Flag for recording the time that each step of the algorithm takes
time.run.duration <- TRUE

# Parameters related to data
beta.true <- function(t) {c(H0(t), beta1, beta2, beta3, beta4)}
par.space.lb <- rep(-10, n.cov + 1)
par.space.ub <- rep(10, n.cov + 1)
par.space <- matrix(c(par.space.lb, par.space.ub), ncol = 2)

#### Generate the data ####

# Set random seed of this iteration
set.seed(seed)

# Generate data 
data <- generateData_simMain(beta.true, n, n.cov, options, H0.inv,
                             plot.data = (verbose >= 4))

#### Estimate the model using only one time point ####

# Set level of the test
alpha <- 0.95
options$alpha <- alpha

# Find the identified set
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

#### Estimate the model using multiple time points: intersection ####

if (t.eval == 1) {
  t.to.check.vec <- seq(0, 1, length.out = 4)[2:4]
} else if (t.eval == 5) {
  t.to.check.vec <- seq(1, 5, length.out = 3)
}
alpha <- 1 - 0.05/3
options$alpha <- alpha

ident.set.MT <- list()
total.run.time.MT <- list()
total.leg.times.MT <- list()
conv.l.MT <- list()
conv.u.MT <- list()

for (t.idx in 1:length(t.to.check.vec)) {
  
  # Select t of this iteration
  t.to.check <- t.to.check.vec[t.idx]
  
  # Find the identified set
  fis.out <- find.identified.set(c, t.to.check, par.space, data, search.method,
                                 options, verbose, time.run.duration)
  ident.set.MT[[t.idx]] <- fis.out$ident.set
  
  # Analyse the running times
  chrono1 <- fis.out$chronometer1
  chrono2 <- fis.out$chronometer2
  total.run.time.MT[[t.idx]] <- chrono1$get.total.time(force = TRUE) + chrono2$get.total.time(force = TRUE)
  if (search.method == "EAM") {
    total.leg.times.MT[[t.idx]] <- chrono1$accumulate.legs(force = TRUE)[1:3] +
      chrono2$accumulate.legs(force = TRUE)[1:3]
  } else {
    total.leg.times.MT[t.idx] <- list(NULL)
  }
  
  # Convergence of the algorithm
  conv.l.MT[[t.idx]] <- fis.out$converge2
  conv.u.MT[[t.idx]] <- fis.out$converge1
  
}

#### Estimate the model using multiple time points: majority vote ####

if (t.eval == 1) {
  t.to.check.vec <- seq(0, 1, length.out = 4)[2:4]
} else if (t.eval == 5) {
  t.to.check.vec <- seq(1, 5, length.out = 3)
}
alpha <- 1 - 0.05/2
options$alpha <- alpha

ident.set.MV <- list()
total.run.time.MV <- list()
total.leg.times.MV <- list()
conv.l.MV <- list()
conv.u.MV <- list()

for (t.idx in 1:length(t.to.check.vec)) {
  
  # Select t of this iteration
  t.to.check <- t.to.check.vec[t.idx]
  
  # Find the identified set
  fis.out <- find.identified.set(c, t.to.check, par.space, data, search.method,
                                 options, verbose, time.run.duration)
  ident.set.MV[[t.idx]] <- fis.out$ident.set
  
  # Analyse the running times
  chrono1 <- fis.out$chronometer1
  chrono2 <- fis.out$chronometer2
  total.run.time.MV[[t.idx]] <- chrono1$get.total.time(force = TRUE) + chrono2$get.total.time(force = TRUE)
  if (search.method == "EAM") {
    total.leg.times.MV[[t.idx]] <- chrono1$accumulate.legs(force = TRUE)[1:3] +
      chrono2$accumulate.legs(force = TRUE)[1:3]
  } else {
    total.leg.times.MV[t.idx] <- list(NULL)
  }
  
  # Convergence of the algorithm
  conv.l.MV[[t.idx]] <- fis.out$converge2
  conv.u.MV[[t.idx]] <- fis.out$converge1
  
}

#### Store all results ####

# Initialize object that will store the results
ncol.sim.results <- 5*7 + 3 * 7 * as.numeric(search.method == "EAM") + 2
sim.results <- matrix(nrow = 0, ncol = ncol.sim.results)
cols <- c("ident.set.l", "ident.set.u", "conv.l", "conv.u", "total.run.time")
if (search.method == "EAM") {
  cols <- c(cols, "E-step", "A-step", "M-step")
}
cols.base <- cols
for (t.idx in 1:3) {
  cols <- c(cols, paste0(cols.base, sprintf("_int.%d", t.idx)))
}
for (t.idx in 1:3) {
  cols <- c(cols, paste0(cols.base, sprintf("_mv.%d", t.idx)))
}
cols <- c(cols, "seed", "per.cens")
colnames(sim.results) <- cols

# Censoring percentage in the data
per.cens <- 1 - sum(data$Delta)/n

# Store the results
sim.result <- c(ident.set, conv.l, conv.u, total.run.time, total.leg.times,
                ident.set.MT[[1]], conv.l.MT[[1]], conv.u.MT[[1]],
                total.run.time.MT[[1]], total.leg.times.MT[[1]],
                ident.set.MT[[2]], conv.l.MT[[2]], conv.u.MT[[2]],
                total.run.time.MT[[2]], total.leg.times.MT[[2]],
                ident.set.MT[[3]], conv.l.MT[[3]], conv.u.MT[[3]],
                total.run.time.MT[[3]], total.leg.times.MT[[3]],
                ident.set.MV[[1]], conv.l.MV[[1]], conv.u.MV[[1]],
                total.run.time.MV[[1]], total.leg.times.MV[[1]],
                ident.set.MV[[2]], conv.l.MV[[2]], conv.u.MV[[2]],
                total.run.time.MV[[2]], total.leg.times.MV[[2]],
                ident.set.MV[[3]], conv.l.MV[[3]], conv.u.MV[[3]],
                total.run.time.MV[[3]], total.leg.times.MV[[3]],
                seed, per.cens)
sim.results <- rbind(sim.results, sim.result)

#### Save results ####

comb <- list("n" = n, "n.cov" = n.cov, "link.function" = link.function, "DGP" = DGP,
             "search.method" = search.method, "alpha" = alpha, "gs.method" = gs.method,
             "K.bar" = K.bar, "n.if.per.cov" = n.if.per.cov,
             "inst.func.family" = inst.func.family, "degree" = degree,
             "inst.func.family.c" = inst.func.family.c, "B" = B, "t.eval" = t.eval)
defaults <- list(Kbar = 3,
                 Alpha = Inf,
                 IF = "cd",
                 method = 1,
                 Search = "GS",
                 Gc = "spline",
                 degree = 3)
results.idx <- unique(((seed - 4000123) %/% 5) %% 20)
results.dir <- sprintf("%s/res%s", master.dir, results.idx)
check_create.dir(results.dir)
path <- get.file.name(comb, 1, seed, results.dir, shortened = TRUE,
                      defaults = defaults)
write.csv(sim.results, path, row.names = FALSE)




