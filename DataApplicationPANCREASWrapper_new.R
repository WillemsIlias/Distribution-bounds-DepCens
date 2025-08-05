
################################################################################

on.VSC <- grepl("vsc36272", getwd())
makeDesigns <- FALSE

################################################################################

#### Make matrix of combinations to check ####

if (!on.VSC & makeDesigns) {
  
  # Define settings
  t <- c(6, 12, 18)
  criteria.idx <- 1:3
  n.if.per.cov <- 5
  method <- 1:10
  model.idx <- 1:2
  
  # Construct matrix
  settings <- expand.grid(t, criteria.idx, n.if.per.cov, method, model.idx)
  colnames(settings) <- c("t", "crit_idx", "n_if_per_cov", "method", "model_idx")
  settings[settings$method != 1, "t"] <- 6
  settings <- unique(settings)
  write.csv(settings, "DataApplicationPANCREAS_settings.csv", row.names = FALSE)
}

#### Run method ####

if (on.VSC) {
  
  # Load dependencies
  Sys.setenv(TZ = "Europe/Amsterdam")
  lib.loc <- paste0(dirname(dirname(getwd())), "/R/")
  pkg.list <- c("MASS", "lubridate", "SPOT", "EnvStats", "splines2", "copula",
                "gplm", "ks", "PLRModels", "np", "doParallel", "nloptr")
  for (pkg in pkg.list) {
    eval(bquote(require(.(pkg), lib.loc = lib.loc)))
  }
  main.func.dir.name <- dirname(getwd())
  source(paste0(main.func.dir.name, "/simulationFunctions.R"))
  source(paste0(main.func.dir.name, "/Data application functions.R"))
  
  # Get arguments of function
  args <- as.numeric(commandArgs(TRUE))
  
  # Specify necessary paths
  path.data <- "PANCREAS/Pancreas_preprocessed.csv"
  path.crit.list <- "PANCREAS/criteria.list3.pancreas.RData"
  
} else {
  
  # Clear workspace
  rm(list = ls())
  on.VSC <- grepl("vsc36272", getwd())
  
  # Load dependencies
  pkg.list <- c("MASS", "lubridate", "SPOT", "EnvStats", "splines2", "copula",
                "gplm", "ks", "PLRModels", "np")
  for (pkg in pkg.list) {
    eval(bquote(require(.(pkg))))
  }
  source("simulationFunctions.R")
  source("Data application Functions.R")
  
  # Get arguments of function
  args <- c(6, 1, 5, 5, 2)
  path.data <- "Data application/PANCREAS/Pancreas_preprocessed.csv"
  path.crit.list <- "Data application/PANCREAS/criteria.list3.pancreas.RData"
  
}

# Extract arguments
t <- as.numeric(args[1])
criteria.idx <- as.numeric(args[2])
n.if.per.cov <- as.numeric(args[3])
method <- as.numeric(args[4])
model.idx <- as.numeric(args[5])

#### Get the data ####

# Load data set
data <- read.csv(path.data)

# Only keep rows with all non-missing entries
data <- na.omit(data)

# Add intercept column
data$X0 <- rep(1, nrow(data))

# Get subsetting criteria list
load(path.crit.list)

# Subset the data
criteria <- criteria.list3[[criteria.idx]]
sd.out <- subset.data(data, criteria, keep.covs = c("age_at_diag", "tumor_size"))
data.sub <- sd.out[[1]]
name.dict <- sd.out[[2]]

# Remove outlying data points
data.sub <- data.sub[which(data.sub$X1 > quantile(data.sub$X1, 0.025)),]
data.sub <- data.sub[which(data.sub$X1 < quantile(data.sub$X1, 0.975)),]
data.sub <- data.sub[which(data.sub$X2 > quantile(data.sub$X2, 0.025)),]
data.sub <- data.sub[which(data.sub$X2 < quantile(data.sub$X2, 0.975)),]

# Standardize the covariates
data.sub[, "X1"] <- (data.sub[, "X1"] - mean(data.sub[, "X1"]))/sd(data.sub[, "X1"])
data.sub[, "X2"] <- (data.sub[, "X2"] - mean(data.sub[, "X2"]))/sd(data.sub[, "X2"])

#### Run the algorithm and save the results ####

# Some general parameters
idxs.c <- c(1, 2)
n.param <- sum(grepl("X[[:digit:]]+$", colnames(data.sub)))
add.options <- list(n.if.per.cov = n.if.per.cov)
par.space <- matrix(c(rep(-10, n.param), rep(10, n.param)), ncol = 2)

# Set correct model
model.names <- c("AFT_ll", "Cox_wb")
model.name <- model.names[model.idx]
add.options$link.function <- model.name

# Set-up for selected method
typeIerror <- 0.05
if (method == 1) {          # Single t 
  t.to.check <- t
  alpha.seq <- 1 - typeIerror
} else if (method == 2) {   # Three different t, intersection, equal \alpha
  t.to.check <- c(6, 12, 18)
  alpha.seq <- 1 - typeIerror/3
} else if (method == 3) {   # Three different t, intersection, unequal \alpha
  t.to.check <- c(6, 12, 18)
  alpha.seq <- 1 - ((typeIerror/sum(1/(1:3))) * (1/1:3))
} else if (method == 4) {   # Three different t, majority vote
  t.to.check <- c(6, 12, 18)
  alpha.seq <- 1 - typeIerror/2
} else if (method == 5) {   # Five different t, intersection, equal \alpha
  t.to.check <- c(2, 6, 10, 14, 18)
  alpha.seq <- 1 - typeIerror/5
} else if (method == 6) {   # Five different t, intersection, unequal \alpha
  t.to.check <- c(2, 6, 10, 14, 18)
  alpha.seq <- 1 - ((typeIerror/sum(1/(1:5))) * (1/1:5))
} else if (method == 7) {   # Five different t, majority vote
  t.to.check <- c(2, 6, 10, 14, 18)
  alpha.seq <- 1 - 1/2
} else if (method == 8) {   # Nine different t, intersection, equal \alpha
  t.to.check <- c(2, 4, 6, 8, 10, 12, 14, 16, 18)
  alpha.seq <- 1 - typeIerror/9
} else if (method == 9) {   # Nine different t, intersection, unequal \alpha
  t.to.check <- c(2, 4, 6, 8, 10, 12, 14, 16, 18)
  alpha.seq <- 1 - ((typeIerror/sum(1/(1:9))) * (1/1:9))
} else if (method == 10) {   # Nine different t, majority vote
  t.to.check <- c(2, 4, 6, 8, 10, 12, 14, 16, 18)
  alpha.seq <- 1 - 1/2
}

# Initialize object that will store the results
results.list <- list()

# Run the algorithm for each t, at the correct level.
for (t.idx in 1:length(t.to.check)) {
  
  # Select the time point of this iteration
  t <- t.to.check[t.idx]
  
  # Set correct level
  if (length(alpha.seq) == 1) {
    alpha <- alpha.seq
  } else {
    alpha <- alpha.seq[t.idx]
  }
  add.options$alpha <- alpha
  
  # age_at_diag
  idx.param.of.interest <- 2
  results1 <- pi.surv(data.sub, idx.param.of.interest, idxs.c, t, par.space,
                      add.options = add.options)
  
  # tumor_size
  idx.param.of.interest <- 3
  results2 <- pi.surv(data.sub, idx.param.of.interest, idxs.c, t, par.space,
                      add.options = add.options)
  
  # Add results to results list
  results <- rbind(results1, results2)
  results.list[[sprintf("t%s_alpha%s", t, round(alpha, 4))]] <- results
}

# Combine the results. Store in results list
if (method %in% c(2, 3, 5, 6, 8, 9)) {
  results.combined <- cbMV(results.list, 1)
} else if (method %in% c(4, 7, 10)) {
  results.combined <- cbMV(results.list, 0.5)
} else {
  results.combined <- NULL
}
results.list[["combined"]] <- results.combined

# Save the results
file.name <- sprintf("link-%s__t-%d__crit-%d__nifpc-%d__method-%d.csv",
                     model.name, t, criteria.idx, n.if.per.cov, method)
write.csv(results.list, file = file.name, row.names = FALSE)





