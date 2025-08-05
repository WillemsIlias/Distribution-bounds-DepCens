
#### Refined results ####

# Clear workspace
rm(list = ls())

# Load dependencies
library(dplyr)

# Set general parameters
master.dir <- "Simulations_Main_trueBounds"

# Create data frame that will store the results
res.df <- data.frame(DGP = numeric(),
                     beta1.lb = numeric(),
                     beta1.ub = numeric(),
                     beta2.lb = numeric(),
                     beta2.ub = numeric(),
                     finished = logical())

# True values
beta0 <- 1
beta1 <- 1
beta2 <- -1

# Read and process all the refined results
DGP_min <- 21
DGP_max <- 32
for (DGP in DGP_min:DGP_max) {
  for (run in 0:4) {
    
    # DGP + seed of the simulation results of this iteration
    DGP.run <- DGP + run*1000
    
    # Get file name of refined results
    file.name <- sprintf("%s/eval_array%d_refined.Rdata", master.dir, DGP.run)
    
    # Skip this iteration if file does not exist
    if (!file.exists(file.name)) {
      next
    }
    
    # Load the results
    load(file.name)
    
    # Compute identified intervals for beta1 and beta2
    beta1.lb <- min(min(feas.beta[, 2]), beta1)
    beta1.ub <- max(max(feas.beta[, 2]), beta1)
    beta2.lb <- min(min(feas.beta[, 3]), beta2)
    beta2.ub <- max(max(feas.beta[, 3]), beta2)
    
    # Determine if this simulation finished running
    finished <- by.val < 0.05
    
    # Store all results
    res.df <- rbind(res.df,
                    c(DGP, beta1.lb, beta1.ub, beta2.lb, beta2.ub,finished))
  }
}
colnames(res.df) <- c("DGP", "beta1.lb", "beta1.ub", "beta2.lb", "beta2.ub",
                      "finished")

# Compute averages
res.df %>% group_by(DGP) %>% summarise(mean.beta1.lb = mean(beta1.lb),
                                       mean.beta1.ub = mean(beta1.ub),
                                       mean.beta2.lb = mean(beta2.lb),
                                       mean.beta2.ub = mean(beta2.ub))

#### Crude results ####

# Clear workspace
rm(list = ls())

# Load dependencies
# ...

# Set general parameters
master.dir <- "Simulations_Main_trueBounds"

# Create data frame that will store the results
res.df <- data.frame(DGP = numeric(),
                     beta1.lb = numeric(),
                     beta1.ub = numeric(),
                     beta2.lb = numeric(),
                     beta2.ub = numeric())

# Read and process all the refined results
DGP_min <- 21
DGP_max <- 32

# Read and process all the crude results
for (DGP in DGP_min:DGP_max) {
  
  # Get file name of refined results
  file.name <- sprintf("%s/eval_array%d_intermediate.Rdata", master.dir, DGP)
  
  # Load the results
  load(file.name)
  
  # Set true values for beta
  beta1 <- 1
  beta2 <- -1
  
  # Compute set of feasible betas
  tol <- 1e-15
  by.val <- 0.01
  feas.beta.idxs <- arrayInd(which(eval.array <= tol), .dim = dim(eval.array))
  if (nrow(feas.beta.idxs) == 0) {
    stop(sprintf("No initial feasible point found (DGP = %d).", DGP))
  }
  beta0.grid <- seq(-10, 10, by = by.val)
  beta1.grid <- setdiff(seq(-1, 1, by = 5*by.val), 0) + beta1
  beta2.grid <- setdiff(seq(-1, 1, by = 5*by.val), 0) + beta2
  feas.beta <- t(matrix(apply(feas.beta.idxs, 1, function(row) {c(beta0.grid[row[1]], beta1.grid[row[2]], beta2.grid[row[3]])}),
                        ncol = nrow(feas.beta.idxs)))

  
  # Compute identified intervals for beta1 and beta2
  beta1.lb <- min(feas.beta[, 2])
  beta1.ub <- max(feas.beta[, 2])
  beta2.lb <- min(feas.beta[, 3])
  beta2.ub <- max(feas.beta[, 3])
  
  # Store all results
  res.df <- rbind(res.df, c(DGP, beta1.lb, beta1.ub, beta2.lb, beta2.ub))
}
colnames(res.df) <- c("DGP", "beta1.lb", "beta1.ub", "beta2.lb", "beta2.ub")
res.df

