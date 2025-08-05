
# Clear workspace
rm(list = ls())

# Load dependencies
# ...

# Set main directory
master.dir <- "Simulations_Add_timeIndep"

# Read the results
results.df <- data.frame(seed = numeric(), n = numeric(), nifpc = numeric(),
                         DGP = numeric(), link.idx = numeric(),
                         t.vec.idx = numeric(), lb = numeric(), ub = numeric())
for (file.name in list.files(master.dir)) {
  
  # Make full file name
  full.file.name <- sprintf("%s/%s", master.dir, file.name)
  
  # Read the file
  file <- read.csv(full.file.name)
  
  # Parse file name
  comps <- as.list(as.numeric(strsplit(gsub(".csv", "", file.name), "__")[[1]]))
  names(comps) <- c("seed", "n", "nifpc", "DGP", "link.idx", "t.vec.idx")
  
  # Compute the identified interval
  lb <- min(file[file$t.stat <= file$crit.val, "theta"])
  ub <- max(file[file$t.stat <= file$crit.val, "theta"])
  
  # Store the results
  row <- c(comps, "lb" = lb, "ub" = ub)
  results.df <- rbind(results.df, row)
}

# List of values for t used.
t.vec.list <- list(
  c(1, 3, 5),
  c(0.3333, 0.6666, 1),
  c(2, 4, 6, 8, 10)
)
