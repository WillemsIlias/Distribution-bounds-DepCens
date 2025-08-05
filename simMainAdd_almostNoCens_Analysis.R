
# Clear workspace
rm(list = ls())

# Load dependencies
source("simulationAnalysisFunctions.R")
source("lowLevelFunctions.R")
source("simulationAnalysisFunctions.R")

#### Set all variables ####

# General settings
idx.param.of.interest <- 2
master.dir <- "Simulations_Main_almostNoCens"

# Set true parameter vector
H0 <- function(t) {t}
H0.inv <- function(t) {t}
beta1 <- 1
beta2 <- -1
beta3 <- 2
beta4 <- -1
beta.true <- function(t) {c(H0(t), beta1, beta2, beta3, beta4)}

#### If necessary, extract results from subfiles ####

check_create.dir(sprintf("%s/%s", master.dir, "LF-AFT_ll"))
check_create.dir(sprintf("%s/%s", master.dir, "LF-Cox_wb"))

dirs.main.dir <- list.dirs(master.dir, recursive = FALSE)
res.subdirs <- dirs.main.dir[grepl("res", dirs.main.dir)]
for (dir.name in res.subdirs) {
  
  # For each folder in the subdirectory...
  for (folder in list.dirs(dir.name, recursive = FALSE, full.names = FALSE)) {
    
    # Construct full folder path
    full.folder.path <- sprintf("%s/%s", dir.name, folder)
    
    if (length(list.files(full.folder.path)) != 180) {
      warning(sprintf("%s: not enough simulations!", full.folder.path))
    }
    
    # Move all files from the folder to the equally named folder in master 
    # directory.
    for (file in list.files(full.folder.path)) {
      
      # Get full file path
      full.file.path <- sprintf("%s/%s", full.folder.path, file)
      
      # Get new file path
      new.file.path <- sprintf("%s/%s/%s", master.dir, folder, file)
      
      # Move the file
      file.rename(full.file.path, new.file.path)
    }
  }
}

#### Analyze the results ####

# Get a list of all simulation designs
shortened <- TRUE
defaults <- list(Kbar = 3,
                 Alpha = 0.95,
                 IF = "cd",
                 Search = "GS",
                 method = 1,
                 Gc = "spline",
                 degree = 3)
sim.settings <- get.sim.settings(master.dir, defaults)

# For each design, read the simulations
sim.results <- NULL
for (setting.idx in 1:nrow(sim.settings)) {
  print(setting.idx)
  setting <- sim.settings[setting.idx, ]
  sim.files <- get.sim.files(master.dir, setting, shortened, defaults)
  summary <- summarize.files(sim.files, beta.true(setting[, "t.eval"]),
                             idx.param.of.interest, setting.idx = setting.idx)
  if (is.na(setting[, "nifpc"])) {
    setting[, "nifpc"] <- setting[, "nbpc"]
  }
  setting <- setting[, which(colnames(setting) != "nbpc")]
  sim.params <- setting[, c("DGP", "n", "n.cov", "IF", "nifpc", "degree", "Gc", "Alpha")]
  sim.results <- rbind(sim.results, cbind(sim.params, summary))
}
sim.results$DGP <- as.numeric(sim.results$DGP)
sim.results <- sim.results[order(sim.results$avg.width), ]

#### Produce a LaTeX table of the results ####

# Results for AFT model
sim.results.subset.AFT_LL <- sim.results[sim.results$DGP %in% 33:35, ]
results2latex(sim.results.subset.AFT_LL,
              "Almost no censoring. Results AFT model",
              "tab: X1 normal, almost no censoring, AFT")

# Results for Cox model
sim.results.subset.Cox_wb <- sim.results[sim.results$DGP %in% 36:38, ]
results2latex(sim.results.subset.Cox_wb,
              "Almost no censoring. Results Cox model",
              "tab: X1 normal, almost no censoring, Cox")





