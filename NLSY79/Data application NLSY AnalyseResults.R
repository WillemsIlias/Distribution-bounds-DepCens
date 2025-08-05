
# Clear workspace
rm(list = ls())

# Load dependencies
library(stringr)

# Set results directory
results.dir <- "Results"

# Read the results
results <- data.frame(model = character(), stratum = numeric(),
                      lb.age = numeric(), ub.age = numeric(),
                      lb.educ = numeric(), ub.educ = numeric())
for (subdir in list.dirs(results.dir, recursive = FALSE)) {
  
  # Parse subdirectory name
  stratum <- as.numeric(substr(subdir, nchar(subdir), nchar(subdir)))
  
  for (model in c("AFT", "COX")) {
    
    # Get names of all files pertaining to the selected model
    model.files <- list.files(subdir)[grepl(model, list.files(subdir))]
    
    for (file.name in model.files) {
      
      # Read file
      file <- read.csv(paste0(subdir, "/", file.name))
      
      # Parse file name
      components <- str_split(gsub(".csv", "", file.name), "_")[[1]]
      model <- components[1]
      
      if (components[3] == "ageAtStart") {
        lb.age <- as.numeric(file$lower)
        ub.age <- as.numeric(file$upper)
      } else if (components[3] == "education") {
        lb.educ <- as.numeric(file$lower)
        ub.educ <- as.numeric(file$upper)
      }
    }
    
    # Add result to results data frame
    results <- rbind(results, c(model = model, stratum = stratum, lb.age = lb.age,
                     ub.age = ub.age, lb.educ = lb.educ, ub.educ = ub.educ))
  }
}

# Give proper column names
colnames(results) <- c("model", "stratum", "lb.age", "ub.age", "lb.educ", "ub.educ")

#### Make LaTeX tabular of the results ####

# AFT model results
results.AFT <- round(apply(results[results$model == "AFT", 2:6], 2, as.numeric), 2)

results.Cox <- round(apply(results[results$model == "COX", 2:6], 2, as.numeric), 2)

results.all <- cbind(results.AFT, results.Cox[, -1])

make.intervals <- function(row) {
  sprintf("%s & [%0.2f, %0.2f] & [%0.2f, %0.2f] & [%0.2f, %0.2f] & [%0.2f, %0.2f]", row[1], row[2], row[3],
          row[4], row[5], row[6], row[7], row[8], row[9])
}

results.matrix <- paste(apply(results.all, 1, make.intervals), collapse = "\\\\\n")
message(results.matrix)





  