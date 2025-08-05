# Clear workspace
rm(list = ls())

# Load dependencies
library(ggplot2)
library(gridExtra)
library(survival)
main.func.dir.name <<- paste0(dirname(dirname(getwd())))
source(paste0(main.func.dir.name, "/simulationFunctions.R"))
source(paste0(main.func.dir.name, "/Data application functions.R"))

# Set some parameters
master.dir <- "PANCREAS results new"
sim.results <- data.frame(link = character(), t = character(), crit = character(),
                          nifpc = character(), method = character(), var = numeric(),
                          lower = numeric(), upper = numeric(),
                          conv.l = numeric(), conv.u = numeric())

#### Read all results ####

cov.names <- c("Age", "Tumor size")
for (file.name in list.files(master.dir)) {
  
  # Get the content of the file
  file <- read.csv(sprintf("%s/%s", master.dir, file.name))
  
  # Parse file name
  comps <- strsplit(strsplit(gsub(".csv", "", file.name), "__")[[1]], "-")
  names(comps) <- lapply(comps, function(elem) {elem[1]})
  comps <- lapply(comps, function(elem) {elem[-1]})
  
  # Append results to file name info
  method.idxs.intersection <- c(2, 3, 5, 6, 8, 9)
  method.idxs.majvote <- c(4, 7, 10)
  for (row.idx in 1:nrow(file)) {
    if (comps[["method"]] == 1) {
      res <- as.list(as.numeric(file[row.idx, ]))
      names(res) <- c("lower", "upper", "conv.l", "conv.u")
      row <- c(comps, var = cov.names[row.idx], res)
      sim.results <- rbind(sim.results, row)
    } else if (as.numeric(comps[["method"]]) %in% method.idxs.intersection) {
      res <- list(lower = file[row.idx, "combined.lower"],
                  upper = file[row.idx, "combined.upper"],
                  conv.l = file[row.idx, "combined.conv.l"],
                  conv.u = file[row.idx, "combined.conv.u"])
      row <- c(comps, var = cov.names[row.idx], res)
      sim.results <- rbind(sim.results, row)
    } else if (comps[["method"]] %in% method.idxs.majvote) {
      
      # Obtain precomputed combined bounds
      res <- list(lower = file[row.idx, "combined.lower"],
                  upper = file[row.idx, "combined.upper"],
                  conv.l = file[row.idx, "combined.conv.l"],
                  conv.u = file[row.idx, "combined.conv.u"])
      
      # If bounds are Inf (model misspecification at one point), compute the
      # combined bounds manually.
      if (res$conv.l == 2 | res$conv.u == 2) {
        lbs <- unname(unlist(file[row.idx, grepl("0.5.lower", colnames(file))]))
        ubs <- unname(unlist(file[row.idx, grepl("0.5.upper", colnames(file))]))
        temp <- list()
        for (i in 1:length(lbs)) {
          temp[[i]] <- data.frame("lower" = lbs[i], "upper" = ubs[i],
                                  "conv.l" = 1, "conv.u" = 1)
        }
        res <- cbMV(temp, 0.5)
        res <- list(lower = res[, "lower"],
                    upper = res[, "upper"],
                    conv.l = res[, "conv.l"],
                    conv.u = res[, "conv.u"])
      }
      
      # Add to results
      row <- c(comps, var = cov.names[row.idx], res)
      sim.results <- rbind(sim.results, row)
    }
  }
}
sim.results[!(colnames(sim.results) %in% c("link", "var"))] <- apply(sim.results[!(colnames(sim.results) %in% c("link", "var"))], 2, as.numeric)
sim.results <- sim.results[order(sim.results[, "method"]), ]


#### Fit a parametric model assuming independence for reference ####

## Load the data

# Load data set
path.data <- "pancreas_preprocessed.csv"
data <- read.csv(path.data)

# Only keep rows with all non-missing entries
data <- na.omit(data)

# Add intercept column
data$X0 <- rep(1, nrow(data))

# Get subsetting criteria list
path.crit.list <- "criteria.list3.pancreas.Rdata"
load(path.crit.list)

## Select criterion index

model.params.list <- list()
for (lnk in c("Cox_wb", "AFT_ll")) {
  model.params <- matrix(nrow = 0, ncol = 4)
  colnames(model.params) <- c("age_at_diag", "tumor_size", "lambda", "alpha")
  for (criteria.idx in 1:3) {
    
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
    
    ## Parametric model
    
    # Fit the model
    data.sub$Y <- data.sub$Y + 0.1
    if (lnk == "AFT_ll") {
      fit <- survreg(Surv(Y, Delta) ~ X1 + X2, data = data.sub, dist = "loglogistic")
      alpha <- 1/fit$scale
      lambda <- exp(-(fit$coefficients[1]/fit$scale))
      
      # Store the coefficients (note that we need to take the negative, since the
      # coefficients of the model in log-linear representation are returned).
      model.params <- rbind(model.params, -fit$coefficients[2:3])
    } else if (lnk == "Cox_wb") {
      fit <- coxph(Surv(Y, Delta) ~ X1 + X2, data = data.sub)
      model.params <- rbind(model.params, c(fit$coefficients, NA, NA))
    }
  }
  model.params.list[[lnk]] <- model.params
}

#### Select model under investigation ####

link <- "AFT_ll"
model.params <- model.params.list[[link]]

#### Plot the results for single point test + 1 combination strategy ####

# Select appropriate subset of the data
sim.results.link <- sim.results[sim.results$link == link, ]

# Colors for variables
age.col <- "gray60"
ts.col <- "grey20"

# Make all plots
titles <- paste("", paste(c("local", "regional", "distant"), "cancer"))
plt.list <- list()
for (crit in 1:3) {
  
  # Select stratum criterion and subset data set.
  res.crit <- sim.results.link[which(sim.results.link$crit == crit), ]
  colnames(res.crit)[colnames(res.crit) == "var"] <- "Variable"
  
  # Further subset data to single time point analyses and one time-independent
  # strategy.
  res.crit <- res.crit[res.crit$method %in% c(1, 2), ]
  res.crit[res.crit$method != 1, "t"] <- 0
  
  # For jittering identified intervals of the coefficients horizontally.
  jitter.size <- 0.4
  res.crit$t <- ifelse(res.crit$Variable == "Age", res.crit$t - jitter.size,
                       res.crit$t + jitter.size)
  
  # Initialize plot
  plt <- ggplot(mapping = aes(x = t, colour = Variable), data = res.crit) +
    scale_color_manual(name = "Variable", 
                       values = c("Age" = age.col, "Tumor size" = ts.col))
  
  # Line at y = 0 and x = 0
  plt <- plt + geom_segment(aes(x = -10, y = 0, xend = 40, yend = 0), colour = "black")
  plt <- plt + geom_segment(aes(x = 0, y = -10, xend = 0, yend = 10), colour = "black")
  
  # Plot identified intervals
  plt <- plt + geom_segment(aes(y = lower, xend = t, yend = upper), size = 1)
  plt <- plt + coord_cartesian(xlim = c(-1, 22), ylim = c(-0.5, 1.5))
  
  # Plot the results from the fully parametric model
  plt <- plt + geom_hline(yintercept = model.params[crit, 1],
                          linetype = "dashed", color = age.col, size = 1)
  plt <- plt + geom_hline(yintercept = model.params[crit, 2],
                          linetype = "dashed", color = ts.col, size = 1)
  
  # Some layout changes
  plt <- plt + labs(y = "effect size",
                    x = "time",
                    title = titles[crit],
                    colour = "Variable")
  plt <- plt + theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 15))
  plt <- plt + theme(plot.title = element_text(size = 12))
  plt <- plt + theme(legend.position = "bottom")
  plt <- plt + scale_x_continuous(breaks = c(6, 12, 18))
  
  # If the plot is not the left-most one, remove y axis label
  if (crit > 1) {
    plt <- plt + labs(y = "")
  }
  
  # If the plot is not the right-most one, remove the legend
  if (crit < 3) {
    plt <- plt + theme(legend.position = "none")
  }
  
  print(plt)
  
  plt.list[[crit]] <- plt
}

# Extract legend
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}
mylegend <- g_legend(plt.list[[3]])

# Put plots side-by-side and save the figure
plt.full <- grid.arrange(arrangeGrob(plt.list[[1]], plt.list[[2]],
                                     plt.list[[3]] + theme(legend.position = "none"),
                                     nrow = 1),
                         mylegend, nrow = 2, heights = c(10, 1))
plt.full
ggsave(plot = plt.full,
       filename = sprintf("results_all_%s.png", link),
       device = "png")


#### Plot the results comparing all combination strategies (Cox) ####

link <- "Cox_wb"
model.params <- model.params.list[[link]]

# Select appropriate subset of the data
sim.results.link <- sim.results[sim.results$link == link, ]

# Make all plots
titles <- paste("", paste(c("local", "regional", "distant"), "cancer"))
plt.list <- list()
for (crit in 1:3) {
  
  # Select stratum criterion and subset data set.
  res.crit <- sim.results.link[which(sim.results.link$crit == crit), ]
  colnames(res.crit)[colnames(res.crit) == "var"] <- "Variable"
  
  # when model was misspecified, the identified interval should be empty.
  res.crit[res.crit$lower == -Inf, ] <- 0
  res.crit[res.crit$upper == Inf, ] <- 0
  
  # Further subset data to single time point analyses and one time-independent
  # strategy.
  res.crit <- res.crit[res.crit$method %in% 2:10, ]
  res.crit[, "t"] <- res.crit[, "method"]
  
  # Order the combination methods on the plot
  new.order <- c(1, 1, 5, 9, 2, 6, 10, 3, 7, 11)
  res.crit[, "t"] <- new.order[res.crit[, "t"]]
  
  # Add empty space for readibility
  
  # For jittering identified intervals of the coefficients horizontally.
  jitter.size <- 0.15
  res.crit$t <- ifelse(res.crit$Variable == "Age", res.crit$t - jitter.size,
                       res.crit$t + jitter.size)
  
  # Initialize plot
  plt <- ggplot(mapping = aes(x = t, colour = Variable), data = res.crit) +
    scale_color_manual(name = "Variable", 
                       values = c("Age" = age.col, "Tumor size" = ts.col))
  
  # Line at y = 0 and x = 0
  plt <- plt + geom_segment(aes(x = -10, y = 0, xend = 40, yend = 0), colour = "black")

  # Plot identified intervals
  plt <- plt + geom_segment(aes(y = lower, xend = t, yend = upper), size = 1)
  plt <- plt + coord_cartesian(xlim = c(0, 12), ylim = c(-0.5, 1.5))
  
  # Plot the results from the fully parametric model
  plt <- plt + geom_hline(yintercept = model.params[crit, 1],
                          linetype = "dashed", color = age.col, size = 1)
  plt <- plt + geom_hline(yintercept = model.params[crit, 2],
                          linetype = "dashed", color = ts.col, size = 1)
  
  # Some layout changes
  plt <- plt + labs(y = "effect size (Cox)",
                    x = NULL,
                    title = titles[crit],
                    colour = "Variable")
  plt <- plt + theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 15))
  plt <- plt + theme(plot.title = element_text(size = 12))
  plt <- plt + theme(legend.position = "bottom")
  plt <- plt + scale_x_continuous(breaks = c(1, 2, 3, 5, 6, 7, 9, 10, 11))
  
  # If the plot is not the left-most one, remove y axis label
  if (crit > 1) {
    plt <- plt + labs(y = "")
  }
  
  # If the plot is not the right-most one, remove the legend
  if (crit < 3) {
    plt <- plt + theme(legend.position = "none")
  }
  
  print(plt)
  
  plt.list[[crit]] <- plt
}
plt.list.cox <- plt.list

# Extract legend
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}
mylegend <- g_legend(plt.list[[3]])

# Put plots side-by-side and save the figure
plt.full.cox <- grid.arrange(arrangeGrob(plt.list[[1]], plt.list[[2]],
                                         plt.list[[3]] + theme(legend.position = "none"),
                                         nrow = 1),
                             mylegend, nrow = 2, heights = c(10, 1))
plt.full.cox

#### Plot the results comparing all combination strategies (AFT) ####

link <- "AFT_ll"
model.params <- model.params.list[[link]]

# Select appropriate subset of the data
sim.results.link <- sim.results[sim.results$link == link, ]

# Make all plots
titles <- paste("", paste(c("local", "regional", "distant"), "cancer"))
plt.list <- list()
for (crit in 1:3) {
  
  # Select stratum criterion and subset data set.
  res.crit <- sim.results.link[which(sim.results.link$crit == crit), ]
  colnames(res.crit)[colnames(res.crit) == "var"] <- "Variable"
  
  # when model was misspecified, the identified interval should be empty.
  res.crit[res.crit$lower == -Inf, ] <- 0
  res.crit[res.crit$upper == Inf, ] <- 0
  
  # Further subset data to single time point analyses and one time-independent
  # strategy.
  res.crit <- res.crit[res.crit$method %in% 2:10, ]
  res.crit[, "t"] <- res.crit[, "method"]
  
  # Order the combination methods on the plot
  new.order <- c(1, 1, 5, 9, 2, 6, 10, 3, 7, 11)
  res.crit[, "t"] <- new.order[res.crit[, "t"]]
  
  # Add empty space for readibility
  
  # For jittering identified intervals of the coefficients horizontally.
  jitter.size <- 0.15
  res.crit$t <- ifelse(res.crit$Variable == "Age", res.crit$t - jitter.size,
                       res.crit$t + jitter.size)
  
  # Initialize plot
  plt <- ggplot(mapping = aes(x = t, colour = Variable), data = res.crit) +
    scale_color_manual(name = "Variable", 
                       values = c("Age" = age.col, "Tumor size" = ts.col))
  
  # Line at y = 0 and x = 0
  plt <- plt + geom_segment(aes(x = -10, y = 0, xend = 40, yend = 0), colour = "black")
  
  # Plot identified intervals
  plt <- plt + geom_segment(aes(y = lower, xend = t, yend = upper), size = 1)
  plt <- plt + coord_cartesian(xlim = c(0, 12), ylim = c(-0.5, 1.5))
  
  # Plot the results from the fully parametric model
  plt <- plt + geom_hline(yintercept = model.params[crit, 1],
                          linetype = "dashed", color = age.col, size = 1)
  plt <- plt + geom_hline(yintercept = model.params[crit, 2],
                          linetype = "dashed", color = ts.col, size = 1)
  
  # Some layout changes
  plt <- plt + labs(y = "effect size (AFT)",
                    x = NULL,
                    title = NULL,
                    colour = "Variable")
  plt <- plt + theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 15))
  plt <- plt + theme(plot.title = element_text(size = 12))
  plt <- plt + theme(legend.position = "bottom")
  plt <- plt + scale_x_continuous(breaks = c(1, 2, 3, 5, 6, 7, 9, 10, 11))
  
  # If the plot is not the left-most one, remove y axis label
  if (crit > 1) {
    plt <- plt + labs(y = "")
  }
  
  # If the plot is not the right-most one, remove the legend
  if (crit < 3) {
    plt <- plt + theme(legend.position = "none")
  }
  
  print(plt)
  
  plt.list[[crit]] <- plt
}
plt.list.aft <- plt.list

# Extract legend
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}
mylegend <- g_legend(plt.list[[3]])

# Put plots side-by-side and save the figure
plt.full.aft <- grid.arrange(arrangeGrob(plt.list[[1]], plt.list[[2]],
                                         plt.list[[3]] + theme(legend.position = "none"),
                                         nrow = 1),
                             mylegend, nrow = 2, heights = c(10, 1))
plt.full.aft

#### Create full plot ####

# Extract legend
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}
mylegend <- g_legend(plt.list[[3]])

# Put plots side-by-side and save the figure
plt.full.all <- grid.arrange(arrangeGrob(plt.list.cox[[1]],
                                         plt.list.cox[[2]],
                                         plt.list.cox[[3]] + theme(legend.position = "none"),
                                         plt.list.aft[[1]],
                                         plt.list.aft[[2]],
                                         plt.list.aft[[3]] + theme(legend.position = "none"),
                                         nrow = 2, heights = c(11, 10)),
                             mylegend, nrow = 2, heights = c(10, 1))
plt.full.all

ggsave(plot = plt.full.all,
       filename = "results_comparison_all.png",
       device = "png")

# Some explanation on the combination methods:

# 1 - Three different time points, equal Bonferroni
# 2 - Five different time points, equal Bonferroni
# 3 - Nine different time points, equal Bonferroni

# 4 - Three different time points, unequal Bonferroni
# 5 - Five different time points, unequal Bonferroni
# 6 - Nine different time points, unequal Bonferroni

# 7 - Three different time points, majority vote
# 8 - Five different time points, majority vote
# 9 - Nine different time points, majority vote





#### Analyze the time indep. results using extended test ####

# Read the results
TI.dir <- "PANCREAS results TI fine"
sim.results.TI <- data.frame(link = character(), t = character(),
                             crit = character(), nifpc = character(),
                             method = character(), lb = numeric(),
                             ub = numeric())

for (file.name in list.files(TI.dir)) {
  
  # Get the content of the file
  file <- read.csv(sprintf("%s/%s", TI.dir, file.name))
  
  # Parse file name
  comps <- strsplit(strsplit(gsub(".csv", "", file.name), "__")[[1]], "-")
  names(comps) <- lapply(comps, function(elem) {elem[1]})
  comps <- lapply(comps, function(elem) {elem[-1]})
  
  # Obtain the lower and upper bounds
  lb <- min(file[file$t.stat <= file$crit.val, "theta"])
  ub <- max(file[file$t.stat <= file$crit.val, "theta"])
  
  # Append results to file name info
  row <- c(comps, lb = lb, ub = ub)
  sim.results.TI <- rbind(sim.results.TI, row)
  
  # Clear lb and ub variables
  rm('lb'); rm('ub')
}

sim.results.TI$t <- as.numeric(sim.results.TI$t)

# Select the link function
link <- "Cox_wb"

# Plot the results
sim.results.TI.sub <- sim.results.TI[sim.results.TI$link == link, ]
for (crit in sort(unique(sim.results.TI.sub$crit))) {
  
  # Select stratum criterion and subset data set.
  res.crit <- sim.results.TI.sub[which(sim.results.TI.sub$crit == crit), ]
  res.crit$t <- ifelse(res.crit$t == 2, 1, 5)
  
  # Select additional results for this stratum
  res.crit.base <- sim.results[sim.results$crit == crit & sim.results$link == link & sim.results$method %in% c("2", "4", "5", "7") & sim.results$var == "Age", c(1, 2, 3, 4, 5, 7, 8)]
  res.crit.base$t <- ifelse(res.crit.base$method == "2", 2, 
                     ifelse(res.crit.base$method == "4", 3,
                     ifelse(res.crit.base$method == "5", 6, 7)))
  colnames(res.crit.base) <- colnames(res.crit)
  
  res.crit <- rbind(res.crit, res.crit.base)
  
  # Initialize plot
  plt <- ggplot(mapping = aes(x = t), data = res.crit)
  
  # Line at y = 0 and x = 0
  plt <- plt + geom_segment(aes(x = -10, y = 0, xend = 40, yend = 0), colour = "black")
  plt <- plt + geom_segment(aes(x = 0, y = -10, xend = 0, yend = 10), colour = "black")
  
  # Plot identified intervals
  plt <- plt + geom_segment(aes(y = lb, xend = t, yend = ub))
  plt <- plt + coord_cartesian(xlim = c(-1, 10), ylim = c(-0.5, 1.5))
  
  print(plt)
}

