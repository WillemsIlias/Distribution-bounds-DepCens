
# Load dependencies
library(stringr)

#### Helper functions ####

#' @title Returns the last file/directory name of a given path.
#' 
#' @description This function returs the last file/directory name of a given
#' path. E.g. if one would supply the argument "dir1/dir2/filename", this
#' function returns "filename".
get.trailing.name <- function(path.vct) {
  
  # Initialize object that will store all of the trailing names
  trailing.names <- NULL
  
  # For each path in the given vector of paths...
  for (path in path.vct) {
    
    # Split the path into its components
    path.components <- strsplit(path, "/")[[1]]
    
    # Select and store the last component
    trailing.names <- c(trailing.names, path.components[length(path.components)])
  }
  
  # Return the results
  trailing.names
}

#' @title This function creates the path name based on the supplied arguments
#' 
make.path.name <- function(...) {
  paste(unlist(list(...)), collapse = "/")
}

#' @title This function parses the 'general information' (geninfo) part of a
#' name.
parse.geninfo <- function(geninfo) {
  
  # Get components of general info
  geninfo.components <- strsplit(strsplit(geninfo, "-")[[1]][2], "__")[[1]]
  geninfo.components <- as.numeric(geninfo.components)
  
  # Order of the components
  variables.ordered <- c("n", "n.cov", "DGP", "B", "t.eval", "seed", "sim.iter.nbr")
  
  # In an earlier version, 'DGP' was not recorded.
  if (length(geninfo.components) == (length(variables.ordered) - 1)) {
    variables.ordered <- variables.ordered[-3]
  }
  
  # Return the results
  geninfo.parsed <- as.list(geninfo.components)
  names(geninfo.parsed) <- variables.ordered
  geninfo.parsed
}

#' @title Search and - if found - parse the information belonging to a given
#' keyword
#'
#' @description This function searches for a given keyword in the vector of
#' components and parses the corresponding information.
search.and.parse <- function(components, keyword) {
  idx <- which(startsWith(components, keyword))
  if (length(idx) == 1) {
    return(strsplit(components[idx], "-")[[1]][2])
  } else {
    return(NULL)
  }
}

#' @title This function parses a directory/file name and returns the information
#' contained within it
#'
#' @description
#' ...
#' 
#' @param name String to parse
#' @param defaults default values for missing values. If \code{defaults = NULL},
#' missing values will not be returned.
#' 
parse.name <- function(name, defaults = NULL) {
  
  # Remove the device name (.csv) at the end
  name <- gsub(".csv", "", name)
  
  # Split the given name at the separators
  components <- strsplit(name, "__")[[1]]
  
  # Initialize a list that will store all the results
  info <- list()
  
  # Check if 'geninfo' is contained in the name. If so, parse it separately.
  if (any(startsWith(components, "geninfo"))) {
    
    # Get index of component that contains the name 'geninfo'
    geninfo.idx <- which(startsWith(components, "geninfo"))
    
    # Reconstruct general information from components
    geninfo <- paste(components[geninfo.idx:length(components)], collapse = "__")
    
    # Pass to the function that parses general information
    info <- c(info, parse.geninfo(geninfo))
  }
  
  # Search the string for all keywords
  keywords <- c("Search", "IF", "Kbar", "Alpha", "method", "nbpc", "nifpc",
                "degree", "Gc", "LF", "Param")
  keywords.handled <- NULL
  for (keyword in keywords) {
    keyword.info <- search.and.parse(components, keyword)
    if (!is.null(keyword.info)) {
      keywords.handled <- c(keywords.handled, keyword)
      names.info <- names(info)
      info <- c(info, keyword.info)
      names(info) <- c(names.info, keyword)
    }
  }
  
  # Look for defaults values of keywords that weren't present in the string
  for (keyword in setdiff(keywords, keywords.handled)) {
    if (keyword %in% names(defaults)) {
      info[keyword] <- defaults[keyword]
    }
  }
  
  # Return the results
  return(info)
}

#' @title Appends a simulation setting to the matrix that collects all
#' simulation settings.
#' 
#' @description This function appends the given simulation setting to the matrix
#' of simulation settings. This function will add extra columns to the matrix if
#' necessary.
append.sim.setting <- function(sim.settings, sim.info) {
  
  # Row index at which to append the information, if necessary
  row.idx <- nrow(sim.settings) + 1
  
  # Flag whether sim.info should be appended
  append <- FALSE
  
  # Remove seed info from sim.info
  sim.info[["seed"]] <- NULL
  
  # Check if all elements of sim.info are columns in sim.settings. If not, the
  # missing columns should be created and sim.info appended
  if (!all(names(sim.info) %in% colnames(sim.settings))) {
    append <- TRUE
  }
  
  # If the previous if-clause did not run...
  if (!append) {
    
    # Check whether sim.info is contained in sim.settings
    if (!(any(apply(X = matrix(sim.info %in% sim.settings, nrow = nrow(sim.settings)), 1, all)))) {
      append <- TRUE
    }
  }
  
  # If sim.info should be appended, do so.
  if (append) {
    for (elem in names(sim.info)) {
      sim.settings[row.idx, elem] <- sim.info[[elem]]
    }
  }
  
  sim.settings
}

#' @title This function reconstructs the directory name from a given simulation
#' setting.
reconstruct.dirname <- function(setting, shortened, defaults) {
  
  # Reconstruct full dir name
  dir.name <- paste(c(paste0("LF-", setting[, "LF"]),
                      paste0("Search-", setting[, "Search"]),
                      paste0("IF-", setting[, "IF"]),
                      paste0("Kbar-", setting[, "Kbar"])),
                    collapse = "__")
  if (!is.na(setting[, "Alpha"])) {
    dir.name <- paste(c(dir.name, paste0("Alpha-", setting[, "Alpha"])),
                      collapse = "__")
  }
  if ("Param" %in% colnames(setting)) {
    dir.name <- sprintf("%s__Param-%s", dir.name, setting[, "Param"])
  }
  
  # If necessary, remove defaults
  if (shortened) {
    components <- str_split(dir.name, "__")[[1]]
    components.shortened <- NULL
    for (component in components) {
      if (!(str_split(component, "-")[[1]][1] %in% names(defaults))) {
        components.shortened <- c(components.shortened, component)
      }
    }
    
    dir.name <- paste(components.shortened, collapse = "__")
  }
  
  # Return dir name
  dir.name
}

#' @title This function returns a regular expression for the names of the
#' simulation files corresponding to the provided setting.
get.filename.regex <- function(setting, shortened, defaults) {
  
  # Initialize variables
  seed.regex <- "[[:digit:]]+"
  file.regex <- NULL
  
  # If the simulation was run using grid search, add the gridsearch method to
  # file name
  if (!is.na(setting[, "method"])) {
    file.regex <- paste(c(file.regex,
                          paste0("method-", setting[, "method"])),
                        collapse = "__")
  }
  
  # Depending on which type of instrumental function family was used, add the
  # necessary information to file name
  if ("nifpc" %in% colnames(setting)) {
    if (!is.na(setting[, "nifpc"])) {
      file.regex <- paste(c(file.regex,
                            paste0("nifpc-", setting["nifpc"])),
                          collapse = "__")
    }
    if (!is.na(setting[, "Gc"])) {
      file.regex <- paste(c(file.regex,
                            paste0("Gc-", setting["Gc"])),
                          collapse = "__")
    }
    if (!is.na(setting[, "degree"])) {
      file.regex <- paste(c(file.regex,
                            paste0("degree-", setting["degree"])),
                          collapse = "__")
    }
  }
  if ("nbpc" %in% colnames(setting)) {
    if (!is.na(setting[, "nbpc"])) {
      file.regex <- paste(c(file.regex,
                            paste0("nbpc-", setting["nbpc"])),
                          collapse = "__")
    }
  }
  
  # Add general info
  file.regex <- 
    paste(c(file.regex,
            paste0("geninfo-", paste(c(setting[, "n"],
                                        setting[, "n.cov"],
                                        if ("DGP" %in% colnames(setting)) setting[, "DGP"] else NULL,
                                        setting[, "B"],
                                        setting[, "t.eval"],
                                        seed.regex,
                                        setting[, "sim.iter.nbr"]),
                                      collapse = "__"))),
            collapse = "__")
  
  # If necessary, remove defaults
  if (shortened) {
    components <- str_split(file.regex, "__")[[1]]
    components.shortened <- NULL
    for (component in components) {
      if (!(str_split(component, "-")[[1]][1] %in% names(defaults))) {
        components.shortened <- c(components.shortened, component)
      }
    }
    
    file.regex <- paste(components.shortened, collapse = "__")
  }
  
  # Add file suffix (comma separated file)
  file.regex <- paste(c(file.regex, "csv"), collapse = ".")
  
  # Return file regex
  file.regex
}

#### LaTeX formatting functions ####

#' @title Create LaTeX environment
#' 
#' @description
#' ...
#' 
#' @param env.name Name of LaTeX environment
#' @param arguments Arguments to the environment, put between curly brackets.
#' @param options Options to the environment, put between square brackets.
#' 
make.env <- function(env.name, arguments, options, ...) {
  
  # Create variable that will store the LaTeX code to generate the environment
  env <- NULL
  
  # Environment header
  header <- sprintf("\\begin{%s}", env.name)
  if (!is.null(arguments)) {
    header <- sprintf("%s{%s}", header, arguments)
  }
  if (!is.null(options)) {
    header <- sprintf("%s[%s]", header, options)
  }
  
  # Environment body
  body <- paste("   ", unlist(list(...)))
  
  # Environment footer
  footer <- sprintf("\\end{%s}", env.name)
  
  # return everything
  c(header, body, footer)
}

latex.caption <- function(caption) {
  sprintf("\\caption{%s}", caption)
}

latex.label <- function(caption) {
  sprintf("\\label{%s}", caption)
}

get.results.matrix <- function(sim.results.subset, nifpc1 = NULL, nifpc2 = NULL) {
  
  # Obtain the 2 unique values of nifpc
  if (is.null(nifpc1)) {
    nifpc1 <- sort(unique(sim.results.subset$nifpc))[1]
  }
  if (is.null(nifpc2)) {
    nifpc2 <- sort(unique(sim.results.subset$nifpc))[2]
  }
  nifpc.vec <- c(nifpc1, nifpc2)
  
  results.matrix <- matrix(NA, nrow = 3, ncol = 10)
  for (nifpc.idx in 1:2) {
    
    # Subset to the results pertaining to this iteration
    sim.results.nifpci <- sim.results.subset[sim.results.subset$nifpc == nifpc.vec[nifpc.idx], ]
    
    # Order based in increasing sample size
    sim.results.nifpci <- sim.results.nifpci[order(sim.results.nifpci$n), ]
    
    # Extract necessary variables
    if (nifpc.idx == 1) {
      results.matrix[, 1:5] <- as.matrix(sim.results.nifpci[, c("avg.l", "avg.u", "var.width", "sig", "coverage")])
    } else {
      results.matrix[, 6:10] <- as.matrix(sim.results.nifpci[, c("avg.l", "avg.u", "var.width", "sig", "coverage")])
    }
  }
  
  # Round everything two 2 decimal places
  results.matrix <- format(round(results.matrix, digits = 2), nsmall = 2) 
  
  # Return the matrix
  results.matrix
}

get.results.matrix.MT <- function(sim.results.subset) {
  
  results.matrix <- sim.results.subset[1, c("avg.l", "avg.u", "sig", "coverage",
                                            "avg.l.int", "avg.u.int", "sig.int", "coverage.int",
                                            "avg.l.MV", "avg.u.MV", "sig.MV", "coverage.MV")]
  
  # Round everything two 2 decimal places
  results.matrix <- format(round(results.matrix, digits = 2), nsmall = 2)
  
  # Return the matrix
  results.matrix
}

get.sim.setting.latex <- function(DGP) {
  
  if (DGP == 21) {
    return("Indep.\\\\ $\\sim 30\\%$ cens.")
  } else if (DGP == 22) {
    return("Pos. dep.\\\\ $\\sim 30\\%$ cens.")
  } else if (DGP == 23) {
    return("Neg. dep.\\\\ $\\sim 30\\%$ cens.")
  } else if (DGP == 24) {
    return("Indep.\\\\ $\\sim 65\\%$ cens.")
  } else if (DGP == 25) {
    return("Pos. dep.\\\\ $\\sim 65\\%$ cens.")
  } else if (DGP == 26) {
    return("Neg. dep.\\\\ $\\sim 65\\%$ cens.")
  } else if (DGP == 27) {
    return("Indep.\\\\ $\\sim 30\\%$ cens.")
  } else if (DGP == 28) {
    return("Pos. dep.\\\\ $\\sim 30\\%$ cens.")
  } else if (DGP == 29) {
    return("Neg. dep.\\\\ $\\sim 30\\%$ cens.")
  } else if (DGP == 30) {
    return("Indep.\\\\ $\\sim 65\\%$ cens.")
  } else if (DGP == 31) {
    return("Pos. dep.\\\\ $\\sim 65\\%$ cens.")
  } else if (DGP == 32) {
    return("Neg. dep.\\\\ $\\sim 65\\%$ cens.")
  } else if (DGP == 33) {
    return("Indep.\\\\ $\\sim 2\\%$ cens.")
  } else if (DGP == 34) {
    return("Pos. dep.\\\\ $\\sim 2\\%$ cens.")
  } else if (DGP == 35) {
    return("Neg. dep.\\\\ $\\sim 2\\%$ cens.")
  } else if (DGP == 36) {
    return("Indep.\\\\ $\\sim 2\\%$ cens.")
  } else if (DGP == 37) {
    return("Pos. dep.\\\\ $\\sim 2\\%$ cens.")
  } else if (DGP == 38) {
    return("Neg. dep.\\\\ $\\sim 2\\%$ cens.")
  }
}

matrix2latex <- function(results.matrix.list) {
  
  # Set useful variable
  n.vec <- c(500, 1000, 2000)
  latex.matrix.full <- NULL
  
  for (rml.name in names(results.matrix.list)) {
    
    # Get results matrix of this iteration
    results.matrix <- results.matrix.list[[rml.name]]
    
    # Create the latex code.
    latex.matrix.part <- NULL
    for (i in 1:nrow(results.matrix)) {
      part1.bounds <- sprintf("[%s, %s]", results.matrix[i, 1], results.matrix[i, 2])
      part1 <- paste(c(part1.bounds, results.matrix[i, 3:5]), collapse = " & ")
      
      part2.bounds <- sprintf("[%s, %s]", results.matrix[i, 6], results.matrix[i, 7])
      part2 <- paste(c(part2.bounds, results.matrix[i, 8:10]), collapse = " & ")
      row <- sprintf("& %s & %s && %s\\\\", n.vec[i], part1, part2)
      if (i == 1) {
        setting.info.latex <- sprintf("\\multirow{3}{*}{\\shortstack[c]{%s}}",
                                      get.sim.setting.latex(as.numeric(rml.name)))
        row <- sprintf("%s %s", setting.info.latex, row)
      }
      latex.matrix.part <- rbind(latex.matrix.part, row)
    }
    
    latex.matrix.full <- rbind(latex.matrix.full, latex.matrix.part, "\\midrule")
  }
  
  # Remove last redundant midrule
  latex.matrix.full <- latex.matrix.full[-length(latex.matrix.full)]
}

matrix2latex.MT <- function(results.matrix.list) {
  
  # Set useful variable
  latex.matrix.full <- NULL
  
  for (rml.name in names(results.matrix.list)) {
    
    # Get results matrix of this iteration
    results.matrix <- results.matrix.list[[rml.name]]
    
    # Parse rml.name
    DGP <- as.numeric(str_split(rml.name, "_")[[1]][1])
    t <- as.numeric(str_split(rml.name, "_")[[1]][2])
    
    # Create the latex code.
    latex.matrix.part <- NULL
    for (i in 1:nrow(results.matrix)) {
      part1.bounds <- sprintf("[%s, %s]", results.matrix[i, 1], results.matrix[i, 2])
      part1 <- paste(c(part1.bounds, results.matrix[i, 3:4]), collapse = " & ")
      
      part2.bounds <- sprintf("[%s, %s]", results.matrix[i, 5], results.matrix[i, 6])
      part2 <- paste(c(part2.bounds, results.matrix[i, 7:8]), collapse = " & ")
      
      part3.bounds <- sprintf("[%s, %s]", results.matrix[i, 9], results.matrix[i, 10])
      part3 <- paste(c(part3.bounds, results.matrix[i, 11:12]), collapse = " & ")
      
      row <- sprintf("& %s & %s && %s && %s\\\\", DGP, part1, part2, part3)
      
      if (DGP == 21 | DGP == 27) {
        setting.info.latex <- sprintf("\\multirow{6}{*}{\\shortstack[c]{%s}}", t)
        row <- sprintf("%s %s", setting.info.latex, row)
      }
      latex.matrix.part <- rbind(latex.matrix.part, row)
      latex.matrix.full <- rbind(latex.matrix.full, latex.matrix.part)
    }
    
    # If at the end of section with t = 1, put midrule
    if (DGP %in% c(26, 32) & t == 1) {
      latex.matrix.full <- rbind(latex.matrix.full, "\\midrule")
    }
  }
  
  # Return the results
  latex.matrix.full
}

print.latex <- function(latex.code) {
  message(paste(latex.code, collapse = "\n"))
}


#### Main functions ####

# Obtain all simulation settings (based on file names)
get.sim.settings <- function(master.dir, defaults) {
  
  # Try to recover simulation settings from design matrix
  removed._test <- FALSE
  if (grepl("_test", master.dir)) {
    master.dir <- gsub("_test", "", master.dir)
    removed._test <- TRUE
  }
  if (master.dir == "Simulations_Main") {
    design.matrix.dir <- "simMainDesigns"
    design.matrix.subdir <- "simMainDesigns"
  } else if (master.dir == "Simulations_MainMiss") {
    design.matrix.dir <- "simMainMissDesigns"
    design.matrix.subdir <- "simMainMissDesigns"
  } else if (master.dir == "Simulations_Main_almostNoCens") {
    design.matrix.dir <- "simMainDesigns_almostNoCens"
    design.matrix.subdir <- "simMainDesigns"
  } else if (master.dir == "Simulations_Main_depCov") {
    design.matrix.dir <- "simMainDesigns_dependentCovariates"
    design.matrix.subdir <- "simMainDesigns"
  } else if (master.dir == "Simulations_Main_moreIF") {
    design.matrix.dir <- "simMainDesigns_moreIF"
    design.matrix.subdir <- "simMainDesigns"
  } else if (master.dir == "Simulations_Main_MT") {
    design.matrix.dir <- "simMainDesigns_MT"
    design.matrix.subdir <- "simMainDesigns"
  } else if (master.dir == "Simulations_Main_manyCov") {
    design.matrix.dir <- "simMainDesigns_manyCov"
    design.matrix.subdir <- "simMainDesigns"
  }
  if (removed._test) {
    master.dir <- paste0(master.dir, "_test")
  }
  
  # If the design matrix can be found...
  if (design.matrix.dir %in% list.dirs(getwd(), recursive = FALSE, full.names = FALSE)) {
    
    # Read design matrix
    sim.design.matrix <- read.csv(sprintf("%s/%s1.csv", design.matrix.dir,
                                          design.matrix.subdir))
    
    # Remove seed information
    sim.design.matrix <- sim.design.matrix[, -which(colnames(sim.design.matrix) == "seed")]
    
    # Initialize data frame to store the results
    sim.settings <- data.frame()
    
    # Create data frame of simulation designs
    for (row.idx in 1:nrow(sim.design.matrix)) {
      
      # Extract all information
      if ("parametric" %in% colnames(sim.design.matrix)) {
        LF <- ifelse(sim.design.matrix[row.idx, "link_function"] == 1, "Cox_wb", "AFT_ll")
      } else {
        LF <- ifelse(sim.design.matrix[row.idx, "link_function"] == 1, "AFT_ll", "Cox_wb")
      }
      Search <- defaults$Search
      IF <- defaults$IF
      Kbar <- defaults$Kbar
      Alpha <- defaults$Alpha
      method <- defaults$method
      degree <- defaults$degree
      sim.iter.nbr <- 1
      Gc <- defaults$Gc
      n <- sim.design.matrix[row.idx, "n"]
      n.cov <- ifelse("n_cov" %in% colnames(sim.design.matrix), sim.design.matrix[row.idx, "n_cov"], 2)
      nifpc <- sim.design.matrix[row.idx, "n_if_per_cov"]
      DGP <- sim.design.matrix[row.idx, "DGP"]
      B <- 600
      t.eval <- sim.design.matrix[row.idx, "t_eval"]
      Param <- NULL
      if ("parametric" %in% colnames(sim.design.matrix)) {
        Param <- sim.design.matrix[row.idx, "parametric"]
      }
      if ("inst_func_family_idx" %in% colnames(sim.design.matrix)) {
        if (sim.design.matrix[row.idx, "inst_func_family_idx"] == 1) {
          IF <- "box"
        } else if (sim.design.matrix[row.idx, "inst_func_family_idx"] == 2) {
          IF <- "spline"
        } else if (sim.design.matrix[row.idx, "inst_func_family_idx"] == 3) {
          IF <- "cd"
        } else if (sim.design.matrix[row.idx, "inst_func_family_idx"] == 4) {
          IF <- "cdmc"
        }
      }
      
      sim.info <- list(
        LF = LF,
        Param = Param,
        n = n,
        DGP = DGP,
        nifpc = nifpc,
        Search = Search,
        IF = IF,
        Kbar = Kbar,
        Alpha = Alpha,
        method = method,
        degree = degree,
        Gc = Gc,
        n.cov = n.cov,
        B = B,
        sim.iter.nbr = sim.iter.nbr,
        t.eval = t.eval
      )
      sim.info <- sim.info[!sapply(sim.info, is.null)]
      
      # Write results to data frame
      for (elem in names(sim.info)) {
        sim.settings[row.idx, elem] <- sim.info[[elem]]
      }
    }
    
    # Return the results
    return(unique(sim.settings))
  }
  
  # Initialize matrix that will store all simulation settings.
  sim.settings <- data.frame()
  
  # Get all sub-directories of the master directory
  sub.dirs <- setdiff(get.trailing.name(list.dirs(master.dir)), master.dir)
  
  # For each file in each subdirectory, obtain the simulation setting and append
  # it to the matrix that collects all simulation settings.
  for (sub.dir in sub.dirs) {
    
    # Parse directory name
    dir.info <- parse.name(sub.dir, defaults)
    files <- list.files(make.path.name(master.dir, sub.dir))
    
    # For each file in the directory...
    for (i in 1:length(files)) {
      
      # Update user
      if (i %% 100 == 0) {
        message(paste0(round(100 * i/length(files), 2), "% completion"))
      }
      
      # Parse file name
      file <- files[i]
      sim.info <- c(dir.info, parse.name(file, defaults))
      
      # If necessary, append setting to matrix
      sim.settings <- append.sim.setting(sim.settings, sim.info)
    }
  }
  
  # Return all unique simulation settings
  unique(subset(sim.settings, select = -seed))
}

#' @title  Obtain all file names for a specific simulation setting
#' 
#' @description
#' Reconstructs the file names based on the simulation design
#' 
#' @param master.dir
#' @param setting
#' @param shortened Boolean variable that indicates whether a shortened vversion
#' of the path name should be provided (omitting all defaults).
#' @param defaults
#' 
get.sim.files <- function(master.dir, setting, shortened, defaults) {
  
  # Reconstruct the directory and file name from the given settings
  dir.name <- reconstruct.dirname(setting, shortened, defaults)
  file.regex <- get.filename.regex(setting, shortened, defaults)
  
  # In the appropriate directory, find all file names that match the regex
  files.in.dir <- list.files(make.path.name(master.dir, dir.name))
  sim.files <- files.in.dir[grepl(file.regex, files.in.dir)]
  
  # Return the results
  paste(master.dir, dir.name, sim.files, sep = "/")
}

# Summarize the files supplied as argument to this function
summarize.files <- function(sim.files, beta.true, idx.param.of.interest,
                            setting.idx = "", files.df = NULL) {
  
  #### Obtain data frame of all simulation results ####
  
  if (is.null(files.df)) {
    files.df <- read.csv(sim.files[1])
    for (file.name in sim.files[-1]) {
      files.df <- rbind(files.df, read.csv(file.name))
    }
  }
  
  #### Perform some checks ####
  
  # There should be no duplicate seeds
  if (length(unique(files.df$seed)) != nrow(files.df)) {
    stop("Duplicate random seed(s) detected!")
  }
  
  # Due to a slight bug in a previous version of the code, it was possible for
  # the identified intervals to exceed the specified bounds. If this happened,
  # clip the bounds to their specified values.
  files.df$ident.set.l <- pmax(files.df$ident.set.l, -10)
  files.df$ident.set.u <- pmin(files.df$ident.set.u, 10)
  
  #### Compute summary statistics ####
  
  # Record the number of times no initial value was found (in which case one
  # would conclude that the model is misspecified).
  msp <- table(as.numeric((files.df$conv.l == 2) & (files.df$conv.u == 2)))["1"]
  if (is.na(msp)) {
    msp <- 0
  }
  
  # Record the number of times an empty instrumental function occured (in which
  # case the method could not be executed).
  eir <- table(as.numeric((files.df$conv.l == 3) & (files.df$conv.u == 3)))["1"]
  if (is.na(eir)) {
    eir <- 0
  }
  
  # Subset the data to only converged iterations. Store information about the
  # iterations that did not converge
  conv <- table(as.numeric((files.df$conv.l == 1) & (files.df$conv.u == 1)))["1"]
  files.df <- files.df[which((files.df$conv.l == 1) & (files.df$conv.u == 1)), ]
  
  # For all simulations that either produced an interval or concluded model
  # misspecification, compute the proportion of model misspecifications
  msp.prop <- msp / (msp + conv)
  
  # Get average lower and upper bounds
  avg.l <- mean(files.df$ident.set.l)
  avg.u <- mean(files.df$ident.set.u)
  
  # Get average and variance of width of the interval
  widths <- files.df$ident.set.u - files.df$ident.set.l
  avg.width <- mean(widths)
  var.width <- var(widths)
  
  # Compute the percentage of times the bounds do not contain zero
  sig <- sum(!((files.df$ident.set.l < 0) & (0 < files.df$ident.set.u)))/(conv)
  
  # Get coverage percentage
  coverage <- sum((files.df$ident.set.l <= beta.true[idx.param.of.interest]) &
                  (files.df$ident.set.u >= beta.true[idx.param.of.interest]))/(conv)
  
  # Get percentage of times the independence model parameter is contained in the
  # identified interval
  indep.coverage <- sum((files.df$ident.set.l <= files.df$coef.ident.model) &
                        (files.df$ident.set.u >= files.df$coef.ident.model))/(conv)
  
  # Get average total run time
  avg.runtime <- mean(files.df$total.run.time)
  
  # If applicable, get time spent on each step of EAM algorithm.
  avg.Estep <- NULL
  avg.Astep <- NULL
  avg.Mstep <- NULL
  if (all(c("E-step", "A-step", "M-step") %in% colnames(files.df))) {
    avg.Estep <- mean(files.df$`E-step`)
    avg.Astep <- mean(files.df$`A-step`)
    avg.Mstep <- mean(files.df$`M-step`)
  }
  
  # Get average censoring percentage
  avg.per.cens <- mean(files.df$per.cens)
  
  # Correlation of interval width on censoring percentage
  cor.cens.width <- cor(files.df$per.cens, widths)
  
  #### Return the results ####
  
  summary.vals <- c(avg.l, avg.u, avg.width, var.width, sig, coverage, avg.per.cens,
                    cor.cens.width,  avg.runtime, conv, msp, eir, indep.coverage, 
                    avg.Estep, avg.Astep, avg.Mstep)
  summary.names <- c("avg.l", "avg.u", "avg.width", "var.width", "sig", "coverage",
                     "avg.per.cens", "cor.cens.width",  "avg.runtime", "conv",
                     "msp", "eir", "indep.cov", "avg.Estep", "avg.Astep",
                     "avg.Mstep")[1:length(summary.vals)]
  summary <- matrix(nrow = 1, ncol = length(summary.vals))
  summary[1, ] <- summary.vals
  colnames(summary) <- summary.names
  
  summary
}

# Function to construct combined bounds based on majority vote
combine.bounds.MV <- function(lbs, ubs) {
  
  # Get bounds of each part
  ths <- sort(unique(c(lbs, ubs)))
  parts <- matrix(rep(ths, c(1, rep(2, length(ths) - 2), 1)), ncol = 2, byrow = TRUE)
  
  # For each part, obtain the number of votes
  parts <- cbind(parts, rep(0, nrow(parts)))
  for (part.idx in 1:nrow(parts)) {
    parts[part.idx, 3] <- sum((lbs < mean(parts[part.idx, 1:2])) & (ubs > mean(parts[part.idx, 1:2])))
  }
  
  # Subset to parts getting majority (> 0.5) vote
  bounds.MV <- parts[parts[, 3] > floor(length(lbs)/2), 1:2, drop = FALSE]
  
  # Recombine neighbouring bounds
  row.idx <- 1
  while (row.idx < nrow(bounds.MV)) {
    if (bounds.MV[row.idx, 2] == bounds.MV[row.idx + 1, 1]) {
      bounds.MV[row.idx, 2] <- bounds.MV[row.idx + 1, 2]
      bounds.MV <- bounds.MV[setdiff(1:nrow(bounds.MV), row.idx + 1), , drop = FALSE]
    } else {
      row.idx <- row.idx + 1
    }
  }
  
  # Give proper column names
  colnames(bounds.MV) <- c("lb", "ub")
  
  # Return the result
  bounds.MV
}

# Summarize the files coming from the simulations regarding multiple testing.
summarize.files.MT <- function(sim.files, beta.true, idx.param.of.interest,
                               setting.idx = "") {
  
  #### Obtain data frame of all simulation results ####
  
  files.df <- read.csv(sim.files[1])
  for (file.name in sim.files[-1]) {
    files.df <- rbind(files.df, read.csv(file.name))
  }
  
  #### Results without multiple testing ####
  
  summary.noMT <- summarize.files(sim.files, beta.true, idx.param.of.interest,
                                  setting.idx, files.df)
  
  #### Results with multiple testing ####
  
  # Construct combined bounds: intersection.
  lb.intersect <- pmax(files.df$ident.set.l_int.1, files.df$ident.set.l_int.2, files.df$ident.set.l_int.3)
  ub.intersect <- pmin(files.df$ident.set.u_int.1, files.df$ident.set.u_int.2, files.df$ident.set.u_int.3)
  
  # Construct combined convergence information: intersection.
  conv.l.intersect <- ifelse((files.df$conv.l_int.1 == 1 | files.df$conv.l_int.2 == 1 | files.df$conv.l_int.3 == 1), 1, 0)
  conv.u.intersect <- ifelse((files.df$conv.u_int.1 == 1 | files.df$conv.u_int.2 == 1 | files.df$conv.u_int.3 == 1), 1, 0)
  conv.intersect <- sum(conv.l.intersect & conv.u.intersect)
  msp.intersect <- sum((pmin(files.df$conv.l_int.1, files.df$conv.l_int.2, files.df$conv.l_int.3) == 2) &
                         (pmin(files.df$conv.u_int.1, files.df$conv.u_int.2, files.df$conv.u_int.3) == 2))
  
  # Subset results to those that converged
  lb.intersect <- lb.intersect[conv.l.intersect & conv.u.intersect]
  ub.intersect <- ub.intersect[conv.l.intersect & conv.u.intersect]
  
  # Construct results for intersected bounds
  avg.l.intersect <- mean(lb.intersect)
  avg.u.intersect <- mean(ub.intersect)
  widths.intersect <- ub.intersect - lb.intersect
  var.width.intersect <- var(widths.intersect)
  
  sig.intersect <- sum(!((lb.intersect < 0) & (0 < ub.intersect)))/(conv.intersect)
  coverage.intersect <- sum((lb.intersect <= beta.true[idx.param.of.interest]) &
                              (ub.intersect >= beta.true[idx.param.of.interest]))/(conv.intersect)
  
  #### Results with majority voting ####
  
  # Construct combined bounds: majority vote
  bounds.MV <- matrix(nrow = 0, ncol = 2)
  for (i in 1:nrow(files.df)) {
    lbs <- c(files.df[i, "ident.set.l_mv.1"], files.df[i, "ident.set.l_mv.2"], files.df[i, "ident.set.l_mv.3"])
    ubs <- c(files.df[i, "ident.set.u_mv.1"], files.df[i, "ident.set.u_mv.2"], files.df[i, "ident.set.u_mv.3"])
    bounds.MV.i <- combine.bounds.MV(lbs, ubs)
    if (nrow(bounds.MV.i) == 1) {
      bounds.MV <- rbind(bounds.MV, bounds.MV.i)
    } else {
      stop("Combined set via majority vote is not an interval!")
    }
  }
  lb.MV <- bounds.MV[, 1]
  ub.MV <- bounds.MV[, 2]
  
  # Construct combined convergence information: intersection.
  conv.l.MV <- ifelse((files.df$conv.l_mv.1 == 1 | files.df$conv.l_mv.2 == 1 | files.df$conv.l_mv.3 == 1), 1, 0)
  conv.u.MV <- ifelse((files.df$conv.u_mv.1 == 1 | files.df$conv.u_mv.2 == 1 | files.df$conv.u_mv.3 == 1), 1, 0)
  conv.MV <- sum(conv.l.MV & conv.u.MV)
  msp.MV <- sum((pmin(files.df$conv.l_mv.1, files.df$conv.l_mv.2, files.df$conv.l_mv.3) == 2) &
                         (pmin(files.df$conv.u_mv.1, files.df$conv.u_mv.2, files.df$conv.u_mv.3) == 2))
  
  # Subset results to those that converged
  lb.MV <- lb.MV[conv.l.MV & conv.u.MV]
  ub.MV <- ub.MV[conv.l.MV & conv.u.MV]
  
  # Construct results for intersected bounds
  avg.l.MV <- mean(lb.MV)
  avg.u.MV <- mean(ub.MV)
  widths.MV <- ub.MV - lb.MV
  var.width.MV <- var(widths.MV)
  
  sig.MV <- sum(!((lb.MV < 0) & (0 < ub.MV)))/(conv.MV)
  coverage.MV <- sum((lb.MV <= beta.true[idx.param.of.interest]) &
                              (ub.MV >= beta.true[idx.param.of.interest]))/(conv.MV)

  #### Return the results ####
  
  summary.vals <- c(summary.noMT[, "avg.l"], summary.noMT[, "avg.u"],
                    summary.noMT[, "var.width"], summary.noMT[, "sig"],
                    summary.noMT[, "coverage"], summary.noMT[, "conv"],
                    summary.noMT[, "msp"], summary.noMT[, "eir"],
                    avg.l.intersect, avg.u.intersect, var.width.intersect,
                    sig.intersect, coverage.intersect, conv.intersect, msp.intersect,
                    avg.l.MV, avg.u.MV, var.width.MV,
                    sig.MV, coverage.MV, conv.MV, msp.MV,
                    summary.noMT[, "avg.per.cens"])
  summary.names <- c("avg.l", "avg.u", "var.width", "sig", "coverage", "conv", "msp", "eir",
                     "avg.l.int", "avg.u.int", "var.width.int", "sig.int", "coverage.int", "conv.int", "msp.int",
                     "avg.l.MV", "avg.u.MV", "var.width.MV", "sig.MV", "coverage.MV", "conv.MV", "msp.MV",
                     "avg.per.cens")
  summary <- matrix(nrow = 1, ncol = length(summary.vals))
  summary[1, ] <- summary.vals
  colnames(summary) <- summary.names
  
  summary
}

#' @title Produce LaTeX table of results
#' 
#' @description
#' This function produces LaTeX code that generates a table of the provided
#' simulation results.
#' 
#' @param sim.results.subset Simulation results
#' @param nifpc.options Numbers (2) of instrumental functions for which to
#' include the results in the output matrix. Default is
#' \code{nifpc.options = NULL}, in which case it is assumed that only 2 values
#' for it are possible (obtained from \code{sim.results.subset}).
#' 
results2latex <- function(sim.results.subset, caption, label,
                          nifpc.options = NULL) {
  
  #### Extract results ####
  
  idxs.numeric.variables <- which(!(colnames(sim.results.subset) %in% c("Gc", "IF")))
  sim.results.subset[, idxs.numeric.variables] <- apply(sim.results.subset[,idxs.numeric.variables], 2, as.numeric)
  
  if (is.null(nifpc.options)) {
    if (length(unique(sim.results.subset$nifpc)) == 2) {
      nifpc1 <- sort(unique(sim.results.subset$nifpc))[1]
      nifpc2 <- sort(unique(sim.results.subset$nifpc))[2]
    } else {
      stop("When number of unique nifpc is not equal to 2, nifpc.options should be specified.")
    }
  } else {
    nifpc1 <- nifpc.options[1]
    nifpc2 <- nifpc.options[2]
  }
  
  unique.DGPs <- unique(sim.results.subset$DGP)
  results.matrix.list <- list()
  for (DGP in sort(unique.DGPs)) {
    results.matrix.list[[as.character(DGP)]] <- get.results.matrix(sim.results.subset[sim.results.subset$DGP == DGP, ],
                                                                   nifpc1, nifpc2)
  }
  
  #### Make table ####
  
  tabular.arguments <- paste(rep("c", 11), collapse = "")
  toprule <- "\\toprule"
  title1 <- sprintf("& & \\multicolumn{4}{c}{nifpc = %s} & \\phantom{c} & \\multicolumn{4}{c}{nifpc = %s}\\\\",
                    nifpc1, nifpc2)
  title.sep <- "\\cmidrule{3-6} \\cmidrule{8-11}"
  title2 <- paste0("Setting & n ", paste(rep("& Bounds & Var & Sig & Cov", 2), collapse = " & "), "\\\\")
  tabular.body <- matrix2latex(results.matrix.list)
  bottomrule <- "\\bottomrule"
  midrule <- "\\midrule"
  
  latex.code <- make.env("table", NULL, "h",
                         "\\centering",
                         make.env("tabular", tabular.arguments, NULL,
                                  toprule,
                                  title1,
                                  title.sep,
                                  title2,
                                  midrule,
                                  tabular.body,
                                  bottomrule),
                         latex.caption(caption),
                         latex.label(label)
  )
  
  #### Print table to console ####
  
  print.latex(latex.code)
}

#' @title Produce LaTeX table of results for multiple testing.
#' 
#' @description
#' This function produces LaTeX code that generates a table of the provided
#' simulation results (specific function for multiple testing simulations).
#' 
#' @param sim.results.subset Simulation results
#' @param nifpc.options Numbers (2) of instrumental functions for which to
#' include the results in the output matrix. Default is
#' \code{nifpc.options = NULL}, in which case it is assumed that only 2 values
#' for it are possible (obtained from \code{sim.results.subset}).
#' 
results2latex.MT <- function(sim.results.subset, caption, label) {
  
  #### Extract results ####
  
  idxs.numeric.variables <- which(!(colnames(sim.results.subset) %in% c("Gc", "IF")))
  sim.results.subset[, idxs.numeric.variables] <- apply(sim.results.subset[,idxs.numeric.variables], 2, as.numeric)
  
  unique.DGPs <- unique(sim.results.subset$DGP)
  results.matrix.list <- list()
  for (t in sort(unique(sim.results.subset$t.eval))) {
    for (DGP in sort(unique.DGPs)) {
      results.matrix.list[[paste0(as.character(DGP), "_", as.character(t))]] <-
        get.results.matrix.MT(sim.results.subset[sim.results.subset$DGP == DGP &
                                                 sim.results.subset$t.eval == t, ])
    }
  }
  
  #### Make table ####
  
  tabular.arguments <- paste(rep("c", 13), collapse = "")
  toprule <- "\\toprule"
  title1 <- sprintf("& & \\multicolumn{3}{c}{%s} & \\phantom{c} & \\multicolumn{3}{c}{%s} & \\phantom{c} & \\multicolumn{3}{c}{%s}\\\\",
                    "Single point", "Intersection", "Majority vote")
  title.sep <- "\\cmidrule{3-5} \\cmidrule{7-9} \\cmidrule{11-13}"
  title2 <- paste0("t & DGP ", paste(rep("& Bounds & Sig & Cov", 3), collapse = " & "), "\\\\")
  tabular.body <- matrix2latex.MT(results.matrix.list)
  bottomrule <- "\\bottomrule"
  midrule <- "\\midrule"
  
  latex.code <- make.env("table", NULL, "h",
                         "\\centering",
                         make.env("tabular", tabular.arguments, NULL,
                                  toprule,
                                  title1,
                                  title.sep,
                                  title2,
                                  midrule,
                                  tabular.body,
                                  bottomrule),
                         latex.caption(caption),
                         latex.label(label)
  )
  
  #### Print table to console ####
  
  print.latex(latex.code)
}










