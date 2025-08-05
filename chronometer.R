
# Load dependencies
library(lubridate)
library(R6)

#' @title Chronometer object
#' 
#' @description
#' R6 object that mimics a chronometer. It can be started, paused, record legs
#' and stopped.
#'
#' @method show: displays information stored in chronometer
#' @method reset: resets the chronometer
#' @method start: Starts the chronometer and records the starting time
#' @method stop: Stops the chronometer. If legs were being recorded, also
#' computes the last leg time.
#' @method record.leg: Records a leg time
#' 
#' @import lubridate
#' @import R6
Chronometer <- R6Class("Chronometer",
                       
  public = list(
    
    # Print values of the chronometer
    show = function() {
      print(paste0("Start: ", private$start.time.abs))
      if (length(private$leg.times.rel) != 0) {
        for (i in 1:length(private$leg.times.rel)) {
          leg.name <- colnames(private$leg.times.rel)[i]
          print(paste0(leg.name, ": ", private$leg.times.rel[i]))
        }
      }
      print(paste0("Stop: ", private$stop.time.rel))
    },
      
    # Reset chronometer
    reset = function() {
      private$start.time.abs = NULL
      private$stop.time.abs = NULL
      private$leg.times.abs = matrix(nrow = 1, ncol = 0)
      private$stop.time.rel = NULL
      private$leg.times.rel = matrix(nrow = 1, ncol = 0)
    },
  
    # Start chronometer
    start = function() {
      if (!is.null(private$start.time.abs)) {
        stop("Chronometer already started")
      }
      private$start.time.abs <- lubridate::now()
    },
    
    # Stop chronometer
    stop = function(leg.name = NULL) {
      if (is.null(private$start.time.abs)) {
        stop("Chronometer not started yet")
      }
      if (!is.null(private$stop.time.rel)) {
        stop("Chronometer already stopped")
      }
      
      private$stop.time.abs <- lubridate::now()
      private$stop.time.rel <- 
        private$get.time.diff(private$start.time.abs, private$stop.time.abs)
      
      # If legs were being recorded, also record the final leg time
      if (ncol(private$leg.times.rel) != 0) {
        private$append.leg.time.abs(private$stop.time.abs, leg.name)
        private$update.leg.times.rel()
      }
    },
    
    # Record a leg time
    record.leg = function (leg.name = NULL) {
      if (is.null(private$start.time.abs)) {
        stop("Chronometer not started yet")
      }
      if (!is.null(private$stop.time.rel)) {
        stop("Chronometer already stopped")
      }
      
      # Compute and record leg time
      private$append.leg.time.abs(lubridate::now(), leg.name)
      private$update.leg.times.rel()
    },
    
    # Return all of the variables of the chronometer
    get.chronometer.data = function() {
      list(private$start.time.abs,
           private$stop.time.abs,
           private$leg.times.abs,
           private$stop.time.rel,
           private$leg.times.rel
           )
    },
    
    # Return the total time span between start and stop
    get.total.time = function(force = FALSE) {
      if (is.null(private$stop.time.rel)) {
        if (force) {
          return(-Inf)
        } else {
          stop("Chronometer not stopped yet")
        }
      }
      private$stop.time.rel
    },
    
    # Return the total time spent per leg category (based on leg names)
    accumulate.legs = function(force = FALSE) {
      if (ncol(private$leg.times.rel) == 0) {
        if (force) {
          return(rep(-Inf, 4))
        } else {
          stop("No leg times recorded yet.")
        }
      }
      categories <- unique(colnames(private$leg.times.rel))
      rtrn <- matrix(nrow = 1, ncol = length(categories))
      colnames(rtrn) <- categories
      for (cat in categories) {
        cat.idx <- which(categories == cat)
        idxs.col.to.sum <- which(colnames(private$leg.times.rel) == cat)
        rtrn[1, cat.idx] <- sum(private$leg.times.rel[1, idxs.col.to.sum])
      }
      rtrn
    }
  ),
  
  private = list(
    
    # Variables
    start.time.abs = NULL,
    stop.time.abs = NULL,
    leg.times.abs = matrix(nrow = 1, ncol = 0),
    stop.time.rel = NULL,
    leg.times.rel = matrix(nrow = 1, ncol = 0),
    
    # Function to get time differences in seconds
    get.time.diff = function(t1, t2) {
      as.numeric(abs(difftime(t1, t2, units = "secs")))
    },
    
    # Function to append a leg time to the matrix of leg times.
    append.leg.time.abs = function(time, leg.name) {
      old.col.names <- colnames(private$leg.times.abs)
      if (is.null(leg.name)) {
        leg.name <- sprintf("Leg %s", ncol(private$leg.times.abs) + 1)
      }
      private$leg.times.abs <- cbind(private$leg.times.abs, time)
      colnames(private$leg.times.abs) <- c(old.col.names, leg.name)
    },
    
    # Function to update the relative leg times based on the absolute leg times
    update.leg.times.rel = function () {
      time.vct <- c(private$start.time.abs, c(private$leg.times.abs))
      private$leg.times.rel <- matrix(nrow = 1, ncol = ncol(private$leg.times.abs))
      for (i in 1:ncol(private$leg.times.abs)) {
        private$leg.times.rel[1, i] <- private$get.time.diff(time.vct[i], time.vct[i+1])
      }
      colnames(private$leg.times.rel) <- colnames(private$leg.times.abs)
    }
  )
)
