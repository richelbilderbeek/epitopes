#' Smoothing of prediction vectors
#'
#' This function runs a sliding window smoother on prediction vectors.
#'
#' @param x vector of predictions (must be numeric)
#' @param window_size size of the sliding window
#' @param type type of smoothing function to use. Accepts "mode" (for majority 
#'        vote) or "mean".
#'
#' @return vector of smoothed prediction values.
#'
#' @author Felipe Campelo (\email{f.campelo@@aston.ac.uk})
#'
#' @export
#'
smooth_predictions <- function(x, window_size, type = "mode"){
  
  valid_types <- c("mean", "mode")
  assertthat::assert_that(is.numeric(x), 
                          assertthat::is.count(window_size),
                          length(x) >= window_size,
                          is.character(type),
                          length(type) == 1,
                          type %in% valid_types)
  
  xlen <- length(x)
  y    <- x
  for (i in 1:xlen){
    # Define range covered by the window
    i1 <- max(1, min(i - floor(window_size / 2),
                     xlen - window_size + 1))
    i2 <- min(xlen, i1 + window_size - 1)
    
    # Perform smoothing
    if (type == "mean"){
      y[i] <- mean(x[i1:i2])
    } else if (type == "mode"){
      uniqx <- unique(x[i1:i2])
      y[i] <-uniqx[which.max(tabulate(match(x[i1:i2], uniqx)))]
    } else stop("Unrecognised 'type' in smooth_predictions")
  }
  
  return(y)
}