# Function to extract windows from a given row of df
# (local scope only - not exported)
extract_windows <- function(x, window_size, step_size){
  # Initialise dataframe
  wdf <- data.frame(window_seq = rep(NA_character_, x$epitope_len + 2),
                    Class      = x$Class,
                    epitope_id = x$epitope_id,
                    protein_id = x$protein_id,
                    stringsAsFactors = FALSE)

  # Initial position for window
  i1   <- max(1, x$epitope_start - floor(window_size / 2))
  i2   <- min(x$protein_len, i1 + window_size - 1)
  stop <- FALSE
  j    <- 1
  while(!stop){
    wdf$window_seq[j] <- substr(x$protein_seq, i1, i2)
    i1 <- i1 + step_size
    i2 <- i2 + step_size
    j  <- j + 1
    if(i2 > min(x$protein_len,  x$epitope_stop + ceiling(window_size / 2))) {
      stop <- TRUE
    }
  }
  return(wdf[!is.na(wdf$window_seq), ])
}
