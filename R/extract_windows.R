# Function to extract windows from a given row of df
# (local scope only - not exported)
extract_windows <- function(x, window_size, step_size, window_exp){

  # Initialise dataframe
  wdf <- data.frame(window_seq = rep(NA_character_, x$epitope_len + 2),
                    window_exp = rep(NA_character_, x$epitope_len + 2),
                    Class      = x$Class,
                    epitope_id = x$epitope_id,
                    protein_id = x$protein_id,
                    taxid      = x$protein_taxid,
                    host_id    = x$host_id,
                    org_id     = x$org_id,
                    file_id    = x$file_id,
                    stringsAsFactors = FALSE)

  if (window_exp <= 0) wdf$window_exp <- NULL

  # Initial position for window (and expanded window)
  i1   <- max(1, x$epitope_start - floor(window_size / 2))
  i2   <- min(x$protein_len, i1 + window_size - 1)
  i1e  <- max(1, i1 - window_exp)
  i2e  <- min(x$protein_len, i2 + window_exp)
  j    <- 1
  stop <- FALSE
  while(!stop){
    wdf$window_seq[j] <- substr(x$protein_seq, i1, i2)
    if(window_exp > 0) wdf$window_exp[j] <- substr(x$protein_seq, i1e, i2e)
    i1  <- i1 + step_size
    i2  <- i2 + step_size
    i1e <- i1e + step_size
    i2e <- min(x$protein_len, i2e + step_size)
    j   <- j + 1
    if(i2 > min(x$protein_len,  x$epitope_stop + ceiling(window_size / 2))) {
      stop <- TRUE
    }
  }
  return(wdf[!is.na(wdf$window_seq), ])
}
