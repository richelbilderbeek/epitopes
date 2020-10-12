# Little field check function
nullcheck <- function(x) { ifelse(is.null(x), yes = NA, no = x) }


# ======================================================================
# Progress bar function
mypb <- function(i, max_i, t0, npos){
  nb <- max(1, ceiling(max_i / npos))
  if (i == 0){
    pbstr <- paste0("\n  |", paste(rep("_", npos), collapse = ""), "|")
    cat(pbstr, "0% processed. Elapsed time: 0s")
  } else if (i >= (max_i - 1)) {
    pbstr <- paste(rep(">", times = npos), collapse = "")
    td <- Sys.time() - t0
    perc_done <- 100
    cat(sprintf("\r  |%s|%d%% processed. Elapsed time: %2.1f %s",
                pbstr, perc_done, as.numeric(td), attr(td, "units")))
  } else if (!(i %% nb)) {
    nn <- i / nb
    pbstr <- paste(rep(c(">", "_"), times = c(nn, npos - nn)),
                   collapse = "")
    td <- Sys.time() - t0
    perc_done <- round(100 * i / max_i, digits = 0)
    cat(sprintf("\r  |%s|%d%% processed. Elapsed time: %2.1f %s",
                pbstr, perc_done, as.numeric(td), attr(td, "units")))
  }
}


ismc <- function(cl){
  # Check if multilpe cores are to be used. cl is either a cluster object or a
  # numeric integer
  return(length(cl) > 1 || is.numeric(cl) & cl > 1)
}
