#' Calculate frequencies of N-peptides
#'
#' This function is used to calculate the frequencies of N-peptide subsequences
#' from a vector of peptides.
#'
#' @param input either a vector of peptides or a data frame with a variable
#' called "window_seq" containing the peptides.
#' @param N the length of the subsequences to consider.
#'
#' @return Data frame containing the calculated AA percentages. If `input` is
#' a`data.frame` the original input is returned with the AA percentages
#' appended as columns.
#'
#' @author Felipe Campelo (\email{f.campelo@@aston.ac.uk});
#'         Jodie Ashford (\email{ashfojsm@@aston.ac.uk})
#'
#' @export
#'
calc_Npeptide_composition <- function(input, N = 2){

  # ========================================================================== #
  # Sanity checks and initial definitions
  assertthat::assert_that(assertthat::is.count(N))

  if(is.data.frame(input)){
    assertthat::assert_that("window_seq" %in% names(input),
                            is.character(input$window_seq))
    pepvec <- input$window_seq
  } else {
    assertthat::assert_that(is.vector(input),
                            is.character(input))
    pepvec <- input
    input  <- data.frame(window_seq = pepvec)
  }

  aa_codes <- get_aa_codes()
  isvalid <- sapply(pepvec,
                    function(x){
                      all(strsplit(x, split = "")[[1]] %in% aa_codes)
                    })

  if(any(!isvalid)){
    stop("Unrecognised aminoacid code in entry(ies): ",
         paste(which(!isvalid), collapse = ", "))

  }
  # ========================================================================== #

  # Creates a list of tables contaning the dipeptide percentages for each window sequence in df
  Npeptides <- lapply(X = pepvec,
                      FUN = function(x) {
                        start <- 1:(nchar(x) - N + 1)
                        stop  <- start + N - 1
                        table(substring(x, start, stop)) / nchar(x)
                      })

  # binds the rows of all dipeptide sequence tables
  Npeptides_full <- do.call(what = dplyr::bind_rows, args = Npeptides)
  Npeptides_full[is.na(Npeptides_full)] <- 0

  # Generate any columns that may be missing
  myarg <- vector("list", N)
  for (i in 1:N) myarg[[i]] <- aa_codes
  myarg   <- do.call(expand.grid, myarg)
  toadd   <- apply(myarg, MARGIN = 1,
                   paste, collapse = "")
  toadd   <- toadd[!(toadd %in% names(Npeptides_full))]
  newcols <- as.data.frame(matrix(0, nrow = nrow(input), length(toadd)))

  names(newcols) <- toadd
  Npeptides_full <- cbind(Npeptides_full, newcols)

  # change column names
  names(Npeptides_full) <- paste0("perc_of_", names(Npeptides_full))
  Npeptides_full <- Npeptides_full[, order(names(Npeptides_full))]

  return(cbind(input, Npeptides_full))
}
