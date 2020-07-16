#' Calculate cojoint triads
#'
#' This function is used to calculate the cojoint triads for a vector of
#' peptides.
#'
#' @param input either a vector of peptides or a data frame with a variable
#' called "window_seq" containing the peptides.
#' @param ncores number of cores to use when calculating the features.
#'
#' @return Data frame containing the calculated AA percentages. If `input` is
#' a`data.frame` the original input is returned with the AA percentages
#' appended as columns.
#'
#' @author Felipe Campelo (\email{f.campelo@@aston.ac.uk});
#'         Jodie Ashford (\email{ashfojsm@@aston.ac.uk})
#'
#' @importFrom parallel mclapply
#' @export
#'
calc_cojoint_triads <- function(input, ncores = 1){
  # ========================================================================== #
  # Sanity checks and initial definitions
  assertthat::assert_that(assertthat::is.count(ncores))

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

  aa_prop <- readRDS(system.file("extdata", "amino_acid_propensity.rds",
                                 package = "epitopes"))

  # creates a list of conjoint triads for each sequence
  triad_list <- mclapply(pepvec,
                         function(x){
                           # Calculate CT representation
                           xpl <- strsplit(x, split = "")[[1]]
                           idx <- sapply(xpl,
                                         function(y) {
                                           which(aa_prop$One_letter_code == y)
                                         })
                           ct_rep <- paste(aa_prop$CT_group[idx], collapse = "")

                           # Extract triads
                           ct_list <- substring(ct_rep,
                                                first = 1:(nchar(ct_rep) - 2),
                                                last  = 3:nchar(ct_rep))

                           # Return frequencies of eah triad
                           table(ct_list) / length(ct_list)
                         },
                         mc.cores = ncores)


  # Binds the rows of all conjoint triad tables
  triad_full <- do.call(what = dplyr::bind_rows, args = triad_list)
  triad_full[is.na(triad_full)] <- 0

  # Generate any columns that may be missing
  myarg   <- expand.grid(0:6, 0:6, 0:6)
  toadd   <- apply(myarg, MARGIN = 1,
                   paste, collapse = "")
  toadd   <- toadd[!(toadd %in% names(triad_full))]
  newcols <- as.data.frame(matrix(0, nrow = nrow(input), length(toadd)))

  triad_full <- cbind(triad_full, newcols)
  triad_full <- triad_full[, order(names(triad_full))]

  # Changes the names of all columns
  names(triad_full) <- paste0("perc_of_", names(triad_full))


  return(cbind(input, triad_full))
}
