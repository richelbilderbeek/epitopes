#' Calculate aminoacid frequencies
#'
#' This function is used to calculate aminoacid frequencies from a vector of
#' peptides.
#'
#' @param input either a vector of peptides or a data frame with a variable
#' called "window_seq" containing the peptides.
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
amino_acid_percentages <- function(input){

  # ========================================================================== #
  # Sanity checks and initial definitions
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

  # Creates a table of aa percentages for each sequence
  perc_values <- lapply(X = pepvec,
                        FUN = function(x){
                          table(unlist(strsplit(x, split = ""))) / nchar(x)[[1]]
                        })

  # binds all percentage tables
  aa_percent <- do.call(what = dplyr::bind_rows, args = perc_values)
  aa_percent[is.na(aa_percent)] <- 0
  names(aa_percent) <- paste0("percent_", names(aa_percent))

  # adds the aa percentage columns to the data frame
  input <- cbind(input, aa_percent)

  return(input)
}
