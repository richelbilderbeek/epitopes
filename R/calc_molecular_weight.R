#' Calculate molecular weight of peptide sequences
#'
#' This function is used to calculate the molecular weights from a vector
#' of peptides.
#'
#' @param input either a vector of peptides or a data frame with a variable
#' called "window_seq" containing the peptides.
#'
#' @return Data frame containing the calculated molecular weights. If `input` is
#' a`data.frame` the original input is returned with the molecular weights
#' appended as columns.
#'
#' @author Felipe Campelo (\email{f.campelo@@aston.ac.uk});
#'         Jodie Ashford (\email{ashfojsm@@aston.ac.uk})
#'
#' @export
#'
calc_molecular_weight <- function(input){

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

  aa_prop <- readRDS(system.file("extdata", "amino_acid_propensity.rds",
                                 package = "epitopes"))

  total_mw <- sapply(pepvec,
                     FUN = function(x){
                       w <- 0
                       x <- strsplit(x, "")[[1]]
                       for(i in seq_along(x)){
                         w <- w + aa_prop[aa_prop$One_letter_code == x[i],
                                          "Amino_acid_molecular_weight"]
                       }
                       return(w)
                     })

  input$total_molecular_weight <- total_mw

  return(input)
}
