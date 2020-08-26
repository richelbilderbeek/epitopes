#' Calculate sequence entropy of peptides
#'
#' This function is used to calculate the sequence entropies
#' from a vector of peptides. Entropy = -sum(p_k * log2(p_k))
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

calc_sequence_entropy <- function(input){
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

  # Entropy = -SUM(p(e)log2[p(e)])
  cat("\nCalculating sequence entropy...")
  entropy <- sapply(pepvec,
                    function(x){
                      # table the frequency of each aa in a window sequence
                      aa <- table(strsplit(x, split = "")[[1]])
                      pk <- aa / sum(aa)
                      H  <- pk * log2(pk)
                      H[is.nan(H)] <- 0
                      return(-sum(H))
                    })

  return(cbind(input, entropy_of_seq = entropy))
}
