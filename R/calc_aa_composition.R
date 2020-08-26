#' Calculate qualitative aminoacid composition of peptides
#'
#' This function is used to calculate the percent composition of peptides, in
#' terms of nine AA types: "Tiny", "Small", "Aliphatic", "Aromatic",
#' "NonPolar", "Polar", "Charged", "Basic" and "Acidic".
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
calc_aa_composition <- function(input){
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

  # Compute the aa composition of each sequence
  cat("\nCalculating AA composition...")
  suppressMessages({
    tmp <- as.matrix(dplyr::bind_cols(lapply(Peptides::aaComp(pepvec),
                                             as.data.frame)))
  })
  # remove all counts leaving just % behind and transpose matrix
  tmp <- t(tmp[, grep("Mole%", colnames(tmp))]) / 100
  colnames(tmp) <- paste0("perc_", colnames(tmp))
  rownames(tmp) <- NULL

  return(cbind(input, tmp))
}
