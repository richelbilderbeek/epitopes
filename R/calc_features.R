#' Calculate statistical and physicochemical features for peptides
#'
#' This function is used to calculate several distinct families of features for
#' a vector of peptides.
#'
#' To get a list of all routines for calculating the
#' individual types of features use `get_feature_functions()`. Some of the
#' features are calculated based on an aminoacide propensity scale, which was
#' originally obtained from _E.L. Ulrich et al., "BioMagResBank". Nucleic Acids
#' Research 36, D402-D408 (2008) DOI: 10.1093/nar/gkm957_. The data was
#' downloaded from [http://www.bmrb.wisc.edu/ref_info/aadata.dat](http://www.bmrb.wisc.edu/ref_info/aadata.dat)
#'
#' @param input either a vector of peptides or a data frame with a variable
#' called "window_seq" containing the peptides.
#' @param max.N maximum length of N-peptide frequency features to be calculated.
#'              See `calc_Npeptide_composition()` for details.
#' @param ncores number of cores to use for calculating some of the features.
#'
#' @return Data frame containing the calculated features. If `input` is
#' a`data.frame` the original input is returned with the features
#' appended as columns.
#'
#' @author Felipe Campelo (\email{f.campelo@@aston.ac.uk});
#'         Jodie Ashford (\email{ashfojsm@@aston.ac.uk})
#'
#' @importFrom dplyr %>%
#' @export
#'
calc_features <- function(input,
                          max.N = 2,
                          ncores = 1){
  # ========================================================================== #
  # Sanity checks and initial definitions
  assertthat::assert_that(assertthat::is.count(ncores),
                          assertthat::is.count(max.N))

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
    cat("Unrecognised aminoacid code in entry(ies): ",
         paste(which(!isvalid), collapse = ", "))
    cat("\n(Probably due to degenerate aminoacid codes)")
    cat("\nRemoving these entries before proceeding...")
    input <- input[-which(!isvalid), ]
  }
  # ========================================================================== #

  # Calculate features
  cat("\nCalculating features:")
  df <- input %>%
    calc_aa_composition() %>%
    calc_aa_descriptors(ncores = ncores) %>%
    calc_molecular_weight() %>%
    calc_number_of_atoms() %>%
    calc_sequence_entropy() %>%
    calc_cojoint_triads(ncores = ncores)

  # Add Npeptide percentages
  for (i in 1:max.N){
    df <- calc_Npeptide_composition(df, N = i, ncores = ncores)
  }

  return(df[order(df$protein_id,
                  df$center_pos), ])
}
