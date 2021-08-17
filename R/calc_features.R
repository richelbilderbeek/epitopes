#' Calculate features for epitope prediction
#'
#' This function is used to calculate several distinct families of features for
#' epitope prediction. Two major types of features can be calculated:
#'
#' \itemize{
#'     \item _local features_, which are calculated based on the local
#'     neighbourhood of each AA position (column *Info_window* of the windowed
#'     data frame).
#'     \item _global features_, which are calculated using the full sequence
#'     of the protein (listed in column *Info_protein_id*  of the windowed
#'     data frame)
#' }
#'
#' All features are calculated based on the implementations available in package
#' [**protr**](https://cran.r-project.org/package=protr).
#'
#' @param splits.list list object returned by [make_data_splits()].
#'
#' @author Felipe Campelo (\email{f.campelo@@aston.ac.uk})
#'
#' @export
#'
calc_features <- function(splits.list,
                          peptide.features = list(),
                          protein.features = list(),
                          ncpus = 1){
  # ========================================================================== #
  # Sanity checks and initial definitions


  # ========================================================================== #





}
