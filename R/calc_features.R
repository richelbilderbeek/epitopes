#' Calculate features for epitope prediction
#'
#' This function is used to calculate several distinct families of features for
#' epitope prediction.
#'
#' Two major types of features can be calculated:
#' \itemize{
#'     \item _local features_, which are calculated based on the local
#'     neighbourhood of each AA position (column *Info_window* of the windowed
#'     data frame).
#'     \item _global features_, which are calculated using the full sequence
#'     of the protein (listed in column *Info_protein_id*  of the windowed
#'     data frame).
#' }
#'
#' All features are calculated based on the implementations available in package
#' [**protr**](https://cran.r-project.org/package=protr). As of August 2021,
#' the following groups of features were available (see the
#' [package vignette](https://cran.r-project.org/web/packages/protr/vignettes/protr.html#3_package_overview)
#' for details):
#' \itemize{
#'    \item "AAC" - Amino acid composition
#'    \item "DC"  - Dipeptide composition
#'    \item "TC"  - Tripeptide composition
#'    \item  Autocorrelation:
#'    \itemize{
#'        \item "MoreauBroto" - Normalized Moreau-Broto autocorrelation
#'        \item "Moran" - Moran autocorrelation
#'        \item "Geary" - Geary autocorrelation
#'    }
#'    \item  CTD descriptors:
#'    \itemize{
#'        \item "CTDC" - Composition
#'        \item "CTDT" - Transition
#'        \item "CTDD" - Distribution
#'    }
#'    \item "CTriad" - Conjoint triad descriptors
#'    \item  Quasi-sequence-order descriptors:
#'    \itemize{
#'        \item "SOCN" - Sequence-order-coupling number
#'        \item "QSO" - Quasi-sequence-order descriptors
#'    }
#'    \item  Pseudo-amino acid composition:
#'    \itemize{
#'        \item "PAAC" - Pseudo-amino acid composition (PseAAC)
#'        \item "APAAC" - Amphiphilic pseudo-amino acid composition (APseAAC)
#'    }
#'    \item  Profile-based descriptors:
#'    \itemize{
#'        \item "PSSM"
#'        \item "PSSMAcc"
#'        \item "PSSMFeature"
#'    }
#'    \item Proteochemometric Modeling descriptors:
#'    \itemize{
#'        \item "ScalesGap" - Scales-based descriptors derived by Principal
#'        Components Analysis
#'        \item "ProtFPGap" - Scales-based descriptors derived by amino acid
#'        properties from AAindex (a.k.a. Protein Fingerprint)
#'        \item "DescScales" - Scales-based descriptors derived by 20+ classes
#'        of 2D and 3D molecular descriptors (Topological, WHIM, VHSE, etc.)
#'        \item "FAScales" - Scales-based descriptors derived by Factor Analysis
#'        \item "MDSScales" - Scales-based descriptors derived by
#'        Multidimensional Scaling
#'        \item "BLOSUM" - BLOSUM and PAM matrix-derived descriptors
#'    }
#' }
#'
#' **NOTE**: in all feature groups except those ending in "Gap", invalid AA codes
#' (B, J, O, U, X, Y) are removed from the strings prior to feature calculation.
#' In those ending with "Gap" these codes are replaced by a gap indicator, "-".
#'
#' Each feature group may be used either at the local or global level -
#' which does not mean it _should_ be. The **protr** documentation provides
#' the following warning: "*Users need to intelligently evaluate the underlying*
#' *details of the descriptors provided, instead of using protr with their data*
#' *blindly, especially for the descriptor types with more flexibility. It*
#' *would be wise for the users to use some negative and positive control*
#' *comparisons where relevant to help guide interpretation of the results.*".
#' Users should therefore be savvy when choosing which features to use for epitope
#' prediction, and the choice should ideally be guided by domain expertise.
#' Certain feature groups may not make sense as local features,
#' as the (usually very short) local substrings will not allow the features to
#' be informative (e.g., tripeptide composition, "TC"); or may require
#' specific overriding of standard parameters (e.g., parameter
#' `lambda` for "PAAC" must be smaller than the length of the shortest local
#' string).
#'
#' **IMPORTANT**: if the splitting level intended to be used for modelling is
#' "peptide" (this can be checked in `peptides.list$splits.attrs$split_level`)
#' then protein-level features should be avoided, as their use can result in
#' data leakage across splits and contaminate performance calculations.
#'
#' @section **Feature Lists**:
#' Input lists `local.features` and `global.features` are used to define
#' which features are calculated at either level. These input parameters must be
#' named lists, where each element is itself a list named after the
#' feature group the user wants to compute, containing elements corresponding to
#' the desired parameters to be passed down to the feature calculation routines.
#' For more information on the parameters needed for any feature group, check
#' `?protr::extractXYZ`, replacing `XYZ` by the group abbreviation (see
#' **Details** or the documentation of the **protr** package for the list of
#' available feature groups).
#'
#' For instance, suppose the user wants to calculate the local features
#' _Amino Acid Composition_, "AAC" (no parameters required) and
#' _Pseudo Amino Acid Composition_, "PAAC" (assume that the user wants to set
#' override parameter `lambda` with `lambda = 5`). In this case,
#' `local.features` should be:
#'
#' ```
#' local.features = list(AAC = list(), PAAC = list(lambda = 5))
#' ```
#'
#' @param peptides.list list object returned by [make_data_splits()].
#' @param local.features list of features to be calculated
#' at the local neighbourhood (`peptides.list$df$Info_window`) level. Notice
#' that invalid AA codes (B, J, O, U, X and Y) are removed prior to the
#' calculation of local features. See **Feature Lists** for details.
#' @param global.features lists of features to be calculated
#' at the peptide/protein level. See **Feature Lists** for details.
#' @param ncpus positive integer, number of cores to use.
#' @param overwrite logical: if a feature group already exists in
#' `peptides.list$df` (local) or `peptides.list$proteins` (global), should it be
#' overwritten? If it already exists and `overwrite == FALSE` then a message is
#' printed and calculation is skipped for that group.
#'
#' @return Updated `peptides.list` object, with local features added as columns
#' to `peptides.list$df`, and global features added as columns to
#' `peptides.list$proteins`.
#'
#' @author Felipe Campelo (\email{f.campelo@@aston.ac.uk})
#'
#' @export
#'
#' @importFrom dplyr %>%
#' @importFrom rlang .data
#'
#'
calc_features <- function(peptides.list,
                          local.features = list(),
                          global.features = list(),
                          ncpus = 1,
                          overwrite = FALSE){
  # ========================================================================== #
  # Sanity checks and initial definitions
  assertthat::assert_that(is.list(peptides.list),
                          all(c("df", "proteins") %in% names(peptides.list)),
                          is.list(local.features),
                          is.list(global.features),
                          assertthat::is.count(ncpus))
  # ========================================================================== #

  # Calculate local features
  if(length(local.features) > 0) {
    for (i in seq_along(local.features)){
      y <- call_protr(SEQs      = peptides.list$df$Info_window,
                      feat.list = local.features[i],
                      txt.opts  = c("local", "df"),
                      dfnames   = names(peptides.list$df),
                      ncpus     = ncpus)

      if(is.data.frame(y)) {
        peptides.list$df <- cbind(peptides.list$df, y)
      }
    }
    class(peptides.list) <- unique(c(class(peptides.list), "local.features"))
  }

  # Calculate global features
  if(length(global.features) > 0) {
    for (i in seq_along(global.features)){
      y <- call_protr(SEQs      = peptides.list$proteins$TSeq_sequence,
                      feat.list = global.features[i],
                      txt.opts  = c("global", "proteins"),
                      dfnames   = names(peptides.list$proteins),
                      ncpus     = ncpus,
                      overwrite = overwrite)

      if(is.data.frame(y)) {
        peptides.list$proteins <- cbind(peptides.list$proteins, y)
      }
    }
    class(peptides.list) <- unique(c(class(peptides.list), "global.features"))
  }

  return(peptides.list)
}
