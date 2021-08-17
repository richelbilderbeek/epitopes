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
#'     data frame)
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
#' Each of these groups may be used either at the peptide or the protein level -
#' which does not mean it _should_ be. The **protr** documentation provides
#' the following warning: "*Users need to intelligently evaluate the underlying*
#' *details of the descriptors provided, instead of using protr with their data*
#' *blindly, especially for the descriptor types with more flexibility. It*
#' *would be wise for the users to use some negative and positive control*
#' *comparisons where relevant to help guide interpretation of the results.*".
#' Users should be savvy when choosing which features to use for epitope
#' prediction, and the choice should ideally be guided by domain expertise.
#' **IMPORTANT**: if the splitting level intended to be used for modelling is
#' "peptide" (this can be checked in `peptides.list$splits.attrs$split_level`)
#' then protein-level features should be avoided, as they can result in data
#' leakage across splits and contaminate performance calculations.
#'
#' @section **Feature Lists**:
#' Input lists `peptide.features` and `protein.features` are used to define
#' which features are calculated at either level. These input parameters must be
#' named list vectors, where each element is itself a list with the name of the
#' feature group the user wants to compute, and elements corresponding to the
#' desired parameters to be passed down to the feature calculation routines.
#' For more information on the parameters needed for any feature group, check
#' `?protr::extractXYZ`, replacing `XYZ` by the group abbreviation (see
#' **Details** or the documentation of the **protr** package for the list of
#' available feature groups).
#'
#' For instance, suppose the user wants to calculate the
#' _Amino Acid Composition_, "AAC" (no parameters required) and the
#' _Pseudo Amino Acid Composition_, "PAAC" (up to four optional parameters -
#' assume that the user wants to set parameter `lambda = 50`), at the peptide
#' level. In this case, `peptide.features` should be:
#'
#' ```
#' peptide.features = list(AAC = list(), PAAC = list(lambda = 50))
#' ```
#'
#' @param peptides.list list object returned by [make_data_splits()].
#' @param peptide.features,protein.features lists of features to be calculated
#' at the peptide/protein level. See **Feature Lists** for details.
#' @param ncpus positive integer, number of cores to use.
#'
#' @author Felipe Campelo (\email{f.campelo@@aston.ac.uk})
#'
#' @export
#'
#' @return Updated `peptides.list` object, with features calculated for element
#' `peptides.list$df`.
#'
#'
calc_features <- function(peptides.list,
                          peptide.features = list(),
                          protein.features = list(),
                          ncpus = 1){
  # ========================================================================== #
  # Sanity checks and initial definitions


  # ========================================================================== #




}
