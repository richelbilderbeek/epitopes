#' Calculate a set of mean AA descriptors
#'
#' This function is used to calculate 66 aminoacid descriptors based on
#' function [Peptides::aaDescriptors()]. These descriptors apply to each
#' individual aminoacid in a peptide, and are summarised for the peptide using a
#' simple mean. The features returned are provided in the `Details`.
#'
#'   \itemize{
#'      \item\code{crucianiProperties}
#'      \itemize{
#'          \item PP1: Polarity
#'          \item PP2: Hydrophobicity
#'          \item PP3: H-bonding
#'      }
#'      \item\code{kideraFactors} (the four first features are essentially
#'      pure physical properties; the remaining six are linear combinations of
#'      several properties):
#'      \itemize{
#'          \item KF1: Helix/bend preference
#'          \item KF2: Side-chain size
#'          \item KF3: Extended structure preference
#'          \item KF4: Hydrophobicity
#'          \item KF5 - KF10: linear combinations of other characteristics
#'      }
#'      \item\code{zScales} (based on physicochemical properties of the AAs
#'      including NMR data and thin-layer chromatography data):
#'      \itemize{
#'          \item Z1: Lipophilicity
#'          \item Z2: Steric properties (Steric bulk/Polarizability)
#'          \item Z3: Electronic properties (Polarity / Charge)
#'          \item Z4-Z5: relate electronegativity, heat of formation,
#'          electrophilicity and hardness.
#'      }
#'      \item\code{FASGAI} (based on physicochemical properties of the AAs
#'      including NMR data and thin-layer chromatography data):
#'      \itemize{
#'          \item F1: Hydrophobicity index
#'          \item F2: Alpha and turn propensities
#'          \item F3: Bulky properties
#'          \item F4: Compositional characteristic index
#'          \item F5: Local flexibility
#'          \item F6: Electronic properties
#'      }
#'      \item\code{tScales} (based on 67 common topological descriptors of
#'      amino acids. These topological descriptors are based on the connectivity
#'      table of amino acids alone, and to not explicitly consider 3D properties
#'      of each structure):
#'      \itemize{
#'          \item T1 - T5
#'      }
#'      \item\code{VHSE scales} (principal component score Vectors of
#'      Hydrophobic, Steric, and Electronic properties. Derived from PCA on
#'      independent families of 18 hydrophobic properties, 17 steric properties,
#'      and 15 electronic properties):
#'      \itemize{
#'          \item VHSE1 - VHSE2: Hydrophobic properties
#'          \item VHSE3 - VHSE4: Steric properties
#'          \item VHSE5 - VHSE8: Electronic properties
#'      }
#'      \item\code{ProtFP descriptors}: (constructed from a large initial
#'      selection of indices obtained from the AAindex database for all 20
#'      naturally occurring amino acids.):
#'      \itemize{
#'          \item ProtFP1-ProtFP8
#'      }
#'      \item\code{stScales} (proposed by Yang et al.'2010, taking 827
#'      properties into account which are mainly constitutional, topological,
#'      geometrical, hydrophobic, electronic, and steric properties of AAs):
#'      \itemize{
#'          \item ST1 - ST8
#'      }
#'      \item\code{BLOSUM indices} (derived of physicochemical properties that
#'      have been subjected to a VARIMAX analyses and an alignment matrix of the
#'      20 natural AAs using the BLOSUM62 matrix):
#'      \itemize{
#'          \item BLOSUM1 - BLOSUM10
#'      }
#'      \item\code{MS-WHIM scores} (derived from 36 electrostatic potential
#'      properties derived from the three-dimensional structure of the 20
#'      natural amino acids):
#'      \itemize{
#'          \item MSWHIM1 - MSWHIM3
#'      }
#' }
#'
#' See the documentation for [Peptides::aaDescriptors()] for further information
#' and references.
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
#' @importFrom data.table rbindlist
#' @importFrom parallel mclapply
#' @export
#'
calc_aa_descriptors <- function(input, ncores = 1){
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

  # Prepare feature names
  col_names <- colnames(Peptides::aaDescriptors("K"))

  # Calculate features for all windows
  cat("\n Calculating AA descriptors:")
  X <- data.table::rbindlist(pbmcapply::pbmclapply(pepvec,
                                                   function(x){
                                                     x <- strsplit(x, split = "")[[1]]
                                                     indfeat <- t(rowMeans(sapply(x,Peptides::aaDescriptors)))
                                                     colnames(indfeat) <- col_names
                                                     return(as.data.frame(indfeat))
                                                   },
                                                   mc.preschedule = TRUE,
                                                   mc.cores = ncores))

  return(cbind(input, X))
}
