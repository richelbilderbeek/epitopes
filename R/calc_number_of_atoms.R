#' Calculate number of each type of atom in peptide sequences
#'
#' This function is used to calculate the number of specific atoms from a vector
#' of peptides.
#'
#' @param input either a vector of peptides or a data frame with a variable
#' called "window_seq" containing the peptides.
#'
#' @return Data frame containing the calculated numbers of atoms. If `input` is
#' a`data.frame` the original input is returned with the numbers of atoms
#' appended as columns.
#'
#' @author Felipe Campelo (\email{f.campelo@@aston.ac.uk});
#'         Jodie Ashford (\email{ashfojsm@@aston.ac.uk})
#'
#' @importFrom dplyr %>%
#' @export
#'
calc_number_of_atoms <- function(input){

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

  elements = c('Number_of_carbon_atoms_in_aa',
               'Number_of_hydrogen_atoms_in_aa',
               'Number_of_nitrogen_atoms_in_aa',
               'Number_of_oxygen_atoms_in_aa',
               'Number_of_sulphur_atoms_in_aa')

  # adds the total number of atoms of each element for each sequence to the
  # data frame
  for(element in elements){
    atom_counts <- sapply(X = pepvec,
                          FUN = function(x){
                            # Get AA counts
                            aa_counts <- table(strsplit(x, "")[[1]])

                            a_count   <- aa_prop %>%
                              # for every aa in the table (of the sequence split string)
                              dplyr::filter(.data$One_letter_code %in% names(aa_counts)) %>%
                              # retrieve the number of atoms for that aa from the propensity table
                              dplyr::select(.data$One_letter_code, !!!element) %>%
                              # and the number of that aa in the sequence
                              dplyr::mutate(Count = aa_counts)
                            sum(a_count[, element] * a_count$Count)
                          })

    input[[substr(element, 0, (nchar(element)-6))]] <- atom_counts
  }

  return(input)

}
