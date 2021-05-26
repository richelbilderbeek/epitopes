#' Extract peptides and windowed representation for epitope prediction.
#'
#' Extract relevant peptides from the consolidated protein-epitope data
#' (generated using [consolidate_data()]) and build column with neightbourhood
#' of each labeled aminoacid residue.
#'
#' @param df data frame of consolidated protein-epitope data, returned by
#' [consolidate_data()].
#' @param min_peptide positive integer, shortest peptide to be considered
#' @param max_peptide positive integer, longest peptide to be considered
#' @param window_size positive integer, size of the local neighbourhood to be
#' considered.
#' @param save_folder path to folder for saving the results.
#'
#' @return List containing two data frames:
#'
#' \itemize{
#'    \item A data frame containing the labeled positions of `df` (filtered by
#' *min_peptide* and *max_peptide*) with a new column expressing the local
#' neighbourhood of each position (to be used for feature calculation)
#'    \item A data frame containing the individual peptides and classes
#' }
#'
#' @author Felipe Campelo (\email{f.campelo@@aston.ac.uk})
#'
#' @export
#'
#' @importFrom dplyr %>%
#' @importFrom rlang .data

extract_peptides <- function(df,
                             min_peptide = 8, max_peptide = 30,
                             window_size = (2 * min_peptide) - 1,
                             save_folder = NULL){

  # ========================================================================== #
  # Sanity checks and initial definitions
  assertthat::assert_that(is.data.frame(df),
                          assertthat::is.count(min_peptide),
                          assertthat::is.count(max_peptide),
                          min_peptide <= max_peptide,
                          assertthat::is.count(window_size),
                          is.null(save_folder) | (is.character(save_folder)),
                          is.null(save_folder) | length(save_folder) == 1)

  # ========================================================================== #
  # Identify contiguous labelled peptides in each protein
  cat("\nIdentifying contiguous labeled regions...")
  df <- df %>%
    dplyr::group_by(.data$Info_protein_id) %>%
    dplyr::mutate(IsBreak = find_breaks(.data$Class),
                  Info_PepID   = ifelse(is.na(.data$Class), NA,
                                        paste0(.data$Info_protein_id, ":",
                                               cumsum(.data$IsBreak)))) %>%
    dplyr::group_by(.data$Info_PepID) %>%
    dplyr::mutate(Info_peptide_length = (!is.na(.data$Class)) * dplyr::n()) %>%
    dplyr::ungroup()


  # Extract individual contiguous peptides of length between min_peptide and
  # max_peptide
  cat("\nExtracting labelled peptides...")
  peptides <- df %>%
    dplyr::filter(!is.na(.data$Class),
                  .data$Info_peptide_length >= min_peptide,
                  .data$Info_peptide_length <= max_peptide) %>%
    dplyr::group_by(.data$Info_PepID) %>%
    dplyr::summarise(Info_organism_id    = dplyr::first(.data$Info_organism_id),
                     Info_protein_id     = dplyr::first(.data$Info_protein_id),
                     Info_start_pos      = dplyr::first(.data$Info_pos),
                     Info_end_pos        = dplyr::last(.data$Info_pos),
                     Info_peptide        = paste(.data$Info_AA, collapse = ""),
                     Info_peptide_length = dplyr::n(),
                     Class               = dplyr::first(.data$Class),
                     .groups = "drop")


  # Update df to include local neighbourhood and remove NA regions and those
  # that do not comply with min_peptide and max_peptide
  cat("\nExtracting windows...")
  df <- df %>%
    dplyr::group_by(.data$Info_protein_id) %>%
    dplyr::mutate(Info_window = make_windows(.data$Info_AA,
                                             .data$Class,
                                             window_size)) %>%
    dplyr::filter(!is.na(.data$Class),
                  .data$Info_peptide_length >= min_peptide,
                  .data$Info_peptide_length <= max_peptide) %>%
    dplyr::select(.data$Info_PepID, dplyr::everything(),
                  -.data$IsBreak, -.data$Info_peptide_length, .data$Class)


  # Check save folder and save files
  if(!is.null(save_folder)) {
    if(!dir.exists(save_folder)) dir.create(save_folder, recursive = TRUE)
    saveRDS(peptides, paste0(normalizePath(save_folder), "/peptides.rds"))
    saveRDS(df, paste0(normalizePath(save_folder), "/windowed_df.rds"))
  }

  return(list(df = df, peptides = peptides))
}


