#' Assemble a sliding window representation of epitope data
#'
#' This function is used to turn epitope and protein data into a dataframe where
#' each row corresponds to a fixed-length window that slides over the epitope
#' data (and its surrounding aminoacid residues in the protein).
#'
#' Window size, window expansion value and step size can be set either by the
#' user or automatically. In the latter case the window is set as
#' `max(3, (2 * min_epit) - 1)`, the step size as
#' `min(2, floor(min_epit / 2))`, and the window expansion value as
#' `min_epit - 1`.
#'
#' The initial position of the window is set as the first position where
#' the majority of the AAs covered by the window belong to the target region.
#' Similarly, the last position of the window is determined as the last one
#' where the majority of the AAs covered by the window still belong to the
#' target region.
#'
#' The window expansion value indicates how many AA extra positions are captured
#' from each side of the window. If set to zero the resulting dataframe will not
#' contain the `window_exp` column.
#'
#' As an example, assume an epitope "QGPGAPQGPGAP" contained within a given
#' protein, `window_size = 9`, `window_exp = 0` and `step_size = 2`.
#' Let "X" represent any AA not belonging to the target region, and | be a start
#' and end positions covered by the window. This routine would then generate the
#' following windows (each in a distinct row of the resulting dataframe):
#'
#' ...XXXXXQGPGAPQGPGAPXXXXX...\cr
#' . . . | . . . . . . . . . . . . | . . . . . . . . . . . . . . . . . . . XXXXQGPGA\cr
#' . . . . . . | . . . . . . . . . . . . | . . . . . . . . . . . . . . . . XXQGPGAPQ\cr
#' . . . . . . . . . | . . . . . . . . . . . . | . . . . . . . . . . . . . QGPGAPQGP\cr
#' . . . . . . . . . . . . | . . . . . . . . . . . . | . . . . . . . . . . PGAPQGPGA\cr
#' . . . . . . . . . . . . . . . | . . . . . . . . . . . . | . . . . . . . APQGPGAPX\cr
#' . . . . . . . . . . . . . . . . . . | . . . . . . . . . . . . | . . . . QGPGAPXXX\cr
#'
#'
#' @param epitopes data frame of epitope data (returned by [get_linear_bcell_epitopes()].
#' @param proteins data frame of epitope data (returned by [retrieve_protein_data()].
#' @param save_file file name for saving the resulting data frame.
#' @param min_epit positive integer, smallest epitope to be considered
#' @param max_epit positive integer, largest epitope to be considered
#' @param only_exact logical, should only "Exact epitopes" be considered?
#' @param window_size positive integer, size of window to use (see `Details`)
#' @param window_exp non-negative integer, how many extra AAs to capture
#'                   on either side of the window (see `Details`)
#' @param step_size positive integer, step size to use (see `Details`)
#' @param min_prot_len shortest protein length to be considered
#' @param max_prot_len longest protein length to be considered
#'
#' @return A data frame containing the extracted windows is returned invisibly.
#' Each row of the resulting data frame will have the epitope ID, protein ID,
#' windowed sequence, and the class associated with the epitope ID.
#'
#' @author Felipe Campelo (\email{fcampelo@@ufmg.br},
#' \email{f.campelo@@aston.ac.uk})
#'
#' @importFrom dplyr "%>%"
#' @importFrom rlang .data
#'
#' @export
#'
#' @examples
#' \dontrun{
#' epitopes <- data.frame(epitope_id = c("1234", "4321"),
#'               molecule_id = c("P.001", "A.002"),
#'               start_pos = c(14, 39),
#'               end_pos = c(18, 50),
#'               seq = c("LKLAT", "QGPGAPQGPGAP"),
#'               epit_struc_def = c("Exact Epitope",
#'                 "Epitope containing region/antigenic site"),
#'               qual_measure = factor(c("Negative", "Positive"),
#'                 levels = c("Negative", "Positive", "Positive-Low",
#'                   "Positive-High", "Positive-Intermediate")),
#'               epit_len = c(5, 12),
#'               host_id = c(NA,NA),
#'               sourceOrg_id = c(NA,NA),
#'               file_id = c(NA,NA),
#'               stringsAsFactors = FALSE)
#'
#' proteins <- data.frame(molecule_id = c("P.001", "A.002"),
#'               TSeq_sequence = c("TYGACPKYVKQNTLKLATGMRNVPEKQT",
#'                 "LPKEEKKDDPPKDPKKDDPPKEAQNKLNQPVVADENVDQGPGAPQGPGAPQGPGAPQGPGAPQGPGAPQGPGAPQ"),
#'               TSeq_taxid = c(NA,NA),
#'               TSeq_orgname = c(NA,NA),
#'               stringsAsFactors = FALSE)
#'
#' x <- assemble_windowed_dataframe(epitopes, proteins,
#'      save_file = "./xyz.rds")
#' }

assemble_windowed_dataframe <- function(epitopes, proteins, save_file,
                                        min_epit     = 5,
                                        max_epit     = 20,
                                        only_exact   = FALSE,
                                        window_size  = NULL,
                                        window_exp   = NULL,
                                        step_size    = NULL,
                                        min_prot_len = 1,
                                        max_prot_len = Inf){

  # ========================================================================== #
  # Sanity checks and initial definitions
  assertthat::assert_that(assertthat::is.count(min_epit),
                          assertthat::is.count(max_epit),
                          assertthat::is.count(min_prot_len),
                          assertthat::is.count(max_prot_len),
                          is.null(window_size) | assertthat::is.count(window_size),
                          is.null(step_size) | assertthat::is.count(step_size),
                          is.null(window_exp) | assertthat::is.count(window_exp + 1),
                          is.logical(only_exact), length(only_exact) == 1,
                          is.data.frame(epitopes),
                          is.data.frame(proteins),
                          is.character(save_file),
                          min_epit <= max_epit,
                          min_prot_len <= max_prot_len)

  if(is.null(window_size)) window_size <- max(3, (2 * min_epit) - 1)
  if(is.null(window_exp))  window_exp  <- min_epit - 1
  if(is.null(step_size))   step_size   <- min(2, floor(min_epit / 2))

  # Check save file extension and create error file name
  if(!identical(substr(save_file, nchar(save_file) - 3,
                       nchar(save_file)), ".rds")) {
    save_file <- paste0(save_file, ".rds")
  }
  mydir   <- normalizePath(dirname(save_file))
  errfile <- paste0(mydir, "/df_errors.rds")


  # ========================================================================== #
  # Initial preprocessing

  # Join epitopes and proteins dataframes, preliminary feature transformation
  df <- epitopes %>%
    dplyr::left_join(proteins, by = "molecule_id") %>%
    dplyr::transmute(epitope_id    = as.character(.data$epitope_id),
                     epitope_seq   = as.character(.data$seq),
                     epitope_start = as.numeric(.data$start_pos),
                     epitope_stop  = as.numeric(.data$end_pos),
                     epitope_len   = as.numeric(.data$epit_len),
                     epitope_def   = as.character(.data$epit_struc_def),
                     protein_id    = as.character(.data$molecule_id),
                     protein_seq   = as.character(.data$TSeq_sequence),
                     protein_len   = nchar(.data$protein_seq),
                     protein_taxid = as.character(.data$TSeq_taxid),
                     host_id       = as.character(.data$host_id),
                     org_id        = as.character(.data$sourceOrg_id),
                     org_name      = as.character(.data$TSeq_orgname),
                     file_id       = as.character(.data$file_id),
                     Class         = .data$qual_measure)


  # Record/remove entries without protein and/or epitope data,
  # then filter by epitope/protein sizes and recode the target class
  errlist <- df$epitope_id[which(is.na(df$protein_seq) |
                                   is.na(df$epitope_start) |
                                   is.na(df$epitope_stop))]
  df <- df %>%
    dplyr::filter(!is.na(.data$protein_seq),
                  !is.na(.data$epitope_start),
                  !is.na(.data$epitope_stop),
                  .data$epitope_len >= min_epit,
                  .data$epitope_len <= max_epit,
                  .data$protein_len >= min_prot_len,
                  .data$protein_len <= max_prot_len) %>%
    dplyr::mutate(Class = forcats::fct_recode(.data$Class,
                                              Positive = "Positive-Low",
                                              Positive = "Positive-High",
                                              Positive = "Positive-Intermediate"))

  # If needed, keep only "Exact Epitopes"
  if(only_exact) df <- df[df$epitope_struc_def == "Exact Epitope", ]

  # Record/remove entries where the epitope is not where it should be
  epit_def_not_found <- is.na(df$epitope_def)
  epit_misplaced <- mapply(function(ep, pr, i1, i2){ep != substr(pr, i1, i2)},
                           ep  = df$epitope_seq,
                           pr  = df$protein_seq,
                           i1  = df$epitope_start,
                           i2  = df$epitope_stop)

  errlist <- c(errlist, df$epitope_id[which(epit_def_not_found)])
  if (length(epit_misplaced) > 0) {
    errlist <- c(errlist, df$epitope_id[which(epit_misplaced)])
    df <- df[!epit_misplaced & !epit_def_not_found, ]
  } else {
    df <- df[!epit_def_not_found, ]
  }




  # ========================================================================== #
  # Generate dataframe by sliding windows
  # extract_windows() is an internal function defined in "extract_windows.R"
  cat("\nExtracting windows:")
  windows_df <- pbapply::pblapply(X = purrr::pmap(as.list(df), list),
                                  FUN = extract_windows,
                                  window_size = 9,
                                  step_size   = 2,
                                  window_exp  = 3)

  cat("\nAssembling dataframe:")
  windows_df <- data.frame(data.table::rbindlist(windows_df))

  # Save resulting dataframe and error IDs to file
  saveRDS(windows_df, file = save_file)
  saveRDS(errlist,    file = errfile)

  invisible(windows_df)
}
