#' Assemble a sliding window representation of epitope data
#'
#' This function is used to turn epitope and protein data into a dataframe where
#' each row corresponds to a fixed-length window that slides over the epitope
#' data (and its surrounding aminoacid residues in the protein).
#'
#' Window size and step size can be either set by the user or determined
#' implicitly by the routine. In the latter case the window is set as
#' `max(3, (2 * min_epit) - 1)`, and the step size as
#' `min(2, floor(min_epit / 2))`.
#'
#' The initial position of the window is determined as the first position where
#' the majority of the AAs covered by the window belong to the target region.
#' Similarly, the last position of the window is determined as the last one
#' where the majority of the AAs covered by the window still belong to the
#' target region.
#'
#' As an example, assume an epitope "QGPGAPQGPGAP" contained within a given
#' protein, a window size of 9 and a step size of 2. Let "X" represent any AA
#' not belonging to the target region, and | be a start and end positions covered
#' by the window. This routine would then generate the following windows (each
#' in a distinct row of the resulting dataframe):
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
#'
#' @export
#'

assemble_windowed_dataframe <- function(epitopes, proteins, save_file,
                                        min_epit = 5,
                                        max_epit = 30,
                                        window_size = NULL,
                                        step_size   = NULL,
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
                          is.data.frame(epitopes),
                          is.data.frame(proteins),
                          is.character(save_file),
                          min_epit <= max_epit,
                          min_prot_len <= max_prot_len)

  if(is.null(window_size)) window_size <- max(3, (2 * min_epit) - 1)
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
    dplyr::transmute(epitope_id    = as.character(epitope_id),
                     epitope_seq   = as.character(seq),
                     epitope_start = as.numeric(start_pos),
                     epitope_stop  = as.numeric(end_pos),
                     epitope_len   = as.numeric(epit_len),
                     epitope_def   = as.character(epit_struc_def),
                     protein_id    = as.character(molecule_id),
                     protein_seq   = as.character(TSeq_sequence),
                     protein_len   = nchar(protein_seq),
                     protein_taxid = as.character(TSeq_taxid),
                     host_id       = as.character(host_id),
                     org_id        = as.character(sourceOrg_id),
                     org_name      = as.character(TSeq_orgname),
                     file_id       = as.character(file_id),
                     Class         = as.character(qual_measure))


  # Record/remove entries without protein and/or epitope data,
  # then filter by epitope/protein sizes and recode the target class
  errlist <- df$epitope_id[which(is.na(df$protein_seq) |
                                   is.na(df$epitope_start) |
                                   is.na(df$epitope_stop))]
  df <- df %>%
    dplyr::filter(!is.na(df$protein_seq),
                  !is.na(df$epitope_start),
                  !is.na(df$epitope_stop),
                  epitope_len >= min_epit,
                  epitope_len <= max_epit,
                  protein_len >= min_prot_len,
                  protein_len <= max_prot_len) %>%
    dplyr::mutate(Class = forcats::fct_recode(Class,
                                              Positive = "Positive-Low",
                                              Positive = "Positive-High",
                                              Positive = "Positive-Intermediate"))

  # Record/remove entries where the epitope is not where it should be
  epit_def_not_found <- is.na(df$epitope_def)
  epit_not_in_place  <- mapply(FUN = function(ep, pr, i1, i2){ep != substr(pr, i1, i2)},
                               ep  = df$epitope_seq,
                               pr  = df$protein_seq,
                               i1  = df$epitope_start,
                               i2  = df$epitope_stop)
  errlist <- c(errlist,
               df$epitope_id[which(epit_not_in_place)],
               df$epitope_id[which(epit_def_not_found)])

  df <- df[!epit_not_in_place & !epit_def_not_found, ]


  # ========================================================================== #
  # Generate dataframe by sliding windows

  # Function to extract windows from a given row of df
  extract_windows <- function(x, window_size, step_size){
    # Initialise dataframe
    wdf <- data.frame(window_seq = rep(NA_character_, x$epitope_len + 2),
                      Class      = x$Class,
                      epitope_id = x$epitope_id,
                      protein_id = x$protein_id,
                      stringsAsFactors = FALSE)

    # Initial position for window
    i1   <- max(1, x$epitope_start - floor(window_size / 2))
    i2   <- min(x$protein_len, i1 + window_size - 1)
    stop <- FALSE
    j <- 1
    while(!stop){
      wdf$window_seq[j] <- substr(x$protein_seq, i1, i2)
      i1 <- i1 + step_size
      i2 <- i2 + step_size
      j  <- j + 1
      if(i2 > min(x$protein_len,  x$epitope_stop + ceiling(window_size / 2))) {
        stop <- TRUE
      }
    }
    return(wdf[!is.na(wdf$window_seq), ])
  }

  windows_df <- do.call(rbind,
                        lapply(seq_along(df$epitope_id),
                               function(i, df, w, s){
                                 extract_windows(x = df[i, ], w, s)
                               },
                               df = df,
                               w = window_size,
                               s = step_size))

  # Add file save lines here


  invisible(windows_df)

}
