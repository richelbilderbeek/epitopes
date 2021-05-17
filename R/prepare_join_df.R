#' Consolidate epitope and protein data
#'
#' Label protein positions based on linear B-cell epitope information extracted
#' from IEDB.
#'
#' This routine is used to consolidate epitope information, by:
#'  \itemize{
#'    \item Filtering `proteins` such that only entries that appear at least
#'    once in `epitopes$protein_id` are retained
#'    \item Using the `epitopes` data (fields `protein_id`, `start_pos` and
#'    `end_pos`) to generate labels for each position of the selected proteins.
#'  }
#'
#' Prior to the consolidation, entries in `epitopes` input are removed if they:
#' \itemize{
#'   \item Lack a valid *protein_id* (i.e., if they have no corresponding
#'   entry in `proteins`)
#'   \item Lack a valid string in *epit_seq*
#'   \item Lack a valid definition in *epit_struc_def*
#'   \item Have sequences shorter than `min_epit` or longer than `max_epit`.
#'   \item Have a mismatch between the sequence in *epit_seq* and the
#'   corresponding sequence between *start_pos* and *end_pos* on the protein
#'   sequence.
#' }
#'
#' Entries in the `proteins` data frame with a value of *TSeq_seqtype* different
#' from "protein" are also removed.
#'
#' @param epitopes data frame of epitope data (returned by [get_LBCE()]).
#' @param proteins data frame of protein data (returned by [get_proteins()]).
#' @param save_folder path to folder for saving the results.
#' @param min_epit positive integer, shortest epitope to be considered
#' @param max_epit positive integer, longest epitope to be considered
#' @param only_exact logical, should only sequences labelled as "Exact Epitope"
#'        in variable *epit_struc_def* (within `epitopes`) be considered?
#' @param set.positive how to decide whether an observation should be of the
#'        "Positive" (+1) or "Negative" (-1) class? Use "any" to set a sequence as positive if
#'        $n_positive > 0$, "mode" to set it if $n_positive >= n_negative$,
#'        or "all" to set it if $n_negative == 0$. Defaults to "mode".
#' @param ncpus positive integer, number of cores to use
#'
#' @return A data frame containing the labeled proteins.
#'
#' @author Felipe Campelo (\email{f.campelo@@aston.ac.uk})
#'
#' @export
#'
#' @importFrom dplyr %>%
#' @importFrom rlang .data

prepare_join_df <- function(epitopes, proteins,
                            save_folder     = NULL,
                            min_epit        = 8,
                            max_epit        = 25,
                            only_exact      = FALSE,
                            set.positive    = c("any", "mode", "all"),
                            ncpus           = 1){

  # ========================================================================== #
  # Sanity checks and initial definitions
  assertthat::assert_that(is.data.frame(epitopes),
                          is.data.frame(proteins),
                          assertthat::is.count(min_epit),
                          assertthat::is.count(max_epit),
                          min_epit <= max_epit,
                          is.logical(only_exact), length(only_exact) == 1,
                          is.null(save_folder) | (is.character(save_folder)),
                          is.null(save_folder) | length(save_folder) == 1,
                          is.character(set.positive))

  if(all(tolower(set.positive) == "any")) {
    set.positive <- "any"
  } else if(all(tolower(set.positive) == "all")) {
    set.positive <- "all"
  } else set.positive <- "mode"

  # Check save folder and create file names
  if(!is.null(save_folder)) {
    if(!dir.exists(save_folder)) dir.create(save_folder, recursive = TRUE)
    df_file <- paste0(normalizePath(save_folder), "/labelled_prots.rds")
  }

  # ========================================================================== #
  # Initial preprocessing

  # Filter epitopes
  epitopes <- dplyr::filter(epitopes,
                            .data$protein_id %in% unique(proteins$UID),         # must have a corresponding protein
                            !is.na(.data$epit_seq),                             # must have an epit_seq
                            !is.na(.data$epit_struc_def)) %>%                   # must have a valid epit_struc_def
    dplyr::mutate(prot_substr = mapply(function(id, st, en){                    # filter by correct sequence in protein
      substr(proteins[which(proteins$UID == id), "TSeq_sequence"], st, en)},
      id = .data$protein_id, st = .data$start_pos, en = .data$end_pos)) %>%
    dplyr::filter(toupper(.data$epit_seq) == toupper(.data$prot_substr)) %>%
    dplyr::select(-"prot_substr")

  if (only_exact){
    epitopes <- dplyr::filter(epitopes, epit_struc_def == "Exact Epitope")
  }

  # Filter proteins
  df <- dplyr::filter(proteins,
                      .data$UID %in% unique(epitopes$protein_id),
                      .data$TSeq_seqtype == "protein",
                      !is.na(.data$TSeq_sequence))

  # ========================================================================== #
  # Build long data frame

  # Build long protein data frame
  cat("\nBuilding long data frame: proteins")
  df <- mypblapply(split(df, seq(nrow(df))),
                   FUN = function(x){
                     data.frame(Info_organism_id = x$TSeq_taxid,
                                Info_protein_id  = x$UID,
                                Info_pos         = 1:nchar(x$TSeq_sequence),
                                Info_AA          = strsplit(x$TSeq_sequence,
                                                            split = "")[[1]])},
                   ncpus = ncpus) %>%
    dplyr::bind_rows()

  # Build long epitope data frame
  cat("\nBuilding long data frame: epitopes")
  epit_summary <- mypblapply(split(epitopes, seq(nrow(epitopes))),
                             FUN = function(x){
                               data.frame(Info_protein_id   = x$protein_id,
                                          Info_pos          = x$start_pos:x$end_pos,
                                          Info_pubmed_id    = x$pubmed_id,
                                          Info_epitope_id   = x$epitope_id,
                                          Info_sourceOrg_id = x$sourceOrg_id,
                                          Info_host_id      = x$host_id,
                                          nPos              = x$n_Positive,
                                          nNeg              = x$n_Negative)},
                             ncpus = ncpus) %>%
    dplyr::bind_rows() %>%
    #
    # Consolidate information by protein-position
    dplyr::group_by(.data$Info_protein_id, .data$Info_pos) %>%
    dplyr::summarise(Nsources = length(unique(.data$Info_pubmed_id)),
                     across(-.data$Nsources, ~ paste(.x, collapse = ",")),
                     .groups = "drop") %>%
    #
    # Remove duplicated information from specific fields
    dplyr::mutate(Info_pubmed_id    = get_uniques(.data$Info_pubmed_id),
                  Info_epitope_id   = get_uniques(.data$Info_epitope_id),
                  Info_sourceOrg_id = get_uniques(.data$Info_sourceOrg_id),
                  Info_host_id      = get_uniques(.data$Info_host_id))

  # Join epitope information onto long protein data frame
  cat("\nConsolidating data")
  df <- df %>%
    dplyr::left_join(epit_summary, by = c("Info_protein_id", "Info_pos"))

    # TODO:
    # Attribute class to each position
    # nchar(.data$epit_seq) >= min_epit,                  # filter by minimum length
    # nchar(.data$epit_seq) <= max_epit) %>%              # filter by maximum length
    # attribute unique (internal) peptide ID.





    # Set class
    if (set.positive == "any"){
      df$Class <- -1 + 2 * as.numeric(df$n_Positive > 0)
    } else if (set.positive == "all") {
      df$Class <- -1 + 2 * (df$n_Negative == 0)
    } else {
      df$Class <- -1 + 2 * (df$n_Positive >= df$n_Negative)
    }

  class(df) <- c("data.table", "data.frame", "joined_epit_dt")
  attr(df, "min_epit") <- min_epit

  if(!is.null(save_folder)){
    saveRDS(object = df,    file = df_file)
    saveRDS(object = rm.df, file = errfile)
  }

  invisible(df)
}
