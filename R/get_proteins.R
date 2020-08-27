#' Retrieve protein sequences and data from GenBank
#'
#' This function is used to retrieve data from
#' [Genbank's protein database](https://www.ncbi.nlm.nih.gov/genbank/) for
#' given protein IDs. If an ID is not available from Genbank the function will
#' try to retrieve it from [uniprot](https://www.uniprot.org/), based on a
#' query to the address: `https://www.uniprot.org/uniprot/<uid>.fasta`
#' (replacing <uid> by the given ID).
#'
#' @param uids A list of protein IDs provided a character vector.
#' @param save_folder path to folder for saving the results.
#'
#' @return A list containing the extracted proteins and a vector containing any
#' IDs that were not successfully retrieved.
#'
#' @author Felipe Campelo (\email{f.campelo@@aston.ac.uk})
#'
#' @export
#'

get_proteins <- function(uids, save_folder = NULL){

  # ========================================================================== #
  # Sanity checks and initial definitions
  assertthat::assert_that(is.character(uids),
                          length(uids) >= 1,
                          is.null(save_folder) | (is.character(save_folder)),
                          is.null(save_folder) | length(save_folder) == 1)


  # Check save folder and create file names
  if(!is.null(save_folder)) {
    if(!dir.exists(save_folder)) dir.create(save_folder)
    df_file <- paste0(normalizePath(save_folder), "/proteins.rds")
    errfile <- paste0(normalizePath(save_folder), "/proteins_retrieval_errors.rds")
  }

  # Retrieving proteins using individual requests rather than (more efficient)
  # batch requests, to catch and treat efetch() or parsing errors more easily.
  cat("\nRetrieving proteins:\n")
  df <- pbmcapply::pbmclapply(X        = uids,
                              FUN      = retrieve_single_protein,
                              wait     = 0.1,
                              mc.cores = 1)

  errlist <- uids[sapply(df, is.null)]
  df      <- data.frame(data.table::rbindlist(df,
                                              use.names = TRUE,
                                              fill      = TRUE))

  # Save proteins to file and clean up outfile
  if(!is.null(save_folder)){
    saveRDS(object = df,      file = df_file)
    saveRDS(object = errlist, file = errfile)
  }

  invisible(list(proteins = df,
                 errlist  = errlist))
}
