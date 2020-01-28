#' Retrieve protein sequences and data from GenBank
#'
#' This function is used to retrieve data from
#' [GenBank](https://www.ncbi.nlm.nih.gov/genbank/) for given protein IDs.
#'
#' @param uids A list of potein IDs provided a character vector.
#' @param ncpus positive integer, number of cores to use (multi-core
#'        capabilities not yet available for Windows systems.)
#' @param save_folder path to folder for saving the results.
#'
#' @return A data frame containing the extracted proteins is returned invisibly.
#'
#' @author Felipe Campelo (\email{fcampelo@@ufmg.br},
#' \email{f.campelo@@aston.ac.uk})
#'
#' @export
#'

get_proteins <- function(uids,
                         ncpus = 1,
                         save_folder = NULL){

  # ========================================================================== #
  # Sanity checks and initial definitions
  assertthat::assert_that(is.character(uids),
                          length(uids) >= 1,
                          assertthat::is.count(ncpus),
                          is.null(save_folder) | (is.character(save_folder)),
                          is.null(save_folder) | length(save_folder) == 1)


  # Set up parallel processing
  if ((.Platform$OS.type == "windows") & (ncpus > 1)){
    cat("\nAttention: multicore not currently available for Windows.\n
        Forcing ncpus = 1.")
    ncpus <- 1
  } else {
    available.cores <- parallel::detectCores()
    if (ncpus >= available.cores){
      cat("\nAttention: ncpus too large, we only have ", available.cores,
          " cores.\nUsing ", available.cores - 1,
          " cores for run_experiment().")
      ncpus <- available.cores - 1
    }
  }

  # Check save folder and create file names
  if(!is.null(save_folder)) {
    if(!dir.exists(save_folder)) dir.create(save_folder)
    df_file <- paste0(normalizePath(save_folder), "/proteins.rds")
    errfile <- paste0(normalizePath(save_folder), "/proteins_retrieval_errors.rds")
  }

  # Retrieving proteins using individual requests rather than (more efficient)
  # batch requests, to catch and treat efetch() or parsing errors more easily.
  cat("\nRetrieving proteins:\n")
  if(ncpus == 1) {
    df <- pbapply::pblapply(X = uids, FUN = retrieve_single_protein)
  } else {
    df <- pbmcapply::pbmclapply(X = uids, FUN = retrieve_single_protein,
                                mc.cores = ncpus, mc.preschedule = FALSE)
  }

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
