#' Extract linear T Cell epitopes from XML files retrieved from IEDB.
#'
#' This function is used to extract information for **linear T-cell epitopes**
#' from the XML files exported using the functionality provided by
#' [IEDB](https://www.iedb.org/).
#' It assumes that the user has downloaded the _Complete Database Export_ from
#' the _XML Database Export_ field in IEDB's
#' [Database Export](https://www.iedb.org/database_export_v3.php) and extracted
#' it in a given folder, which is passed as an argument to the function.
#'
#' @param data_folder path (either relative or absolute) to the directory
#'        containing the XML files
#' @param ncpus positive integer, number of cores to use (multi-core
#'        capabilities not yet available for Windows systems.)
#' @param save_folder path to folder for saving the results.
#'
#' @return A data frame containing the extracted epitopes is returned invisibly.
#'
#' @author Felipe Campelo (\email{f.campelo@@aston.ac.uk})
#'
#' @export
#'
#'
#' @examples
#' my.dir   <- system.file("extdata/xml_examples", package="epitopes")
#' epitopes <- get_LTCE(my.dir, ncpus = 2)
#'

get_LTCE <- function(data_folder = "./",
                     ncpus = 1,
                     save_folder = NULL){


  # ========================================================================== #
  # Sanity checks and initial definitions
  assertthat::assert_that(is.character(data_folder),
                          length(data_folder) == 1,
                          dir.exists(data_folder),
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


  # =======================================
  # Get file list and initialise variables
  filelist    <- dir(data_folder, pattern = ".xml", full.names = TRUE)
  filelist    <- gsub("//", "/", filelist, fixed = TRUE)

  # Check save folder and create file names
  if(!is.null(save_folder)) {
    if(!dir.exists(save_folder)) dir.create(save_folder)
    df_file <- paste0(normalizePath(save_folder), "/epitopes.rds")
    #errfile <- paste0(normalizePath(save_folder), "/xml_load_errors.rds")
  }

  # ==================================================
  cat("Processing files:\n")

  df <- pbmcapply::pbmclapply(X = filelist, FUN = process_xml_file, type = "T",
                              mc.cores = ncpus, mc.preschedule = FALSE)

  #errlist <- filelist[which(sapply(df, function(x) length(x$Epitope) == 0))]

  cat("\nProcessing resulting list:\n")
  df <- pbmcapply::pbmclapply(df, data.table::rbindlist, mc.cores = 1L)
  df <- data.frame(data.table::rbindlist(df))

  if(!is.null(save_folder)){
    saveRDS(object = df,      file = df_file)
    #saveRDS(object = errlist, file = errfile)
  }

  invisible(df)
}
