#' Extract linear B Cell epitopes from XML files retrieved from IEDB.
#'
#' This function is used to extract information for **linear B-cell epitopes**
#' from the XML files exported using the functionality provided by
#' [IEDB](https://www.iedb.org/).
#' It assumes that the user has downloaded the _Complete Database Export_ from
#' the _XML Database Export_ field in IEDB's
#' [Database Export](https://www.iedb.org/database_export_v3.php) and extracted
#' it in a given folder, which is passed as an argument to the function.
#'
#' @param data_folder path (either relative or absolute) to the directory
#'        containing the XML files
#' @param save_folder path to folder for saving the results.
#'
#' @return A data frame containing the extracted epitopes is returned invisibly.
#'
#' @author Felipe Campelo (\email{fcampelo@@ufmg.br},
#' \email{f.campelo@@aston.ac.uk})
#'
#' @export
#'

get_linear_bcell_epitopes <- function(data_folder = "./",
                                      save_folder = NULL){


  # ========================================================================== #
  # Sanity checks and initial definitions
  assertthat::assert_that(is.character(data_folder),
                          length(data_folder) == 1,
                          dir.exists(data_folder),
                          is.null(save_folder) | (is.character(save_folder)),
                          is.null(save_folder) | length(save_folder) == 1)


  # =======================================
  # Get file list and initialise variables
  start_time  <- Sys.time()
  filelist    <- dir(data_folder, pattern = ".xml",
                     full.names = TRUE)
  filelist    <- gsub("//", "/", filelist, fixed = TRUE)

  # Check save folder and create file names
  if(!is.null(save_folder)) {
    if(!dir.exists(save_folder)) dir.create(save_folder)
    df_file <- normalizePath(paste0(save_folder, "/epitopes.rds"))
    # errfile <- normalizePath(paste0(save_folder, "/xml_load_errors.rds"))
  }

  # ==================================================
  df <- pbapply::pblapply(filelist, process_xml_file)
  df <- df[sapply(df, nrow) > 0]
  df <- data.frame(data.table::rbindlist(df))

  df$start_pos <- as.numeric(df$start_pos)
  df$end_pos   <- as.numeric(df$end_pos)

  if(!is.null(save_folder)){
    saveRDS(object = df,      file = df_file)
    # saveRDS(object = errlist, file = errfile)
  }

  invisible(df)
}
