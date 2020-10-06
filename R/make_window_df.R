#' Assemble a sliding window representation of epitope or protein data
#'
#' Casts epitope or protein data into a dataframe where rows correspond to a
#' fixed-length window centred on consecutive positions of each sequence.
#'
#' The sliding window runs from the first to the last positions of each sequence
#' in the `df` input (which in turn comes from running `prepare_join_df()`).
#' If a window extends beyond the limits of the protein it is padded with the
#' first (or last) letter of the sequence.
#'
#' @param df data frame of epitope data (returned by `prepare_join_df()`) or
#'        protein data (returned by `get_proteins()`).
#' @param save_folder path to folder for saving the results.
#' @param window_size positive integer, size of window to use. If `df` is a
#'        dataframe of class "joined_epitope_dt" (returned from
#'        `prepare_join_df()`) then a standard value for this parameter is
#'        used automatically if `window_size = NULL` (this standard value is
#'        calculated as (2 x `min_epitope` - 1), based on the value of
#'        `min_epitope` used in the call to `prepare_join_df()`.
#' @param ncpus positive integer, number of cores to use
#'
#' @return A data.table object containing the sliding window representation of
#' the input data.
#'
#' @author Felipe Campelo (\email{f.campelo@@aston.ac.uk})
#'
#' @export
#'
#' @importFrom rlang .data
#' @importFrom dplyr %>%

make_window_df <- function(df,
                           save_folder     = NULL,
                           window_size     = NULL,
                           ncpus           = 1){

  # ========================================================================== #
  # Sanity checks and initial definitions
  assertthat::assert_that(is.null(window_size) | assertthat::is.count(window_size),
                          is.data.frame(df),
                          is.null(save_folder) | (is.character(save_folder)),
                          is.null(save_folder) | length(save_folder) == 1,
                          assertthat::is.count(ncpus))


  if(is.null(window_size)){
    if(!("min_epit" %in% names(attributes(df)))){
      stop("df must have a min_epit attribute when window_size == NULL")
    } else {
      assertthat::assert_that(assertthat::is.count(attr(df, "min_epit")))
      window_size <- 2 * attr(df, "min_epit") - 1
    }
  }

  type <- "prot"
  if("joined_epitope_dt" %in% class(df)){
    type <- "epit"
  } else if(!("protein_dt" %in% class(df))){
    cat("\nNote: Cannot determine if 'df' is an epitope or protein dataframe.",
        "\nTrying to treat as 'protein' (default). It may result in errors.")
  }

  # Set up parallel processing
  available.cores <- parallel::detectCores()
  if (ncpus > available.cores){
    cat("\nAttention: cores too large, we only have ", available.cores,
        " cores.\nUsing ", available.cores - 1,
        " cores for get_LBCE().")
    ncpus <- available.cores - 1
  }

  if (.Platform$OS.type == "windows"){
    cl <- parallel::makeCluster(ncpus, setup_timeout = 1)
  } else {
    cl <- ncpus
  }

  # Check save folder and create file names
  if(!is.null(save_folder)) {
    if(!dir.exists(save_folder)) dir.create(save_folder)
    ymd <- gsub("-", "", Sys.Date())
    df_file <- paste0(normalizePath(save_folder), "/windowed_", type, "_df_",
                      ymd, ".rds")
  }

  # ========================================================================== #
  # Generate dataframe by sliding windows
  # extract_windows() is an internal function defined in "extract_windows.R"
  # cat("\nExtracting windows:\n")
  # windows_df <- pbmcapply::pbmclapply(X = purrr::pmap(as.list(df), list),
  #                                     FUN = extract_windows,
  #                                     window_size = window_size,
  #                                     step_size   = step_size,
  #                                     window_exp  = window_exp,
  #                                     mc.cores       = ncpus,
  #                                     mc.preschedule = FALSE)
  #
  # cat("\nAssembling windowed dataframe...")
  # windows_df <- data.frame(data.table::rbindlist(windows_df))
  #
  # # Save resulting dataframe and error IDs to file
  # if(!is.null(save_folder)) {
  #   saveRDS(windows_df, file = df_file)
  #   saveRDS(errlist,    file = errfile)
  # }
  #
  # invisible(list(windows_df = windows_df,
  #                errlist    = errlist))
}
