#' Label protein positions based on epitope data
#'
#' Takes a data.table of data of class *protein_dt* (returned by
#' [get_proteins()]) or *windowed_prot_dt* (returned by [make_window_df()]) and
#' labels each protein position based on known epitope/non-epitope regions
#' documented in an **epitopes** data.table of class *joined_epit_dt* (returned
#' by [prepare_join_df()] or by [filter_epitopes()]).
#'
#' @param epitopes data frame of epitope data, returned by [prepare_join_df()]
#'        or [filter_epitopes()].
#' @param proteins data frame of protein data, returned by [get_proteins()] or
#'        [make_window_df()].
#' @param set.positive how to decide whether a position should be labeled as
#'        "Positive" (+1). Use "any" to set a position as positive if
#'        it is labeled as $+1$ in at least one entry of **epitopes**; "mode" to
#'        set it by majority voting; or "all" to only label a position as
#'        Positive if it has at least one occurrence as $+1$ and none as $-1$.
#'        Unlabelled positions receive **NA** in their *Class* column.
#' @param ncpus number of cores to use
#' @param save_folder path to folder for saving the results.
#'
#' @return A data table of class *windowed_prot_dt* with columns containing the
#' number of positive and negative labels found for each position of each
#' protein, plus a *Class* column calculated according to *set.positive*.
#'
#' @author Felipe Campelo (\email{f.campelo@@aston.ac.uk})
#'
#' @export
#'

label_proteins <- function(proteins, epitopes,
                           set.positive = c("any", "mode", "all"),
                           ncpus = 1,
                           save_folder = NULL){

  # ========================================================================== #
  # Sanity checks and initial definitions
  prot_classes <- c("windowed_prot_dt", "protein_dt")
  epit_classes <- c("joined_epit_dt")
  assertthat::assert_that(is.data.frame(epitopes),
                          is.data.frame(proteins),
                          any(class(proteins) %in% prot_classes),
                          any(class(epitopes) %in% epit_classes),
                          is.null(save_folder) | (is.character(save_folder)),
                          is.null(save_folder) | length(save_folder) == 1,
                          is.character(set.positive),
                          assertthat::is.count(ncpus))

  ncpus <- max(1, min(ncpus, parallel::detectCores() - 1))

  if(all(tolower(set.positive) == "any")) {
    set.positive <- "any"
  } else if(all(tolower(set.positive) == "all")) {
    set.positive <- "all"
  } else set.positive <- "mode"

  # Check save folder and create file names
  if(!is.null(save_folder)) {
    if(!dir.exists(save_folder)) dir.create(save_folder)
    ymd <- gsub("-", "", Sys.Date())
    df_file <- paste0(normalizePath(save_folder), "/", ymd,
                      "df_labelled_prots.rds")
  }

  # Prepare protein data.table for labeling
  if(!("windowed_prot_dt" %in% class(proteins))){
    # Convert proteins into windowed data format
    proteins <- make_window_df(proteins, window_size = 1, ncpus = ncpus)
  }

  names(proteins)[which(names(proteins) == "Info_window_seq")] <- "Info_AA"

  cat("\nExtracting label positions from 'epitopes'\n")
  myf <- function(x){
    st <- as.numeric(x$start_pos)
    en <- as.numeric(x$end_pos)
    data.frame(UID     = x$protein_id,
               pos     = st:en,
               epit_id = x$epitope_id,
               nPos    = x$n_Positive,
               nNeg    = x$n_Negative)
  }
  if (ncpus > 1) {
    cl <- set_mc(ncpus)
    tmp <- pbapply::pblapply(cl  = cl,
                             X   = purrr::pmap(as.list(epitopes), list),
                             FUN = myf)
    close_mc(cl)
  } else {
    tmp <- pbapply::pblapply(cl  = 1,
                             X   = purrr::pmap(as.list(epitopes), list),
                             FUN = myf)
  }
  cat("Done!\n")

  # Aggregate multiply-labelled entries. The variable names are initialised
  # below just to prevent NOTEs on CRAN. The dplyr functions use references to
  # variables internal to df)
  UID <- pos <- nPos <- nNeg <- epit_id <- NULL
  tmp <- data.table::rbindlist(tmp, use.names = TRUE)
  tmp <- tmp[, list(Info_epit_id = paste(unique(epit_id), collapse = ","),
                    Info_nPos    = sum(nPos),
                    Info_nNeg    = sum(nNeg)),
             by = list(UID, pos)]

  # Join epitope labels into proteins dataset
  proteins <- dplyr::left_join(proteins, tmp,
                               by = c("Info_UID" = "UID",
                                      "Info_center_pos" = "pos"))

  # Set class
  if (set.positive == "any"){
    proteins$Class <- -1 + 2 * as.numeric(proteins$Info_nPos > 0)
  } else if (set.positive == "all") {
    proteins$Class <- -1 + 2 * (proteins$Info_nNeg == 0)
  } else {
    proteins$Class <- -1 + 2 * (proteins$Info_nPos >= proteins$Info_nNeg)
  }

  return(proteins)

}
