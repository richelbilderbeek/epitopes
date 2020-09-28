#' Retrieve taxonomic classification tables from NCBI.
#'
#' This function retrieves the taxonomic information of a vector of organism
#' IDs, from the NCBI Taxonomy data base.
#'
#' @param uids vector of organism IDs to retrieve taxonomical information
#' @param save_folder path to folder for saving the results.
#'
#' @return A list containing the information for each element of UIDs.
#'
#' @author Felipe Campelo (\email{f.campelo@@aston.ac.uk})
#'
#' @export
#'
#'
#' @examples
#' uids <- c("6282", # O. volvulus,
#'           "9606") # H. sapiens
#' tax <- get_taxonomy(uids)
#'

get_taxonomy <- function(uids, save_folder = NULL){

  # ========================================================================== #
  # Sanity checks and initial definitions
  ok_classes <- c("NULL", "numeric", "integer", "character")
  assertthat::assert_that(class(uids) %in% ok_classes,
                          length(uids) >= 1,
                          is.null(save_folder) | (is.character(save_folder)),
                          is.null(save_folder) | length(save_folder) == 1)

  # Check save folder and create file names
  if(!is.null(save_folder)) {
    if(!dir.exists(save_folder)) dir.create(save_folder)
    ymd <- gsub("-", "", Sys.Date())
    df_file <- paste0(normalizePath(save_folder), "/00_taxonomy_", ymd, ".rds")
    errfile <- paste0(normalizePath(save_folder),
                      "/00_taxonomy_not_retrieved_", ymd, ".rds")
  }

  errlist <- seq_along(uids)
  reslist <- vector("list", length = length(uids))
  nerr    <- Inf

  while(length(errlist) < nerr && length(errlist) > 0){
    nerr <- length(errlist)
    cat("\n Trying to retrieve", length(errlist), "entries from NCBI (db = taxonomy)\n")
    cc <- 0
    for (idx in errlist){
      if (!(cc %% 25) && cc != 0) cat("", cc,"\n")
      cat(".")
      errk <- FALSE
      tryCatch({
        tt  <- reutils::efetch(as.numeric(uids[idx]),
                               db      = "taxonomy",
                               retmode = "xml")
        ttp <- XML::xmlTreeParse(tt$content, useInternalNodes = TRUE)
        reslist[[idx]]$Taxonomy <- data.frame(
          ScientificName = XML::xpathSApply(ttp,
                                            "//TaxaSet/Taxon/LineageEx/Taxon/ScientificName",
                                            XML::xmlValue),
          Rank = XML::xpathSApply(ttp, "//TaxaSet/Taxon/LineageEx/Taxon/Rank",
                                  XML::xmlValue),
          UID  = XML::xpathSApply(ttp, "//TaxaSet/Taxon/LineageEx/Taxon/TaxId",
                                  XML::xmlValue),
          stringsAsFactors = FALSE)
      },
      warning = function(c) {errk <<- TRUE},
      error   = function(c) {errk <<- TRUE},
      finally = NULL)

      if(!errk){
        reslist[[idx]]$UID <- uids[idx]
      }
      cc <- cc + 1

      # NCBI limits requests to three per second
      Sys.sleep(0.333)
    }
    errlist <- which(sapply(reslist, function(x) {is.null(x$UID)}))
  }

  reslist <- reslist[-errlist]
  errlist <- uids[errlist]

  # Save results to file
  if(!is.null(save_folder)){
    saveRDS(object = reslist, file = df_file)
    saveRDS(object = errlist, file = errfile)
  }

  invisible(reslist)
}
