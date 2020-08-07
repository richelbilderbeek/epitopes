#' Retrieve taxonomic classification tables from NCBI.
#'
#' This function retrieves the taxonomic information of a vector of organism
#' IDs, from the NCBI Taxonomy data base.
#'
#' @param UIDs vector of organism IDs to retrieve taxonomical information
#' @param save_file filename (including path) to save the resulting list.
#'
#' @return A list containing the information for each element of UIDs.
#'
#' @author Felipe Campelo (\email{f.campelo@@aston.ac.uk})
#'
#' @export
#'
#'
#' @examples
#' UIDs <- c("6282", # O. volvulus,
#'           "9606") # H. sapiens
#' tax <- get_taxonomy(UIDs)
#'

get_taxonomy <- function(UIDs, save_file = NULL) {

  # ========================================================================== #
  # Sanity checks and initial definitions
  ok_classes <- c("NULL", "numeric", "integer", "character")
  assertthat::assert_that(class(UIDs) %in% ok_classes,
                          is.null(save_file) || is.character(save_file),
                          is.null(save_file) || length(save_file) == 1)

  cat("\n")
  out <- vector(mode = "list", length = length(UIDs))
  for (i in seq_along(UIDs)){
    cat("\rRetrieving taxonomy for uid", i, "of", length(UIDs))
    errk <- FALSE
    tryCatch({
      tt  <- reutils::efetch(as.numeric(UIDs[i]),
                             db = "taxonomy",
                             retmode = "xml")
      ttp <- XML::xmlTreeParse(tt$content,
                               useInternalNodes = TRUE)
      tmp <- data.frame(
        ScientificName = XML::xpathSApply(ttp,
                                          "//TaxaSet/Taxon/LineageEx/Taxon/ScientificName",
                                          XML::xmlValue),
        Rank = XML::xpathSApply(ttp, "//TaxaSet/Taxon/LineageEx/Taxon/Rank",
                                XML::xmlValue),
        UID  = XML::xpathSApply(ttp, "//TaxaSet/Taxon/LineageEx/Taxon/TaxId",
                                XML::xmlValue))
      },
      warning = function(c) {errk <<- TRUE},
      error   = function(c) {errk <<- TRUE},
      finally = NULL)

    if(errk) {
      tmp <- data.frame(ScientificName = NA, Rank = NA, UID = NA)
    }

    out[[i]]$UID      <- UIDs[i]
    out[[i]]$Taxonomy <- tmp

    # NCBI limits requests to three per second
    Sys.sleep(0.34)

  }

  if(!is.null(save_file)) saveRDS(out, save_file)

  return(out)
}
