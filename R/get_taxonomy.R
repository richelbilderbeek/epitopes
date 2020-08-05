#' Retrieve taxonomic classification tables from NCBI.
#'
#' This function retrieves the taxonomic information of a vector of organism
#' IDs, from the NCBI Taxonomy data base.
#'
#' @param UIDs vector of organism IDs to retrieve taxonomical information
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

get_taxonomy <- function(UIDs) {
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

  return(out)
}
