#' Filter epitopes by taxonomy and host
#'
#' Filters an epitope dataframe (returned by [get_LBCE()] or [get_LTCE()]) by
#' a specific taxon and/or by an specific host.
#'
#' @param epitopes data frame of epitope data (returned by [get_LBCE()] or
#'        [get_LTCE()].
#' @param taxID taxon IDs that we want to filter the epitopes by.
#' @param hostID host IDs that we want to filter the epitopes by.
#' @param removeID organism IDs to be removed from the data frame.
#' @param tax_load_file optional taxonomy file (RDS file generated either by this
#'        function or by [get_taxonomy()]).
#' @param tax_save_file optional file name for saving the taxonomy. Ignored
#'        if `tax_file` is not `NULL`. Must be either `NULL` or a .RDS filename.
#'
#' @return Epitope dataframe filtered by the criteria in `taxID`, `hostID` and
#'         `removeID`.
#'
#' @author Felipe Campelo (\email{f.campelo@@aston.ac.uk})
#'
#' @export
#'

filter_epitopes <- function(epitopes,
                            taxID    = NULL,
                            hostID   = NULL,
                            removeID = NULL,
                            tax_load_file = NULL,
                            tax_save_file = NULL) {

  # ========================================================================== #
  # Sanity checks and initial definitions
  ok_classes <- c("NULL", "numeric", "integer", "character")
  assertthat::assert_that(is.data.frame(epitopes),
                          class(taxID)    %in% ok_classes,
                          class(hostID)   %in% ok_classes,
                          class(removeID) %in% ok_classes,
                          is.null(tax_load_file) || file.exists(tax_load_file),
                          is.null(tax_save_file) || is.character(tax_save_file),
                          is.null(tax_save_file) || length(tax_save_file) == 1)

  if (is.null(tax_load_file)){
    UIDs <- unique(epitopes$sourceOrg_id)
    tax  <- get_taxonomy(UIDs, save_file = tax_save_file)
  } else {
    tax <- readRDS(tax_load_file)
  }

  target_list <- unlist(lapply(tax,
                               function(x){
                                 if (any(taxID %in% x$Taxonomy$UID)) return(x$UID)}))

  idx1 <- idx2 <- rep(TRUE, nrow(epitopes))
  idx3 <- !idx1
  if(!is.null(taxID))    idx1 <- (epitopes$sourceOrg_id %in% target_list)
  if(!is.null(hostID))   idx2 <- (epitopes$host_id %in% hostID)
  if(!is.null(removeID)) idx3 <- (epitopes$sourceOrg_id %in% removeID)

  epitopes <- epitopes[which(idx1 & idx2 & !idx3) , ]


  if ("molecule_id" %in% names(epitopes)){
    epitopes <- epitopes[order(epitopes$molecule_id,
                               epitopes$epitope_id,
                               epitopes$start_pos), ]
  } else if ("protein_id" %in% names(epitopes)) {
    epitopes <- epitopes[order(epitopes$protein_id,
                               epitopes$epitope_id,
                               epitopes$start_pos), ]
  } else {
    epitopes <- epitopes[order(epitopes$epitope_id,
                               epitopes$start_pos), ]
  }

  return(epitopes)

}
