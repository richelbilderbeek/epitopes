#' Filter epitopes by taxonomy and host
#'
#' Filters an epitope dataframe (returned by [get_LBCE()] or [get_LTCE()]) by
#' a specific taxon and/or by an specific host.
#'
#' @param epitopes data frame of epitope data (returned by [get_LBCE()] or [get_LTCE()].
#' @param taxID taxon IDs that we want to filter the epitopes by
#' @param hostID host IDs that we want to filter the epitopes by
#' @param removeID organism IDs to be removed from the data frame
#'
#' @return A data frame filtered by the criteria in taxID or hostID
#'
#' @author Felipe Campelo (\email{f.campelo@@aston.ac.uk})
#'
#' @export
#'

filter_epitopes <- function(epitopes,
                            taxID    = NULL,
                            hostID   = NULL,
                            removeID = NULL) {

  # ========================================================================== #
  # Sanity checks and initial definitions
  ok_classes <- c("NULL", "numeric", "integer", "character")
  assertthat::assert_that(is.data.frame(epitopes),
                          class(taxID)    %in% ok_classes,
                          class(hostID)   %in% ok_classes,
                          class(removeID) %in% ok_classes)

  UIDs <- unique(epitopes$sourceOrg_id)
  tax  <- get_taxonomy(UIDs)

  target_list <- unlist(lapply(tax,
                               function(x){
                                 if (any(taxID %in% x$Taxonomy$UID)) return(x$UID)}))

  if(!is.null(taxID))    idx1 <- (epitopes$sourceOrg_id %in% target_list)
  if(!is.null(hostID))   idx2 <- (epitopes$host_id %in% hostID)
  if(!is.null(removeID)) idx3 <- (epitopes$sourceOrg_id %in% removeID)

  return(epitopes[which(idx1 & idx2 & !idx3) , ])

}
