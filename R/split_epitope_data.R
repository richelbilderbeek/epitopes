#' Split epitope data based on epitope, protein or organism IDs.
#'
#' Takes a data.table of data of class *windowed_epit_dt* (returned by
#' [make_window_df()]) and split it into mutually exclusive sets of
#' observations, based on columns *Info_sourceOrg_id*, *Info_protein_id* or
#' *Info_epitope_id*.
#'
#' If the sum of **split_perc** is less than 100 an extra split is generated
#' with the remaining observations - e.g., `split_perc = c(50, 30)` results in
#' three sets with an approximately 50/30/20% split *of the total observations.*
#' If the sum is greater than 100 the splits are linearly scaled down so that
#' the sum becomes 100. Note that the split percents correspond to the number of
#' observations, not the number of unique IDs.
#'
#' This function will attempt to approximate the desired split levels, but
#' depending on the size of **wdf** set and the desired **split_level** it may
#' not be possible (e.g., if `split_level = "org"` and a single organism
#' corresponds to 90% of the data, one of the splits will necessarily correspond
#' to at least 90% of the data, regardless of the values informed in
#' `split_perc`.
#'
#' @param wdf data table of class *windowed_epit_dt* (returned by
#'        [make_window_df()])
#' @param split_level which level should be used for splitting? Use "org" for
#'        splitting by source organism ID, "prot" by protein ID or "epit" by
#'        epitope ID. When "prot" is used the routine attempts to identify
#'        different protein versions and treat them as a single unit for
#'        splitting purposes.
#' @param split_perc numeric vector of desired splitting percentages. See
#'        Details.
#' @param save_folder path to folder for saving the results.
#'
#' @return A list object containing the splitted data tables.
#'
#' @author Felipe Campelo (\email{f.campelo@@aston.ac.uk})
#'
#' @export
#'

split_epitope_data <- function(wdf,
                               split_level = "prot",
                               split_perc = c(70, 30),
                               save_folder = NULL){
  # ========================================================================== #
  # Sanity checks and initial definitions
  assertthat::assert_that("windowed_epit_dt" %in% class(wdf),
                          is.character(split_level), length(split_level) == 1,
                          split_level %in% c("org", "prot", "epit"),
                          is.numeric(split_perc),
                          all(sapply(split_perc, assertthat::is.count)),
                          is.null(save_folder) | (is.character(save_folder)),
                          is.null(save_folder) | length(save_folder) == 1)

  # Check and adjust split sizes if necessary
  ssp <- sum(split_perc)
  if (ssp < 100){
    split_perc <-c(split_perc, 100 - ssp)
  } else if (ssp > 100){
    split_perc <- 100 * split_perc / ssp
  }
  if (ssp >= 100 & length(split_perc) == 1){
    cat("\nNo splitting performed (split_perc = 100).")
    invisible(wdf)
  }

  nsplits <- length(split_perc)

  # Check save folder and create file names
  if(!is.null(save_folder)) {
    if(!dir.exists(save_folder)) dir.create(save_folder)
    ymd <- gsub("-", "", Sys.Date())
    df_files <- character(nsplits)
    for (i in 1:nsplits){
      df_files[i] <- paste0(normalizePath(save_folder), "/", ymd,
                            "_df_Split", i, ".rds")
    }
  }

  # Determine splits by splitting column
  # (Note that protein IDs ignore the trailing version number)
  id_var <- switch(split_level,
                   org  = wdf$Info_sourceOrg_id,
                   prot = gsub("\\.[1-9]+$", "", wdf$Info_protein_id),
                   epit = wdf$Info_epitope_id)

  # Get the frequencies of occurrence of each ID
  divs <- as.data.frame(sort(table(id_var), decreasing = TRUE))
  divs$Freq <- 100 * divs$Freq / length(id_var)
  divs$attr <- FALSE

  # Determine which IDs go into which splits
  split_ids <- vector("list", nsplits)
  for(i in 1:nsplits) split_ids[[i]] <- list(id_type = split_level,
                                             IDs = character(),
                                             Perc = 0)
  oc <- 1
  sc <- 1
  nc <- 0
  while (oc <= nrow(divs)){
    if(divs$Freq[oc] <= split_perc[sc] - split_ids[[sc]]$Perc){
      split_ids[[sc]]$IDs  <- c(split_ids[[sc]]$IDs, as.character(divs$id_var[oc]))
      split_ids[[sc]]$Perc <- split_ids[[sc]]$Perc + divs$Freq[oc]
      divs$attr[oc] <- TRUE
      oc <- oc + 1
      nc <- 0
    } else {
      nc <- nc + 1
    }
    if (nc > nsplits) break
    sc <- 1 + sc %% nsplits
  }

  # Allocate any remaining ones
  rem <- which(!divs$attr)
  if (length(rem) > 0){
    sidx <- which.max(split_perc - sapply(split_ids, function(x) x$Perc))
    split_ids[[sidx]]$IDs  <- c(split_ids[[sidx]]$IDs, divs$id_var[rem])
    split_ids[[sidx]]$Perc <- split_ids[[sidx]]$Perc + divs$Freq[rem]
  }

  for (i in seq_along(split_ids)){
    split_ids[[i]]$wdf <- wdf[which(id_var %in% split_ids[[i]]$IDs), ]
    if(!is.null(save_folder)) {
      saveRDS(split_ids[[i]]$wdf, file = df_files[i])
    }
  }

  invisible(split_ids)
}


