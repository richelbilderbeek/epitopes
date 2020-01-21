#' Retrieve protein sequences and data from GenBank
#'
#' This function is used to retrieve data from
#' [GenBank](https://www.ncbi.nlm.nih.gov/genbank/) for given protein IDs.
#'
#' @param uids A list of potein IDs provided a character vector.
#' @param save_file file name for saving the resulting data frame.
#'
#' @return A data frame containing the extracted proteins is returned invisibly.
#'
#' @author Felipe Campelo (\email{fcampelo@@ufmg.br},
#' \email{f.campelo@@aston.ac.uk})
#'
#' @export
#'

retrieve_protein_data <- function(uids, save_file){

  # Some initialisation
  start_time  <- Sys.time()
  nmols       <- length(uids)

  # Check save file extension and create error file name
  if(!identical(substr(save_file, nchar(save_file) - 3,
                       nchar(save_file)), ".rds")) {
    save_file <- paste0(save_file, ".rds")
  }
  mydir   <- normalizePath(dirname(save_file))
  errfile <- paste0(mydir, "/protein_retrieval_errors.rds")

  # Prepare dataframe for protein data
  prot_data <- data.frame(molecule_id   = character(),
                          TSeq_seqtype  = character(),
                          TSeq_accver   = character(),
                          TSeq_taxid    = character(),
                          TSeq_orgname  = character(),
                          TSeq_defline  = character(),
                          TSeq_length   = character(),
                          TSeq_sequence = character(),
                          TSeq_sid      = character(),
                          stringsAsFactors = FALSE)

  tmpfail <- prot_data
  tmpfail[1, ] <- NA

  tmpfile <- "tmp.xml"

  # Retrieving proteins using individual requests rather than (more efficient)
  # batch requests, to catch and treat efetch() or parsing errors more easily.
  errk <- FALSE
  errlist <- numeric()
  for (i in seq_along(uids)) {

    # Clear variables from previous iteration
    if(file.exists(tmpfile)) file.remove(tmpfile)
    if(exists("tmp01")) rm(tmp01)

    tryCatch(
      {reutils::efetch(uid     = uids[i],
                       db      = "protein",
                       rettype = "fasta",
                       outfile = tmpfile);
        tmp01 <- XML::xmlToDataFrame(XML::xmlParse(tmpfile),
                                     stringsAsFactors = FALSE)},
      warning = function(c) {errk <<- TRUE},
      error   = function(c) {errk <<- TRUE},
      finally = NULL)

    # If error in retrieval or parsing then catch, record and skip
    if (errk){
      errlist <- c(errlist, i)
      tmp01   <- tmpfail
      errk    <- FALSE
    }

    # Using dplyr::bind_rows() to make sure columns are in the same order
    tmp01$molecule_id <- uids[i]
    prot_data <- dplyr::bind_rows(prot_data, tmp01)

    # Save and echo elapsed time every 25 proteins processed
    if (!(i %% 25)){
      t_elaps  <- signif(Sys.time() - start_time, digits = 3)
      t_expect <- signif(t_elaps * (nmols / i - 1), digits = 3)
      t_units  <- units(t_elaps)
      cat("\n", i, "/", nmols, " proteins processed (",
          length(errlist),  " errors). ",
          "\tElapsed time: ", t_elaps, " ", t_units,
          ".\tETF: ", t_expect, " ", t_units,
          sep = "")

      saveRDS(prot_data, file = save_file)
      saveRDS(errlist,   file = errfile)
    }

    # A short pause between requests
    Sys.sleep(.1)
  }

  # Save proteins to file and clean up outfile
  saveRDS(prot_data, file = save_file)
  saveRDS(errlist,   file = errfile)
  if (file.exists(tmpfile)) file.remove(tmpfile)

  invisible(prot_data)
}
