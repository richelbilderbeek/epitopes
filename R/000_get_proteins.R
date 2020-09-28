#' Retrieve protein sequences and data from GenBank and Uniprot
#'
#' This function is used to retrieve data from
#' [Genbank's protein database](https://www.ncbi.nlm.nih.gov/genbank/) for
#' given protein IDs. If an ID is not available from Genbank the function will
#' try to retrieve it from [uniprot](https://www.uniprot.org/), based on a
#' query to the address: `https://www.uniprot.org/uniprot/<uid>.fasta`
#' (replacing <uid> by the given ID).
#'
#' Queries are processed one by one (rather than in batch) to enable treatment
#' of individual inconsistencies (e.g., wrong UIDs, queries that return a
#' different identifier, etc.). This makes this routine substantially slower,
#' but considerably more robust to errors.
#'
#' @param uids A list of protein IDs provided a character vector.
#' @param save_folder path to folder for saving the results.
#'
#' @return A list object containing a data frame with the extracted proteins
#' plus a vector of IDs that were not successfully retrieved.
#'
#' @author Felipe Campelo (\email{f.campelo@@aston.ac.uk})
#'
#' @export
#'

get_proteins <- function(uids, save_folder = NULL){

  # ========================================================================== #
  # Sanity checks and initial definitions
  assertthat::assert_that(is.character(uids),
                          length(uids) >= 1,
                          is.null(save_folder) | (is.character(save_folder)),
                          is.null(save_folder) | length(save_folder) == 1)

  # Check save folder and create file names
  if(!is.null(save_folder)) {
    if(!dir.exists(save_folder)) dir.create(save_folder)
    ymd <- gsub("-", "", Sys.Date())
    df_file <- paste0(normalizePath(save_folder), "/00_proteins_", ymd, ".rds")
    errfile <- paste0(normalizePath(save_folder),
                      "/00_prots_not_retrieved_", ymd, ".rds")
  }

  errlist <- seq_along(uids)
  reslist <- vector("list", length = length(uids))
  nerr    <- Inf

  # First try retrieving from NCBI/protein
  while(length(errlist) < nerr && length(errlist) > 0){
    nerr <- length(errlist)
    cat("\n Trying to retrieve", length(errlist), "entries from NCBI (db = protein)\n")
    cc <- 0
    for (idx in errlist){
      if (!(cc %% 25) && cc != 0) cat("", cc,"\n")
      cat(".")
      tryCatch({
        x <- reutils::efetch(uid = uids[idx],
                             db      = "protein",
                             rettype = "fasta",
                             retmode = "xml")

        reslist[idx] <- XML::xmlToList(x$get_content())
      },
      warning = function(c) {errk <<- TRUE},
      error   = function(c) {errk <<- TRUE},
      finally = NULL)

      if(!is.null(reslist[[idx]]) &
         !("ERROR" %in% names(reslist[[idx]]))){
        reslist[[idx]]$UID <- uids[idx]
        reslist[[idx]]$DB  <- "NCBI protein"
      }
      cc <- cc + 1
    }
    errlist <- which(sapply(reslist, function(x) {is.null(x$UID)}))
  }

  if (length(errlist) > 0){
    # Try retrieving remaining ids from Uniprot
    nerr <- Inf
    while(length(errlist) < nerr && length(errlist) > 0){
      nerr <- length(errlist)
      cat("\n Trying to retrieve", length(errlist), "entries from Uniprot\n")
      cc <- 0
      for (idx in errlist){
        if (!(cc %% 25) && cc != 0) cat("", cc,"\n")
        cat(".")
        errk <- FALSE
        tryCatch({
          myurl <- paste0("https://www.uniprot.org/uniprot/",
                          uids[idx], ".fasta")
          x     <- read.csv(myurl, header = FALSE)
        },
        warning = function(c) {errk <<- TRUE},
        error   = function(c) {errk <<- TRUE},
        finally = NULL)

        if(!errk){
          seq   <- paste0(x[-1, 1], collapse = "")
          #>db|UniqueIdentifier|EntryName ProteinName OS=OrganismName OX=OrganismIdentifier [GN=GeneName ]PE=ProteinExistence SV=SequenceVersion
          mystr <- strsplit(x[1,], "|", fixed = TRUE)[[1]][3]
          reslist[[idx]]$TSeq_seqtype  <- "protein"
          reslist[[idx]]$TSeq_accver   <- NA
          reslist[[idx]]$TSeq_taxid    <- strsplit(strsplit(mystr, " OX=")[[1]][2],
                                                   " ")[[1]][1]
          reslist[[idx]]$TSeq_orgname  <- strsplit(strsplit(mystr, " OS=")[[1]][2],
                                                   " OX=")[[1]][1]
          reslist[[idx]]$TSeq_defline  <- strsplit(mystr, " OS=")[[1]][1]
          reslist[[idx]]$TSeq_length   <- nchar(seq)
          reslist[[idx]]$TSeq_sequence <- seq

          reslist[[idx]]$UID <- uids[idx]

          dbstr <- strsplit(x[1,], "|", fixed = TRUE)[[1]][1]
          reslist[[idx]]$DB  <- ifelse(dbstr == ">tr",
                                       yes = "UniProtKB/TrEMBL",
                                       no = ifelse(dbstr == ">sp",
                                                   yes = "UniProtKB/Swiss-Prot",
                                                   no  = "Uniprot/other"))

        }
      }
      errlist <- which(sapply(reslist, function(x) {is.null(x$UID)}))
    }
  }

  reslist <- reslist[-errlist]
  errlist <- uids[errlist]

  df <- data.table::rbindlist(reslist,
                              use.names = TRUE,
                              fill      = TRUE)


  # Save results to file
  if(!is.null(save_folder)){
    saveRDS(object = df,      file = df_file)
    saveRDS(object = errlist, file = errfile)
  }

  invisible(list(proteins = df,
                 errlist  = errlist))
}





