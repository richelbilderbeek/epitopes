retrieve_single_protein <- function(uid, wait = 0){

  # Try retrieving from GenBank
  tmpfile <- tempfile("tmp_", tmpdir = ".", fileext = ".xml")
  errk    <- FALSE
  tryCatch(
    {reutils::efetch(uid     = uid,
                     db      = "protein",
                     rettype = "fasta",
                     outfile = tmpfile);
      prot_row <- XML::xmlToDataFrame(XML::xmlParse(tmpfile),
                                      stringsAsFactors = FALSE)},
    warning = function(c) {errk <<- TRUE},
    error   = function(c) {errk <<- TRUE},
    finally = NULL)

  # Clean up tmpfile
  if(file.exists(tmpfile)) file.remove(tmpfile)


  # If error try on uniprot:
  if (!errk) {
    prot_row$source_DB <- "Genbank"
  } else {
    errk <- FALSE
    tryCatch({
      myurl <- paste0("https://www.uniprot.org/uniprot/", uid, ".fasta")
      x     <- read.csv(myurl, header = TRUE)
      x     <- paste0(x[, 1], collapse = "")
      prot_row <- data.frame(TSeq_seqtype  = NA,
                             TSeq_accver   = uid,
                             TSeq_taxid    = NA,
                             TSeq_orgname  = NA,
                             TSeq_defline  = NA,
                             TSeq_length   = nchar(x),
                             TSeq_sequence = x,
                             source_DB = "Uniprot",
                             stringsAsFactors = FALSE)},
      warning = function(c) {errk <<- TRUE},
      error   = function(c) {errk <<- TRUE},
      finally = NULL)
  }

  # If error in retrieval return NULL
  if (errk) return(NULL)

  prot_row$molecule_id <- uid
  if (wait > 0) Sys.sleep(wait * (1 + stats::runif(1)))

  return(prot_row)
}
