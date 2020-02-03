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

  # If error in retrieval return NULL
  if (errk) return(NULL)

  prot_row$molecule_id <- uid
  if (wait > 0) Sys.sleep(wait * (1 + stats::runif(1)))

  return(prot_row)
}
