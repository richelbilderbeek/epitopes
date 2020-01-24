retrieve_single_protein <- function(uid){

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

  return(prot_row)
}
