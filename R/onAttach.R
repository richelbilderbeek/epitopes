.onAttach <- function(...) {
  # Check if BLAST is installed
  errk <- FALSE
  tryCatch({
    invisible(utils::capture.output(
      blast_version <- system("blastp -version", intern = TRUE)[1]))},
    error   = function(c) {errk <<- TRUE},
    warning = function(c) {errk <<- TRUE},
    finally = NULL)

  if (errk) {
    blast_version <- paste0("\nWARNING: blastp not found",
                            "\nPlease follow the instructions in",
                            "\nhttps://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download",
                            "\nto set up BLAST+ on your machine.",
                            "\n====================================================")
  } else {
    blast_version <- paste0(blast_version, " detected",
                            "\n====================================================")
  }

  msg <- paste0("\n================= epitopes package =================",
             "\nNOTE: Data splitting functions require BLAST+.",
             "\nChecking if BLAST+ is available:\n",
             blast_version)

  packageStartupMessage(msg)
}
