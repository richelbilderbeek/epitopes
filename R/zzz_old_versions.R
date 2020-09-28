# retrieve_single_protein <- function(uid, wait = 0){
#
#   # Try retrieving from GenBank
#   tmpfile <- tempfile("tmp_", tmpdir = ".", fileext = ".xml")
#   errk    <- FALSE
#   tryCatch(
#     {reutils::efetch(uid     = uid,
#                      db      = "protein",
#                      rettype = "fasta",
#                      outfile = tmpfile);
#       prot_row <- XML::xmlToDataFrame(XML::xmlParse(tmpfile),
#                                       stringsAsFactors = FALSE)},
#     warning = function(c) {errk <<- TRUE},
#     error   = function(c) {errk <<- TRUE},
#     finally = NULL)
#
#   # Clean up tmpfile
#   if(file.exists(tmpfile)) file.remove(tmpfile)
#
#
#   # If error try on uniprot:
#   if (!errk) {
#     prot_row$source_DB <- "Genbank"
#   } else {
#     errk <- FALSE
#     tryCatch({
#       myurl <- paste0("https://www.uniprot.org/uniprot/", uid, ".fasta")
#       x     <- read.csv(myurl, header = TRUE)
#       x     <- paste0(x[, 1], collapse = "")
#       prot_row <- data.frame(TSeq_seqtype  = NA,
#                              TSeq_accver   = uid,
#                              TSeq_taxid    = NA,
#                              TSeq_orgname  = NA,
#                              TSeq_defline  = NA,
#                              TSeq_length   = nchar(x),
#                              TSeq_sequence = x,
#                              source_DB = "Uniprot",
#                              stringsAsFactors = FALSE)},
#       warning = function(c) {errk <<- TRUE},
#       error   = function(c) {errk <<- TRUE},
#       finally = NULL)
#   }
#
#   # If error in retrieval return NULL
#   if (errk) return(NULL)
#
#   prot_row$molecule_id <- uid
#   if (wait > 0) Sys.sleep(wait * (1 + stats::runif(1)))
#
#   return(prot_row)
# }




#
#
# get_proteins <- function(uids, save_folder = NULL){
#
#   # ========================================================================== #
#   # Sanity checks and initial definitions
#   assertthat::assert_that(is.character(uids),
#                           length(uids) >= 1,
#                           is.null(save_folder) | (is.character(save_folder)),
#                           is.null(save_folder) | length(save_folder) == 1)
#
#
#   # Check save folder and create file names
#   if(!is.null(save_folder)) {
#     if(!dir.exists(save_folder)) dir.create(save_folder)
#     df_file <- paste0(normalizePath(save_folder), "/proteins.rds")
#     errfile <- paste0(normalizePath(save_folder), "/proteins_retrieval_errors.rds")
#   }
#
#   # Retrieving proteins using individual requests rather than (more efficient)
#   # batch requests, to catch and treat efetch() or parsing errors more easily.
#   cat("\nRetrieving proteins:\n")
#   df <- pbmcapply::pbmclapply(X        = uids,
#                               FUN      = retrieve_single_protein,
#                               wait     = 0.1,
#                               mc.cores = 1)
#
#   errlist <- uids[sapply(df, is.null)]
#   df      <- data.frame(data.table::rbindlist(df,
#                                               use.names = TRUE,
#                                               fill      = TRUE))
#
#   # Save proteins to file and clean up outfile
#   if(!is.null(save_folder)){
#     saveRDS(object = df,      file = df_file)
#     saveRDS(object = errlist, file = errfile)
#   }
#
#   invisible(list(proteins = df,
#                  errlist  = errlist))
# }
