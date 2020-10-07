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




# prepare_protein_for_prediction <- function(proteins,
#                                            window_size,
#                                            save_folder  = NULL,
#                                            ncpus        = 1){
#
#   # ========================================================================== #
#   # Sanity checks and initial definitions
#   assertthat::assert_that(assertthat::is.count(window_size),
#                           is.data.frame(proteins),
#                           is.null(save_folder) | (is.character(save_folder)),
#                           is.null(save_folder) | length(save_folder) == 1,
#                           assertthat::is.count(ncpus))
#
#   # Set up parallel processing
#   if ((.Platform$OS.type == "windows") & (ncpus > 1)){
#     cat("\nAttention: multicore not currently available for Windows.\n
#         Forcing ncpus = 1.")
#     ncpus <- 1
#   } else {
#     available.cores <- parallel::detectCores()
#     if (ncpus >= available.cores){
#       cat("\nAttention: ncpus too large, we only have ", available.cores,
#           " cores.\nUsing ", available.cores - 1,
#           " cores for run_experiment().")
#       ncpus <- available.cores - 1
#     }
#   }
#
#   # Check save folder and create file names
#   if(!is.null(save_folder)) {
#     if(!dir.exists(save_folder)) dir.create(save_folder)
#     df_file <- paste0(normalizePath(save_folder), "/df_windowed.rds")
#   }
#
#   # ========================================================================== #
#   # Initial preprocessing
#
#   # Join epitopes and proteins dataframes, preliminary feature transformation
#   df <- proteins %>%
#     dplyr::transmute(protein_id    = as.character(.data$molecule_id),
#                      protein_seq   = as.character(.data$TSeq_sequence),
#                      protein_len   = nchar(.data$TSeq_sequence),
#                      protein_taxid = as.character(.data$TSeq_taxid),
#                      org_name      = as.character(.data$TSeq_orgname))
#
#   # ========================================================================== #
#   # Generate dataframe by sliding windows
#   # get_prot_window() is defined at the bottom of this file.
#   cat("\nExtracting windows:\n")
#   mydf <- pbmcapply::pbmclapply(X = purrr::pmap(as.list(df), list),
#                                 FUN            = get_prot_window,
#                                 ws             = window_size,
#                                 mc.cores       = ncpus,
#                                 mc.preschedule = FALSE)
#
#   cat("\nAssembling windowed dataframe...")
#   mydf <- data.frame(data.table::rbindlist(mydf))
#
#   # Save resulting dataframe and error IDs to file
#   if(!is.null(save_folder)) {
#     saveRDS(mydf, file = df_file)
#   }
#
#   invisible(mydf)
# }
#
# get_prot_window <- function(x, ws){
#   wdf <- data.frame(window_seq = rep(NA_character_, x$protein_len),
#                     taxid      = x$protein_taxid,
#                     protein_id = x$protein_id,
#                     center_pos = rep(NA_integer_, x$protein_len),
#                     stringsAsFactors = FALSE)
#   for (i in 1:x$protein_len){
#     i1   <- max(1, min(i - floor(ws / 2), x$protein_len - ws + 1))
#     i2   <- min(x$protein_len, i1 + ws - 1)
#     wdf$window_seq[i] <- substr(x$protein_seq, i1, i2)
#     wdf$center_pos[i] <- i
#   }
#   return(wdf)
# }



# process_individual_epitope_T <- function(ep, file_id){
#
#   not_valid <- FALSE
#   # If it is not a linear T-Cell epitope, then ignore
#   # -->>> ASSUMPTION: All assays of same type
#   # -->>> ASSUMPTION: All linear T-Cell epitopes appear as
#   # -->>> "FragmentOfANaturalSequenceMolecule" and contain a
#   # -->>> field named "LinearSequence"
#   not_valid <- not_valid | is.null(ep$Assays$TCell)
#   not_valid <- not_valid | is.null(ep$EpitopeStructure$FragmentOfANaturalSequenceMolecule$LinearSequence)
#
#   if(not_valid) return(NULL)
#
#   # ============= ONLY LINEAR B-CELL EPITOPES CROSS THIS LINE ============= #
#
#   # Extract relevant fields.
#   host_id        <- nullcheck(ep$Assays$TCell$Immunization$HostOrganism$OrganismId)
#   sourceOrg_id   <- nullcheck(ep$EpitopeStructure$FragmentOfANaturalSequenceMolecule$SourceOrganismId)
#   epitope_id     <- nullcheck(ep$EpitopeId)
#   molecule_id    <- nullcheck(ep$EpitopeStructure$FragmentOfANaturalSequenceMolecule$SourceMolecule$GenBankId)
#   seq            <- nullcheck(ep$EpitopeStructure$FragmentOfANaturalSequenceMolecule$LinearSequence)
#   evid_code      <- nullcheck(ep$EpitopeEvidenceCode)
#   epit_struc_def <- nullcheck(ep$EpitopeStructureDefines)
#   qual_measure   <- nullcheck(ep$Assays$TCell$AssayInformation$QualitativeMeasurement)
#
#
#   # Extract start and end positions, checking against sequence length
#   stpos  <- as.numeric(nullcheck(ep$ReferenceStartingPosition))
#   endpos <- as.numeric(nullcheck(ep$ReferenceEndingPosition))
#   sqlen  <- ifelse(is.na(seq),
#                    yes = endpos - stpos + 1,
#                    no  = nchar(seq))
#   if (any(is.na(c(stpos, endpos))) || (endpos - stpos + 1 != sqlen)){
#     stpos  <- as.numeric(nullcheck(ep$EpitopeStructure$FragmentOfANaturalSequenceMolecule$StartingPosition))
#     endpos <- as.numeric(nullcheck(ep$EpitopeStructure$FragmentOfANaturalSequenceMolecule$EndingPosition))
#   }
#
#   sqlen  <- ifelse(is.na(seq),
#                    yes = endpos - stpos + 1,
#                    no  = nchar(seq))
#   if (all(!is.na(c(stpos, endpos))) && (endpos - stpos + 1 == sqlen)){
#     start_pos <- stpos
#     end_pos   <- endpos
#   }
#
#   if(!exists("start_pos")) start_pos <- NA
#   if(!exists("end_pos"))   end_pos   <- NA
#
#   return(data.frame(file_id        = file_id,
#                     host_id        = host_id,
#                     sourceOrg_id   = sourceOrg_id,
#                     epitope_id     = epitope_id,
#                     molecule_id    = molecule_id,
#                     start_pos      = start_pos,
#                     end_pos        = end_pos,
#                     epit_len       = sqlen,
#                     seq            = seq,
#                     evid_code      = evid_code,
#                     epit_struc_def = epit_struc_def,
#                     qual_measure   = qual_measure,
#                     stringsAsFactors = FALSE))
# }



# get_LTCE <- function(data_folder = "./",
#                      ncpus = 1,
#                      save_folder = NULL){
#
  #
  #   # ========================================================================== #
  #   # Sanity checks and initial definitions
  #   assertthat::assert_that(is.character(data_folder),
  #                           length(data_folder) == 1,
  #                           dir.exists(data_folder),
  #                           assertthat::is.count(ncpus),
  #                           is.null(save_folder) | (is.character(save_folder)),
  #                           is.null(save_folder) | length(save_folder) == 1)
  #
  #
  #   # Set up parallel processing
  #   if ((.Platform$OS.type == "windows") & (ncpus > 1)){
  #     cat("\nAttention: multicore not currently available for Windows.\n
  #         Forcing ncpus = 1.")
  #     ncpus <- 1
  #   } else {
  #     available.cores <- parallel::detectCores()
  #     if (ncpus >= available.cores){
  #       cat("\nAttention: ncpus too large, we only have ", available.cores,
  #           " cores.\nUsing ", available.cores - 1,
  #           " cores for run_experiment().")
  #       ncpus <- available.cores - 1
  #     }
  #   }
  #
  #
  #   # =======================================
  #   # Get file list and initialise variables
  #   filelist    <- dir(data_folder, pattern = ".xml", full.names = TRUE)
  #   filelist    <- gsub("//", "/", filelist, fixed = TRUE)
  #
  #   # Check save folder and create file names
  #   if(!is.null(save_folder)) {
  #     if(!dir.exists(save_folder)) dir.create(save_folder)
  #     df_file <- paste0(normalizePath(save_folder), "/epitopes.rds")
  #     #errfile <- paste0(normalizePath(save_folder), "/xml_load_errors.rds")
  #   }
  #
  #   # ==================================================
  #   cat("Processing files:\n")
  #
  #   df <- pbmcapply::pbmclapply(X = filelist, FUN = process_xml_file, type = "T",
  #                               mc.cores = ncpus, mc.preschedule = FALSE)
  #
  #   #errlist <- filelist[which(sapply(df, function(x) length(x$Epitope) == 0))]
  #
  #   cat("\nProcessing resulting list:\n")
  #   df <- pbmcapply::pbmclapply(df, data.table::rbindlist, mc.cores = 1L)
  #   df <- data.frame(data.table::rbindlist(df))
  #
  #   if(!is.null(save_folder)){
  #     saveRDS(object = df,      file = df_file)
  #     #saveRDS(object = errlist, file = errfile)
  #   }
  #
  #   invisible(df)
# }

