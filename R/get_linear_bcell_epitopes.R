#' Extract linear B Cell epitopes from XML files retrieved from IEDB.
#'
#' This function is used to extract information for **linear B-cell epitopes**
#' from the XML files exported using the functionality provided by
#' [IEDB](https://www.iedb.org/).
#' It assumes that the user has downloaded the _Complete Database Export_ from
#' the _XML Database Export_ field in IEDB's
#' [Database Export](https://www.iedb.org/database_export_v3.php) and extracted
#' it in a given folder, which is passed as an argument to the function.
#'
#' @param data_folder path (either relative or absolute) to the directory
#'        containing the XML files
#' @param save_file file name for saving the resulting data frame.
#'
#' @return A data frame containing the extracted epitopes is returned invisibly.
#'
#' @author Felipe Campelo (\email{fcampelo@@ufmg.br},
#' \email{f.campelo@@aston.ac.uk})
#'
#' @export
#'

get_linear_bcell_epitopes <- function(data_folder = "./", save_file = NULL){

  # A (little) defence against accidental overwriting of existing files.
  if(file.exists(save_file)){
    msg <- paste0("File ", save_file, "already exists. Overwrite?")
    proceed <- utils::askYesNo(msg, default = FALSE)
    if(is.na(proceed) || !proceed) stop("Execution stopped by user.")
  }

  # ==================================================
  # Get file list and initialise variables
  start_time  <- Sys.time()
  filelist    <- dir(data_folder, pattern = ".xml", full.names = TRUE)
  n_files     <- length(filelist)
  file_idx    <- seq_along(filelist)
  errlist     <- numeric()

  # Check save file extension and create error file name
  if(!identical(substr(save_file, nchar(save_file) - 3,
                       nchar(save_file)), ".rds")) {
    save_file <- paste0(save_file, ".rds")
  }
  mydir   <- normalizePath(dirname(save_file))
  errfile <- paste0(mydir, "/xml_load_errors.rds")

  # Little field check function
  nullcheck <- function(x){
    ifelse(is.null(x), yes = NA, no = x)
  }

  # Initialise a (reasonably large) data frame
  nst <- 50000
  df <- data.frame(file_id        = character(nst),
                   host_id        = character(nst),
                   sourceOrg_id   = character(nst),
                   epitope_id     = character(nst),
                   molecule_id    = character(nst),
                   start_pos      = character(nst),
                   end_pos        = character(nst),
                   seq            = character(nst),
                   evid_code      = character(nst),
                   epit_struc_def = character(nst),
                   qual_measure   = character(nst),
                   stringsAsFactors = FALSE)
  df[,] <- NA_character_

  # ==================================================
  xml_data <- NA
  k <- 1
  errk <- FALSE
  for (i in file_idx){
    cat("\nProcessing file", i, "/", n_files, ": ")
    if(exists("xml_data")) rm(xml_data)

    # Load XML file as a list
    tryCatch({list_data <- XML::xmlToList(XML::xmlParse(filelist[i]))},
             error   = function(c) {errk <<- TRUE},
             warning = function(c) {errk <<- TRUE},
             finally = NULL)

    # If data did not load seamlessly, catch the error and skip
    if (errk){
      errlist <- c(errlist, i)
      errk    <- FALSE
      next
    }

    n_epitopes <- length(list_data$Reference$Epitopes)

    # If there's no epitope in this file, skip
    if (n_epitopes < 1) next

    # For every epitope listed:
    for (j in 1:n_epitopes){
      # Isolate epitope data
      tmp <- list_data$Reference$Epitopes[[j]]

      # If it is not a B-Cell epitope, skip
      # -->>> ASSUMPTION: All assays of same type                           <<<--
      # -->>> ASSUMPTION: BCell epitopes always have with field "BCell"     <<<--
      if (is.null(tmp$Assays$BCell)) {
        cat(".")
        next
      }

      # If it is not linear, skip
      # -->>> ASSUMPTION: All linear B-Cell epitopes appear as              <<<--
      # -->>> "FragmentOfANaturalSequenceMolecule" and contain a            <<<--
      # -->>> field named "LinearSequence"                                  <<<--
      if (is.null(tmp$EpitopeStructure$FragmentOfANaturalSequenceMolecule$LinearSequence)){
        cat(".")
        next
      }

      # Extract relevant fields.
      df$file_id[k]        <- nullcheck(list_data$Reference$ReferenceId)
      df$host_id[k]        <- nullcheck(tmp$Assays$BCell$Immunization$HostOrganism$OrganismId)
      df$sourceOrg_id[k]   <- nullcheck(tmp$EpitopeStructure$FragmentOfANaturalSequenceMolecule$SourceOrganismId)
      df$epitope_id[k]     <- nullcheck(tmp$EpitopeId)
      df$molecule_id[k]    <- nullcheck(tmp$EpitopeStructure$FragmentOfANaturalSequenceMolecule$SourceMolecule$GenBankId)
      df$seq[k]            <- nullcheck(tmp$EpitopeStructure$FragmentOfANaturalSequenceMolecule$LinearSequence)
      df$evid_code[k]      <- nullcheck(tmp$EpitopeEvidenceCode)
      df$epit_struc_def[k] <- nullcheck(tmp$EpitopeStructureDefines)
      df$qual_measure[k]   <- nullcheck(tmp$Assays$BCell$AssayInformation$QualitativeMeasurement)


      # Extract start and end positions, checking against sequence length
      stpos  <- as.numeric(nullcheck(tmp$ReferenceStartingPosition))
      endpos <- as.numeric(nullcheck(tmp$ReferenceEndingPosition))
      sqlen  <- ifelse(is.na(df$seq[k]),
                       yes = endpos - stpos + 1,
                       no  = nchar(df$seq[k]))
      if (any(is.na(c(stpos, endpos))) || (endpos - stpos + 1 != sqlen)){
        stpos  <- as.numeric(nullcheck(tmp$EpitopeStructure$FragmentOfANaturalSequenceMolecule$StartingPosition))
        endpos <- as.numeric(nullcheck(tmp$EpitopeStructure$FragmentOfANaturalSequenceMolecule$EndingPosition))
      }

      sqlen  <- ifelse(is.na(df$seq[k]),
                       yes = endpos - stpos + 1,
                       no  = nchar(df$seq[k]))
      if (all(!is.na(c(stpos, endpos))) && (endpos - stpos + 1 == sqlen)){
        df$start_pos[k] <- stpos
        df$end_pos[k]   <- endpos
      }

      # Echo and increment epitope counter
      cat("!")
      k <- k + 1

      # Check if we need to grow the data frame
      if (k > nrow(df)){
        tmp <- df
        tmp[,] <- NA
        df <- rbind(df, tmp)
        rm(tmp)
      }
    }

    # Save and echo elapsed time every 500 files process
    if (!(i %% 500)){
      saveRDS(object = df,      file = save_file)
      saveRDS(object = errlist, file = errfile)

      t_elaps  <- signif(Sys.time() - start_time, digits = 3)
      t_expect <- signif(as.numeric(t_elaps) * (n_files / i - 1), digits = 3)
      t_units  <- units(t_elaps)
      cat("\nEpitopes found so far: ", k, " (", length(errlist), " errors). ",
          "\tElapsed time: ", as.numeric(t_elaps), " ", t_units,
          "\t ETF: ", t_expect, " ", t_units, sep = "")
    }
  }

  # ==================================================
  # cleanup and final saving of files.
  df           <- df[1:(k - 1), ]
  df$epit_len  <- sapply(X = df$seq, FUN = nchar)
  df$start_pos <- as.numeric(df$start_pos)
  df$end_pos   <- as.numeric(df$end_pos)

  saveRDS(object = df,      file = save_file)
  saveRDS(object = errlist, file = errfile)
}
