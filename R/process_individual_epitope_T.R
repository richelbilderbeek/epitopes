process_individual_epitope_T <- function(ep, file_id){
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
}
