run_blast <- function(BLAST_path, prots, ncpus){

  # Run BLASTp to determine protein similarities
  if(!dir.exists(BLAST_path)) dir.create(BLAST_path, recursive = TRUE)
  fn <- gsub("//", "/", paste0(BLAST_path, "/prots.fasta"), fixed = TRUE)

  seqinr::write.fasta(as.list(prots$TSeq_sequence), names = prots$UID,
                      file.out = fn, as.string = TRUE, nbchar = 100000)

  message("===========================================\nBuilding BLASTp database")
  system(paste0("makeblastdb -in ", fn, " -dbtype prot"))

  message("===========================================\nRunning BLASTp\n")
  system(paste0("blastp -num_threads ", ncpus, " -query ", fn, " -db ", fn, " -seg no ",
                "-outfmt '6 qseqid sseqid length qlen slen nident pident",
                " qcovhsp mismatch gaps qstart qend sstart send evalue score' > ",
                fn, "-BLAST"))

  blast <- utils::read.csv(paste0(fn, "-BLAST"), sep = "\t",
                           header = FALSE,
                           stringsAsFactors = FALSE)

  names(blast) <- c("QueryID", "SubjectID", "Alignment_length",
                    "Query_length", "Subject_length", "Num_matches",
                    "Perc_identity", "Query_coverage", "Num_mismatches",
                    "Num_gaps", "Query_match_start", "Query_match_end",
                    "Subject_match_start", "Subject_match_end",
                    "E_value", "Score")

  blast$QueryID   <- gsub("\\.[1-9]+$", "", blast$QueryID)
  blast$SubjectID <- gsub("\\.[1-9]+$", "", blast$SubjectID)

  return(blast)

}
