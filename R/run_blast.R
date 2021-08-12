run_blast <- function(BLAST_path, proteins, ncpus){

  # Run BLASTp to determine protein similarities
  if(!dir.exists(BLAST_path)) dir.create(BLAST_path, recursive = TRUE)
  fn <- gsub("//", "/", paste0(BLAST_path, "/proteins.fasta"), fixed = TRUE)

  seqinr::write.fasta(as.list(proteins$TSeq_sequence), names = proteins$UID,
                      file.out = fn, as.string = TRUE, nbchar = 100000)

  message("===========================================\nBuilding BLASTp database")
  system(paste0("makeblastdb -in ", fn, " -dbtype prot"))

  message("===========================================\nRunning BLASTp\n")
  system(paste0("blastp -query ", fn, " -db ", fn, " -seg no ",
                "-outfmt '6 qseqid sseqid pident qcovhsp' > ",
                fn, "-BLAST"))

  blast <- utils::read.csv(paste0(fn, "-BLAST"), sep = "\t",
                           header = FALSE,
                           stringsAsFactors = FALSE)

  names(blast) <- c("QueryID", "SubjectID", "Perc_identity", "Query_coverage")

  blast$QueryID   <- gsub("\\.[1-9]+$", "", blast$QueryID)
  blast$SubjectID <- gsub("\\.[1-9]+$", "", blast$SubjectID)
  blast <- blast[-which(blast$QueryID == blast$SubjectID), ]

  return(blast)

}
