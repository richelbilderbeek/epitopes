# Sequence Entropy
feat_SeqEntr <- function(x){
  cnts <- stringr::str_count(x, LETTERS)
  cnts <- cnts[cnts > 0] / sum(cnts)
  -sum(cnts * log2(cnts))
}

# Number of Atoms
feat_NumAtom <- function(x, aa_prop){
  cnts <- stringr::str_count(x, aa_prop$One_letter_code)
  as.data.frame(t(colSums(cnts * aa_prop[, -1])))
}

# AA composition
feat_AAComp <- function(x, AAtypes){
  as.data.frame(
    lapply(AAtypes,
           function(pat, x){
             sum(stringr::str_count(x, pat)) / nchar(x)},
           x = x))
}

# AA Descriptors
feat_AADesc <- function(x){
  x <- strsplit(x, split = "")[[1]]
  f <- sapply(x, Peptides::aaDescriptors)
  f <- data.table::as.data.table(t(rowMeans(f)))
}


# Conjoint Triads
feat_CT <- function(x, aa_prop){
  idx <- stringr::str_locate(paste(aa_prop$One_letter_code,
                                   collapse = ""),
                             strsplit(x, "")[[1]])[, 1]
  x <- paste(aa_prop$CT_group[idx], collapse = "")
  ct <- stringr::str_sub(x,
                         start = 1:(nchar(x) - 2),
                         end   = 3:nchar(x))

  # Return percent occurrence of CTs
  as.data.frame(t(as.matrix(table(ct)))) / length(ct)
}

# Molecular Weight
feat_MolWei <- function(x, aa_prop){
  cnts <- stringr::str_count(x, aa_prop$One_letter_code)
  sum(cnts * aa_prop$Amino_acid_molecular_weight)
}

# N-peptides
feat_NPep <- function(x, N){
  i1 <- 1:(nchar(x) - N + 1)
  i2 <- i1 + N - 1
  y <- stringr::str_sub(x, i1, i2)
  as.data.frame(t(as.matrix(table(y)))) / length(y)
}
