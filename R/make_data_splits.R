#' Split peptide data conditionally on peptide or protein similarity.
#'
#' @param df data frame of consolidated epitope information, returned by
#' [extract_peptides()] or [calculate_features()].
#' @param proteins data frame of protein data (returned by [get_proteins()]).
#' @param save_folder path to folder for saving the results.
#' @param split_perc numeric vector of desired splitting percentages. See
#'        Details.
#' @param coverage_threshold coverage threshold for grouping proteins by
#' similarity, see Details.
#' @param identity_threshold identity threshold for grouping proteins by
#' similarity, see Details.
#' @param ncpus positive integer, number of cores to use
#'
#' @return A list object containing the data splits and additional summary
#' information.
#'
#' @author Felipe Campelo (\email{f.campelo@@aston.ac.uk})
#'
#' @importFrom dplyr %>%
#' @importFrom rlang .data
#'
#' @export
#'

make_data_splits <- function(df,
                             peptides    = NULL,
                             proteins    = NULL,
                             split_level = "protein",
                             split_perc  = c(.75, .25),
                             similarity_threshold = .6,
                             save_folder = NULL,
                             ncpus       = 1,
                             substitution_matrix = "BLOSUM62",
                             SAopts = list()){

  # ========================================================================== #
  # Sanity checks and initial definitions
  # split_level <- tolower(split_level)
  # assertthat::assert_that(is.data.frame(df),
  #                         is.data.frame(proteins),
  #                         split_level %in% c("protein", "peptide"),
  #                         is.numeric(split_perc),
  #                         all(split_perc > 0),
  #                         sum(split_perc) == 100,
  #                         assertthat::is.count(coverage_threshold),
  #                         coverage_threshold >= 0, coverage_threshold <= 100,
  #                         assertthat::is.count(identity_threshold),
  #                         identity_threshold >= 0, identity_threshold <= 100,
  #                         is.character(save_folder), length(save_folder) == 1)

  diss_t <- 1 - similarity_threshold

  # Set up split names
  nsplits     <- length(split_perc)
  split_names <- paste0("split_",
                        sprintf("%02d", 1:nsplits), "_",
                        sprintf("%02d", round(split_perc)))

  message("Performing data split at ", split_level, " level")

  if(split_level == "peptide"){
    X <- peptides %>%
      dplyr::select(IDs  = Info_PepID,
                    SEQs = Info_peptide)
  } else if(split_level == "protein"){
    X <- proteins %>%
      dplyr::select(IDs  = UID,
                    SEQs = TSeq_sequence)
  }

  # Run Smith-Waterman local alignment and build similarity score matrix
  message("Calculating similarities (normalized Smith-Waterman scores)")
  scores <- mypblapply(X   = X$SEQs,
                       FUN = function(x){
                         Biostrings::pairwiseAlignment(pattern = rep(x, length(X$SEQs)),
                                                       subject = X$SEQs,
                                                       substitutionMatrix = substitution_matrix,
                                                       type = "local",
                                                       scoreOnly = TRUE)},
                       ncpus = ncpus) %>%
    do.call(what = rbind)

  # Build denominator matrix: D_{ij} = min(scores_{i,i}, scores{j,j})
  denom <- matrix(pmin(rep(diag(scores), times = nrow(scores)),
                       rep(diag(scores), each = nrow(scores))),
                  nrow  = nrow(scores), byrow = FALSE)

  # Calculate normalized dissimilarity
  diss <- 1 - scores / denom

  message("Extracting clusters (Hierarchical, single linkage)")
  clusters <- stats::hclust(d = stats::as.dist(diss), method = "single")
  X$group  <- stats::cutree(clusters, h = diss_t)

  # Check how many positive / negative examples per group
  if (split_level == "protein") {
    Y <- df %>% dplyr::left_join(X, by = c("Info_protein_id" = "IDs"))
  } else if (split_level == "peptide"){
    Y <- df %>% dplyr::left_join(X, by = c("Info_PepID" = "IDs"))
  }
  Y <- Y %>%
    dplyr::group_by(.data$group) %>%
    dplyr::summarise(nPos = sum(Class == 1),
                     nNeg = sum(Class == -1),
                     N    = dplyr::n(),
                     P    = nPos / N)

  # Define split alllocations
  y <- optimise_splits(Y, Nstar, alpha, SAopts)

  Y$allocation <- y$x


#
#
#
#
#   id.var <- gsub("\\.[1-9]+$", "", df$Info_protein_id)
#   #
#   #   # Run blast
#   #   BLAST_path <- paste0(save_folder, "/BLASTp")
#   #   proteins <- proteins %>%
#   #     dplyr::filter(.data$UID %in% unique(df$Info_protein_id)) %>%
#   #     dplyr::select(.data$UID, .data$TSeq_sequence)
#   #   blast <- run_blast(BLAST_path, proteins, ncpus)
#   #   blast <- blast %>%
#   #     dplyr::rowwise() %>%
#   #     dplyr::mutate(Pair = paste(t(apply(cbind(.data$QueryID, .data$SubjectID),
#   #                                        MARGIN = 1, FUN = sort)),
#   #                                collapse = " "))
#   #
#   #
#   #   # Turn blast results into a dissimilarity matrix
#   #   for (i in seq_along(unique(blast$QueryID))){
#   #     for (j in seq_along(unique(blast$SubjectID))){
#   #
#   #     }
#   #   }
#
#
#
#   # Filter relevant blast entries and variables
#   ID_pairs <- blast %>%
#     dplyr::filter((.data$Query_coverage >= coverage_threshold) | (.data$Perc_identity >= identity_threshold),
#                   !duplicated(t(apply(cbind(.data$QueryID, .data$SubjectID), 1, sort)))) %>%
#     dplyr::select(.data$QueryID, .data$SubjectID) %>%
#     apply(MARGIN = 1, paste, collapse = " ")
#
#   # Get the frequencies of occurrence of each ID
#   divs <- as.data.frame(sort(table(id.var), decreasing = TRUE),
#                         stringsAsFactors = FALSE)
#   divs$Freq <- 100 * divs$Freq / length(id.var)
#
#   # Initialise splits
#   split_ids <- vector("list", nsplits)
#   names(split_ids) <- split_names
#   for(i in 1:nsplits) split_ids[[i]] <- list(id_type = split_level,
#                                              IDs = character(),
#                                              IDgroups = vector("list"),
#                                              Perc = 0)
#
#
# }
#
# # ========================================================================== #
# # Determine the similarity clusters
#
# # Get the frequencies of occurrence of each ID
# divs <- as.data.frame(sort(table(id_var), decreasing = TRUE),
#                       stringsAsFactors = FALSE)
# divs$Freq <- 100 * divs$Freq / length(id_var)
#
# # Determine which IDs go into which splits
# split_ids <- vector("list", nsplits)
# names(split_ids) <- split_names
# for(i in 1:nsplits) split_ids[[i]] <- list(id_type = split_level,
#                                            IDs = character(),
#                                            IDgroups = vector("list"),
#                                            Perc = 0)
# sc <- 1   # split counter
# nc <- 0   # number of attempts to attribute
# while (nrow(divs) > 0){
#   IDgroup <- divs$id_var[1]
#   go <- TRUE
#   cc <- 1
#   while(go){
#     idx <- grep(IDgroup[cc], blast)
#     if(length(idx) > 0){
#       x <- blast[idx] # get entries that have cands[cc]
#       blast <- blast[-idx]
#       x <- gsub(paste0("\\s*", IDgroup[cc], "\\s*"), "", x)
#       IDgroup <- c(IDgroup, x)
#       cc <- cc + 1
#     } else {
#       go <- FALSE
#     }
#   }
#
#   Ptot <- sum(divs$Freq[divs$id_var %in% IDgroup])
#
#   # If the current split (sc) can accommodate all the data in IDgroup
#   if(Ptot <= split_perc[sc] - split_ids[[sc]]$Perc){
#     split_ids[[sc]]$IDs      <- c(split_ids[[sc]]$IDs, IDgroup)
#     split_ids[[sc]]$Perc     <- split_ids[[sc]]$Perc + Ptot
#     split_ids[[sc]]$IDgroups <- c(split_ids[[sc]]$IDgroups, list(IDgroup))
#     divs <- divs[-which(divs$id_var %in% IDgroup), ]
#     nc <- 0
#   } else {
#     nc <- nc + 1
#   }
#   if (nc > nsplits) break
#
#   sc <- 1 + sc %% nsplits
# }
#
# # Allocate any remaining ones
# if (length(divs) > 0){
#   sidx <- which.max(split_perc - sapply(split_ids, function(x) x$Perc))
#   split_ids[[sidx]]$IDs  <- c(split_ids[[sidx]]$IDs, divs$id_var)
#   split_ids[[sidx]]$Perc <- split_ids[[sidx]]$Perc + sum(divs$Freq)
# }
#
# for (i in seq_along(split_ids)){
#   split_ids[[i]]$wdf <- wdf[which(id_var %in% split_ids[[i]]$IDs), ]
# }
#
# return(split_ids)
}


