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
                             similarity_threshold = .7,
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
  if(!("maxit" %in% names(SAopts))) SAopts$maxit <- 2000 * round(log10(length(split_perc) ^ nrow(Y)))
  y <- optimise_splits(Y, Nstar = split_perc, alpha, SAopts)

  Y$allocation <- y$x

  X <- dplyr::left_join(dplyr::select(X, -c("SEQs")),
                        dplyr::select(Y, c("group", "allocation")),
                        by = "group")
  if(split_level == "peptide"){
    df <- dplyr::left_join(df, X, by = c("Info_PepID" = "IDs"))
  } else if(split_level == "protein"){
    df <- dplyr::left_join(df, X, by = c("Info_protein_id" = "IDs"))
  }

  # Build splits
  splits <- lapply(seq_along(split_perc), function(i){dplyr::filter(df, allocation == i)})
  names(splits) <- paste0("split_",
                          sprintf("%02d", seq_along(split_perc)), "_",
                          sprintf("%02d", round(100*split_perc)))
  names(y$xl)          <- names(splits)
  names(y$solstats$Gj) <- names(splits)
  names(y$solstats$pj) <- names(splits)

  return(list(splits = splits,
              allocation = y$xl,
              split_perc = y$solstats$Gj,
              split_balance = y$solstats$pj,
              overall_balance = y$solstats$Pstar))
}
