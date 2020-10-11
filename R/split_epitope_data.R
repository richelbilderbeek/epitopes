#' Split epitope data based on epitope, protein or organism IDs.
#'
#' Takes a data.table of data of class *windowed_epit_dt* (returned by
#' [make_window_df()]) and split it into mutually exclusive sets of
#' observations, based on
#'
#' @param epitopes data frame of epitope data, returned by [prepare_join_df()]
#'        or [filter_epitopes()].
#' @param proteins data frame of protein data, returned by [get_proteins()] or
#'        [make_window_df()].
#' @param set.positive how to decide whether a position should be labeled as
#'        "Positive" (+1). Use "any" to set a position as positive if
#'        it is labeled as $+1$ in at least one entry of **epitopes**; "mode" to
#'        set it by majority voting; or "all" to only label a position as
#'        Positive if it has at least one occurrence as $+1$ and none as $-1$.
#'        Unlabelled positions receive **NA** in their *Class* column.
#' @param ncpus number of cores to use
#' @param save_folder path to folder for saving the results.
#'
#' @return A data table of class *windowed_prot_dt* with columns containing the
#' number of positive and negative labels found for each position of each
#' protein, plus a *Class* column calculated according to *set.positive*.
#'
#' @author Felipe Campelo (\email{f.campelo@@aston.ac.uk})
#'
#' @export
#'
