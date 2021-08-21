#' Fit Random Forest model to epitope data
#'
#' Fits a Random Forest model to epitope data, using the data splits and
#' features previously calculated using [make_data_splits()] and
#' [calc_features()], respectively.
#'
#' Function [make_data_splits()] defines data splits based on (protein or
#' peptide) similarity. The split identifiers are stored in
#' `peptides.list$df$Info_split`. Function [calc_features()] calculates
#' the local and/or global features for each entry. Local features are stored in
#' `peptides.list$df` as columns starting with `feat_local_`, whereas global
#' features are stored in `peptides.list$proteins`, as columns starting with
#' `feat_global_`.
#' **NOTE**: Global features should only be used if the splitting level
#' used was "protein", otherwise they may cause contamination of performance
#' assessment due to data leakage. The splitting level of the data can be
#' checked on `peptides.list$splits.attrs$split_level`.
#'
#' @section Dealing with class imbalance:
#' Parameter `sample.rebalancing` regulates whether the resulting model attempts
#' to compensate class imbalances. If `TRUE` the Random Forest model is subject
#' to cost-sensitive training, which is done internally by setting the
#' parameter `case.weights` in the call to [ranger::ranger()] to a vector where
#' each observation of class _i_ has a weight equal to `1 / K_i`, where
#' `K_i` is the total number of cases of class `i` in the training data.
#'
#' @param peptides.list data frame containing the training data (one or more
#' numerical predictors and one **Class** attribute).
#' @param assessment.mode mode of performance assessment to use. Accepts "CV"
#' (for cross-validation using all splits in `peptides.list$df`) or "holdout"
#' (in which case the holdout split must be informed in `holdout.split`)
#' @param holdout.split name of split to be used as a holdout set.Ignored if
#' `assessment.mode = "CV"`.
#' @param threshold probability threshold for attributing a prediction as
#' *positive*.
#' @param sample.rebalancing logical: should the model try to compensate class
#' imbalances by weighted sampling of examples when training the trees in the
#' random forest? See **Dealing with class imbalance**.
#' @param use.global.features logical: should global features (potentially
#' available in `peptides.list$proteins`) be used? Should be left as the default
#' unless the user knows exactly what they're doing. See **Details**.
#' @param ncpus number of cores to use.
#' @param rnd.seed seed for random number generator
#' @param ... other options to be passed down to [ranger::ranger()].
#'
#' @return List containing the fitted model and several performance indicators.
#'
#' @author Felipe Campelo (\email{f.campelo@@aston.ac.uk})
#'
#' @export
#'
#' @importFrom dplyr %>%

fit_model <- function(peptides.list,
                      assessment.mode = c("CV", "holdout"),
                      holdout.split = NULL,
                      threshold = 0.5,
                      sample.rebalancing = TRUE,
                      use.global.features = ifelse(peptides.list$splits.attrs$split_level == "protein", TRUE, FALSE),
                      ncpus = 1,
                      rnd.seed = NULL,
                      ...){


}
