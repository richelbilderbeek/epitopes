# fits a single random forest and assess the performance on a holdout split
fit_RF_holdout <- function(df,
                           sample.rebalancing,
                           holdout.split,
                           threshold,
                           ncpus,
                           ...){


  # Get number of trees (for inbag setup)
  ntrees <- ifelse("num.trees" %in% ...names(),
                   ...elt(which(...names() == "num.trees")),
                   formals(ranger::ranger)$num.trees)

  # Isolate the holdout split
  df.ho <- df[df$Info_split == holdout.split, ]
  df.tr <- df[df$Info_split != holdout.split, ]

  # Determine sampled cases for each tree
  # (ignoring samples from holdout.split since they have weight = 0)
  if (sample.rebalancing){
    message("Sampling class-balanced observations for each tree...")
    inbag  <- mypblapply(1:ntrees,
                         function(i, clvec){
                           nsmp <- ceiling(min(table(clvec))* 2/3)
                           smp  <- c(sample(x       = which(clvec == 1),
                                            size    = nsmp,
                                            replace = FALSE),
                                     sample(x       = which(clvec == -1),
                                            size    = nsmp,
                                            replace = FALSE))
                           return(as.numeric(seq_along(clvec) %in% smp))},
                         clvec = df.tr$Class,
                         ncpus = ncpus)

  } else {
    message("Sampling observations for each tree...")
    inbag  <- mypblapply(1:ntrees,
                         function(i, clvec){
                           nsmp <- ceiling(length(clvec)* 2/3)
                           smp  <- sample.int(length(clvec),
                                              size = nsmp,
                                              replace = FALSE)
                           return(as.numeric(seq_along(clvec) %in% smp))},
                         clvec = df.tr$Class,
                         ncpus = ncpus)
  }

  message("Fitting Random Forest...")
  myRF <- ranger::ranger(Class ~ .,
                         data =  dplyr::select(df.tr,
                                               dplyr::starts_with("feat_"),
                                               "Class"),
                         classification = TRUE,
                         probability    = TRUE,
                         inbag          = inbag,
                         num.threads    = ncpus,
                         min.node.size  = 40,
                         oob.error      = FALSE,
                         ...)

  message("Optimising threshold...")
  preds <- stats::predict(myRF,
                          data = dplyr::select(df.tr,
                                               dplyr::starts_with("feat_"),
                                               "Class"))

  # Optimise threshold
  myf <- function(thres, mypreds, truth, ncpus){
    pred.class <- smooth_predictions(x = ifelse(mypreds >= thres, 1, -1),
                                     type = "mode", window_size = 5)
    myperf <- calc_performance(truth = truth,
                               pred  = pred.class,
                               prob  = mypreds,
                               posValue = 1,
                               negValue = -1,
                               ncpus = ncpus)
    return(myperf$mcc)
  }

  threshold <- stats::optimise(myf, interval = c(0, 1),
                               mypreds = preds$predictions[, 2],
                               truth   = df.tr$Class,
                               ncpus = ncpus,
                               maximum = TRUE)$maximum

  message("Estimating performance on split: ", holdout.split)
  preds <- stats::predict(myRF,
                          data = dplyr::select(df.ho,
                                               dplyr::starts_with("feat_"),
                                               "Class"))

  pred.class <- ifelse(preds$predictions[, 2] >= threshold, 1, -1) %>%
    smooth_predictions(type = "mode", window_size = 5)

  ho.preds <- df.ho %>%
    dplyr::select(c("Info_PepID", "Info_protein_id", "Info_pos", "Info_split",
                    "Class")) %>%
    dplyr::mutate(pred.prob  = preds$predictions[, 2],
                  pred.class = pred.class)

  myperf <- calc_performance(truth = ho.preds$Class,
                             pred  = ho.preds$pred.class,
                             prob  = ho.preds$pred.prob,
                             posValue = 1,
                             negValue = -1,
                             ret.as.list = TRUE,
                             ncpus = ncpus)

  return(list(ho.preds  = ho.preds,
              threshold = threshold,
              ho.perf   = myperf))

}
