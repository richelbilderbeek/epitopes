#' Calculate epitope classification performance indices
#'
#' This function is used to calculate several distinct BINARY classification
#' performance indicators. `NA` values are not considered in the calculations.
#'
#' @param truth vector of reference values
#' @param pred  vector of predicted values
#' @param posValue value that indicates the "positive" class.
#' @param negValue value that indicates the "negative" class.
#'
#' @return list containing several performance indicators
#'
#' @author Felipe Campelo (\email{f.campelo@@aston.ac.uk})
#'
#' @export
#'

calc_performance <- function(truth, pred, posValue = 1, negValue = -1){

  assertthat::assert_that(length(posValue) == 1,
                          length(negValue) == 1,
                          length(truth) == length(pred),
                          all(stats::na.omit(truth) %in% c(posValue, negValue)),
                          all(stats::na.omit(pred)  %in% c(posValue, negValue)))

  idx <- is.na(truth) | is.na(pred)
  truth <- truth[-idx]
  pred  <- pred[-idx]

  TN <- as.numeric(sum(truth == negValue & pred == negValue))
  FN <- as.numeric(sum(truth == posValue & pred == negValue))
  FP <- as.numeric(sum(truth == negValue & pred == posValue))
  TP <- as.numeric(sum(truth == posValue & pred == posValue))

  mccNum  <- TP * TN - FP * FN
  mccDen  <- (TP + FP) * (TP + FN) * (TN + FP) * (TN + FN)

  out <- data.frame(
    TP   = TP,
    TN   = TN,
    FP   = FP,
    FN   = FN,
    sens = TP / (TP + FN),
    spec = TN / (TN + FP),
    ppv  = TP / (TP + FP),
    npv  = TN / (TN + FN),
    f1   = 2 * TP / (2 * TP + FP + FN),
    mcc  = mccNum / ifelse(mccDen == 0, 1, sqrt(mccDen)),
    accuracy =  (TP + TN) / (TP + TN + FP + FN))

  return(out)
}
