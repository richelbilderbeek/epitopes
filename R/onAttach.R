.onAttach <- function(...) {
  op2 <- options()
  op2$pboptions_oldchar   <- op2$pboptions$char
  op2$pboptions_olduse_lb <- op2$pboptions$use_lb
  op2$pboptions$char   <- ">"
  op2$pboptions$use_lb <- TRUE

  # Prevent data table multithreads from crashing parallel routines
  op2$DTthreads_old <- data.table::getDTthreads(verbose = FALSE)
  data.table::setDTthreads(threads = 1)

  options(op2)
}


.onDetach <- function(...) {
  op2 <- options()

  data.table::setDTthreads(threads = op2$DTthreads_old)

  op2$pboptions$char      <- op2$pboptions_oldchar
  op2$pboptions$use_lb    <- op2$pboptions_olduse_lb
  op2$pboptions_oldchar   <- NULL
  op2$pboptions_olduse_lb <- NULL
  op2$DTthreads_old       <- NULL

  options(op2)
}
