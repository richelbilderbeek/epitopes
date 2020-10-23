.onAttach <- function(...) {
  op2 <- options()
  op2$pboptions$char   <- ">"
  op2$pboptions$use_lb <- TRUE
  options(op2)
}


# .onDetach <- function(...) {
#
# }
