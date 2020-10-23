.onAttach <- function(...) {
  # op2 <- op.epitopes.original <- options()
  # op2$pboptions$char      <- ">"
  # op2$pboptions$txt.width <- 30
  # op2$pboptions$use_lb <- TRUE
  # options(op2)
}


.onDetach <- function(...) {
  # if(exists("op.epitopes.original")) options(op.epitopes.original)
}
