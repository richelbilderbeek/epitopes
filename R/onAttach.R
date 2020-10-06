.onAttach <- function(...) {
  op2 <- op.epitopes.original <- options()
  op2$pboptions$char      <- ">"
  op2$pboptions$txt.width <- 30
  options(op2)
}


.onDetach <- function(...) {
  if(exists("op.epitopes.original")) options(op.epitopes.original)
}
