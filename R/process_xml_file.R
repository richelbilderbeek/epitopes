process_xml_file <- function(filename, type = "B"){

  if(!file.exists(filename)) return(NULL)

  # Initialise error flag
  errk <- FALSE

  # Load XML file as a list
  tryCatch({
    invisible(utils::capture.output(
      list_data <- XML::xmlToList(XML::xmlParse(filename))))},
    error   = function(c) {errk <<- TRUE},
    warning = function(c) {errk <<- TRUE},
    finally = NULL)

  if (errk) return(NULL)
  if (length(list_data$Reference$Epitopes) < 1) return(NULL)

  # ============= ONLY VALID LIST_DATA OBJECTS CROSS THIS LINE ============= #
  if (type == "B")
    outlist <- lapply(list_data$Reference$Epitopes,
                      FUN = process_individual_epitope_B,
                      file_id = nullcheck(list_data$Reference$ReferenceId))
  else if (type == "T"){
    outlist <- lapply(list_data$Reference$Epitopes,
                      FUN = process_individual_epitope_T,
                      file_id = nullcheck(list_data$Reference$ReferenceId))
  } else stop("epitope type not recognised")

  return(outlist)
}
