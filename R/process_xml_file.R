process_xml_file <- function(filename){

  if(!file.exists(filename)) return(NULL)

  # Initialise error flag
  errk <- FALSE

  # Load XML file as a list
  tryCatch({list_data <- XML::xmlToList(XML::xmlParse(filename))},
           error   = function(c) {errk <<- TRUE},
           warning = function(c) {errk <<- TRUE},
           finally = NULL)

  if (errk) return(NULL)
  if (length(list_data$Reference$Epitopes) < 1) return(NULL)

  # ============= ONLY VALID LIST_DATA OBJECTS CROSS THIS LINE ============= #
  outlist <- lapply(list_data$Reference$Epitopes,
                    FUN = process_individual_epitope,
                    file_id = nullcheck(list_data$Reference$ReferenceId))

  return(outlist)
}
