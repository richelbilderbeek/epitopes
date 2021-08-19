#' @importFrom dplyr %>%
#' @importFrom rlang .data
extractAtoms <- function(x){
  y <- data.frame(
    Code1    = c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L",
                 "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"),
    nC = c(3, 3, 4, 5, 9, 2, 6, 6, 6, 6, 5, 4, 5, 5, 6, 3, 4, 5, 11, 9),
    nH = c(7, 7, 7, 9, 11, 5, 9, 13, 14, 13, 11, 8, 9, 10, 14, 7, 9, 11, 12, 11),
    nN = c(1, 1, 1, 1, 1, 1, 3, 1, 2, 1, 1, 2, 1, 2, 4, 1, 1, 1, 2, 1),
    nO = c(2, 2, 4, 4, 2, 2, 2, 2, 2, 2, 2, 3, 2, 3, 2, 3, 3, 2, 2, 3),
    nS = c(0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0))

  cnts <- as.data.frame(as.matrix(table(strsplit(x, split = "")[[1]])))
  cnts$Code1 <- rownames(cnts)
  cnts <- y %>%
    dplyr::left_join(cnts, by = "Code1") %>%
    dplyr::mutate(dplyr::across(starts_with("n"), ~.x*.data$V1)) %>%
    dplyr::select(-c("Code1", "V1")) %>%
    dplyr::summarise(dplyr::across(dplyr::everything(), sum, na.rm = TRUE))

  return(cnts)
}
