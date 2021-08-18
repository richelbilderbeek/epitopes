call_protr <- function(SEQs, feat.list, txt.opts, dfnames, ncpus, overwrite){

  fn <- paste0("extract", names(feat.list))

  # Check if function exists
  if(!(fn %in% ls('package:protr'))){
    warning("Function protr::", fn, "() not found. Please check the feature ",
            "group names in input list '", txt.opts[1], ".features'.\nSkipping...")
    return(FALSE)
  }

  # Check if feature already exists
  feat.exists <- any(grepl(paste0("feat_", names(feat.list)), dfnames))
  if (!overwrite & feat.exists){
    message("Feature group ", txt.opts[1], ":", names(feat.list),
            " already present in peptide.list$", txt.opts[2], ".\nSkipping...")
    return(FALSE)
  }

  message("Calculating ", txt.opts[1], ":", names(feat.list))

  # Remove or replace invalid AA codes, depending on feature:
  if (grepl("Gap$", fn)) {
    SEQs <- sapply(SEQs, function(x){gsub("B|J|O|U|X|Y", "-", x)})
  } else {
    SEQs <- sapply(SEQs, function(x){gsub("B|J|O|U|X|Y", "", x)})
  }

  y <- mypblapply(SEQs,
                  function(x, fn, myargs){
                    myargs$x <- x
                    do.call(eval(parse(text = paste0("protr::", fn))),
                            args = myargs)},
                  fn = fn, myargs = feat.list[[1]],
                  ncpus = ncpus) %>%
    dplyr::bind_rows() %>%
    dplyr::rename_with(~ paste0("feat_", txt.opts[1], "_", names(feat.list), "_", .x))

  return(y)
}
