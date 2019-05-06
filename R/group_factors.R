group_factors = function(.data, ..., keep_factors = FALSE, verbose = TRUE) {
  grouped <- group_by(.data, ...)
  names <- rlang::eval_bare(rlang::expr(paste(!!!group_keys(grouped), sep = " / ")))
  gd =   grouped %>%
    group_split(keep = TRUE) %>%
    rlang::set_names(names)
  if (keep_factors == FALSE){
    if (verbose == TRUE){
      if (sum(lapply(gd[[1]], is.factor)==TRUE)>0){
      message("The columns ", paste0(collapse = " ", names(gd[[1]][ , unlist(lapply(gd[[1]], is.factor)) ])),
            " where deleted. Use 'keep_factors = TRUE' to keep this columns in the grouped data. ")
      }
    }
    gd = lapply(gd, function(x){
      x[ , unlist(lapply(x, is.numeric))] %>% as.data.frame()
    })
  } else{
    gd = lapply(gd, function(x){
         as.data.frame(x)
         })
  }
  return(structure(gd, class = "group_factors"))
}
