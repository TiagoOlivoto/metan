group_factors = function(.data, ..., keep_factors = FALSE) {
  grouped <- group_by(.data, ...)
  names <- rlang::eval_bare(rlang::expr(paste(!!!group_keys(grouped), sep = " / ")))
  gd =   grouped %>%
    group_split(keep = FALSE) %>%
    rlang::set_names(names)
  if (keep_factors == FALSE){
      message("The factors ", paste0(collapse = " ", names(gd[[1]][ , unlist(lapply(gd[[1]], is.factor)) ])),
            " where excluded from the data")
    gd = lapply(gd, function(x){
      x[ , unlist(lapply(x, is.numeric))] %>% as.data.frame()
    })
  } else{
    gd = lapply(gd, function(x){
      x %>% as.data.frame()})
  }
  return(structure(gd, class = "group_factors"))
}

