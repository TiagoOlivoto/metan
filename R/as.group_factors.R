as.group_factors = function(x, verbose = TRUE) {
  grouped <- group_by_if(x, is.factor)
  names <- rlang::eval_bare(rlang::expr(paste(!!!group_keys(grouped), sep = " / ")))
  gd =   grouped %>%
    group_split(keep = FALSE) %>%
    rlang::set_names(names)
  return(structure(gd, class = "group_factors"))
}
