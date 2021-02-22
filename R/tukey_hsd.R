#' Tukey Honest Significant Differences
#' @description
#' `r badge('experimental')`
#'
#' Helper function to perform Tukey post-hoc tests. It is used in [gafem].
#'
#' @param model an object of class `aov` or `lm`.
#' @param ... other arguments passed to the function
#'   [stats::TukeyHSD()]. These include:
#' * **which**: A character vector listing terms in the fitted model for which the
#' intervals should be calculated. Defaults to all the terms.
#' * **ordered**: A logical value indicating if the levels of the factor should be
#' ordered according to increasing average in the sample before taking
#' differences. If ordered is true then the calculated differences in the means
#' will all be positive. The significant differences will be those for which the
#' lwr end point is positive.
#' @param out The format of outputs. If `out = "long"` a 'long' format
#'   (tibble) is returned. If `out = "wide"`, a matrix with the adjusted
#'   p-values for each term is returned.
#' @md
#' @return A tibble data frame containing the results of the pairwise
#'   comparisons (if `out = "long"`) or a "list-columns" with p-values for
#'   each term (if `out = "wide"`).
#' @export
#' @examples
#'\donttest{
#' library(metan)
#' mod <- lm(PH ~ GEN + REP, data = data_g)
#' tukey_hsd(mod)
#' tukey_hsd(mod, out = "wide")
#' }
tukey_hsd <- function(model, ..., out = "long"){
  d <- match.call()
  if(!has_class(model, c("lm", "aov"))){
    stop("object '", d[["model"]], "' must be of class 'lm' or 'aov'.", call. = FALSE)
  }
  mod <-
    aov(model) %>%
    TukeyHSD(...)
  results <-  do.call(rbind,
                      lapply(seq_along(mod), function(i){
                        mod[i] %>%
                          as.data.frame() %>%
                          rownames_to_column("comparison") %>%
                          separate(comparison, into= c("group2", "group1"), sep = "-") %>%
                          add_cols(term = names(mod)[i], .before = group2) %>%
                          set_names("term", "group2", "group1", "estimate", "conf.low", "conf.high", "p.adj") %>%
                          add_cols(sign = stars_pval(p.adj)) %>%
                          reorder_cols(group1, .before = group2)
                      })
  )
  if(out == "wide"){
    results %<>%
      group_by(term) %>%
      doo(~make_mat(., group1, group2, p.adj))
    return(results)
  } else{
    return(results)
  }
}
