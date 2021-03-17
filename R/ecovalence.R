#' Stability analysis based on Wricke's model
#' @description
#' `r badge('stable')`
#'
#' The function computes the ecovalence (Wricke, 1965) for stability analysis.
#'
#'
#' @param .data The dataset containing the columns related to Environments,
#'   Genotypes, replication/block and response variable(s).
#' @param env The name of the column that contains the levels of the
#'   environments.
#' @param gen The name of the column that contains the levels of the genotypes.
#' @param rep The name of the column that contains the levels of the
#'   replications/blocks.
#' @param resp The response variable(s). To analyze multiple variables in a
#'   single procedure use, for example, `resp = c(var1, var2, var3)`.
#' @param verbose Logical argument. If `verbose = FALSE` the code will run
#'   silently.
#' @return An object of class `ecovalence` containing the results for each
#'   variable used in the argument `resp`.
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @references Wricke, G. 1965. Zur berechnung der okovalenz bei sommerweizen
#'   und hafer. Z. Pflanzenzuchtg 52:127-138.
#' @export
#' @examples
#' \donttest{
#' library(metan)
#'out <- ecovalence(data_ge2,
#'                  env = ENV,
#'                  gen = GEN,
#'                  rep = REP,
#'                  resp = PH)
#'}
#'
ecovalence <- function(.data, env, gen, rep, resp, verbose = TRUE) {
  factors  <-
    .data %>%
    select({{env}}, {{gen}}, {{rep}}) %>%
    mutate(across(everything(), as.factor))
  vars <-
    .data %>%
    select({{resp}}, -names(factors)) %>%
    select_numeric_cols()
  factors %<>% set_names("ENV", "GEN", "REP")
  listres <- list()
  nvar <- ncol(vars)
  if (verbose == TRUE) {
    pb <- progress(max = nvar, style = 4)
  }
  for (var in 1:nvar) {
    data <- factors %>%
      mutate(Y = vars[[var]])
    if(has_na(data)){
      data <- remove_rows_na(data)
      has_text_in_num(data)
    }
    data2 <- data %>%
      means_by(ENV, GEN) %>%
      as.data.frame()
    data3 <- mutate(data2,
                    ge = residuals(lm(Y ~ ENV + GEN, data = data2)))
    ge_effect <- make_mat(data3, GEN, ENV, ge)
    Ecoval <- rowSums(ge_effect^2 * nlevels(data$REP), na.rm = TRUE)
    Ecov_perc <- (Ecoval/sum(Ecoval)) * 100
    rank <- rank(Ecoval)
    temp <- cbind(ge_effect, Ecoval, Ecov_perc, rank) %>%
      as_tibble(rownames = NA) %>%
      rownames_to_column("GEN")
    if (verbose == TRUE) {
      run_progress(pb,
                   actual = var,
                   text = paste("Evaluating trait", names(vars[var])))
    }
    listres[[paste(names(vars[var]))]] <- temp
  }
  return(structure(listres, class = "ecovalence"))
}







#' Print an object of class ecovalence
#'
#' Print the `ecovalence` object in two ways. By default, the results
#' are shown in the R console. The results can also be exported to the directory
#' into a *.txt file.
#'
#'
#' @param x The `ecovalence` x
#' @param export A logical argument. If `TRUE`, a *.txt file is exported to
#'   the working directory.
#' @param file.name The name of the file if `export = TRUE`
#' @param digits The significant digits to be shown.
#' @param ... Options used by the tibble package to format the output. See
#'   [`tibble::print()`][tibble::formatting] for more details.
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @method print ecovalence
#' @export
#' @examples
#' \donttest{
#' library(metan)
#' eco <- ecovalence(data_ge2,
#'                   env = ENV,
#'                   gen = GEN,
#'                   rep = REP,
#'                   resp = PH)
#' print(eco)
#' }
print.ecovalence <- function(x, export = FALSE, file.name = NULL, digits = 3, ...) {
  if (!class(x) == "ecovalence") {
    stop("The object must be of class 'ecovalence'")
  }
  opar <- options(pillar.sigfig = digits)
  on.exit(options(opar))
  if (export == TRUE) {
    file.name <- ifelse(is.null(file.name) == TRUE, "ecovalence print", file.name)
    sink(paste0(file.name, ".txt"))
  }
  for (i in 1:length(x)) {
    var <- x[[i]]
    cat("Variable", names(x)[i], "\n")
    cat("---------------------------------------------------------------------------\n")
    cat("Genotypic confidence index\n")
    cat("---------------------------------------------------------------------------\n")
    print(var)
  }
  cat("\n\n\n")
  if (export == TRUE) {
    sink()
  }
}
