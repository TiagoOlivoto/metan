#' Fox's stability function
#'
#' Performs a stability analysis based on the criteria of Fox et al. (1990),
#' using the statistical "TOP third" only. A stratified ranking of the genotypes
#' at each environment is done. The proportion of locations at which the
#' genotype occurred in the top third are expressed in the output.
#'
#'
#' @param .data The dataset containing the columns related to Environments,
#'   Genotypes, replication/block and response variable(s).
#' @param env The name of the column that contains the levels of the
#'   environments.
#' @param gen The name of the column that contains the levels of the genotypes.
#' @param resp The response variable(s). To analyze multiple variables in a
#'   single procedure use, for example, \code{resp = c(var1, var2, var3)}.
#' @param rep \strong{Deprecated argument. It will be retired in the next release.}
#' @param verbose Logical argument. If \code{verbose = FALSE} the code will run
#'   silently.
#' @return An object of class \code{Fox}, which is a list containing the results
#'   for each variable used in the argument \code{resp}. For each variable, a
#'   tibble with the following columns is returned.
#' * \strong{GEN} the genotype's code.
#' * \strong{mean} the mean for the response variable.
#' * \strong{TOP} The proportion of locations at which the
#' @md
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @references Fox, P.N., B. Skovmand, B.K. Thompson, H.J. Braun, and R.
#'   Cormier. 1990. Yield and adaptation of hexaploid spring triticale.
#'   Euphytica 47:57-64.
#'   \href{https://link.springer.com/article/10.1007/BF00040364}{doi:10.1007/BF00040364}.
#'
#' @export
#' @examples
#' \donttest{
#' library(metan)
#' out <- Fox(data_ge2, ENV, GEN, PH)
#' print(out)
#'}
#'
Fox <- function(.data, env, gen, resp, rep = "deprecated", verbose = TRUE) {
  if(rep != "deprecated"){
    warning("`verbose` is deprecated. It will be defunct in the new release.", call. = FALSE)
  }
  factors  <-
    .data %>%
    select({{env}}, {{gen}}) %>%
    mutate_all(as.factor)
  vars <-
    .data %>%
    select({{resp}}, -names(factors)) %>%
    select_numeric_cols()
  factors %<>% set_names("ENV", "GEN")
  listres <- list()
  nvar <- ncol(vars)
  if (verbose == TRUE) {
    pb <- progress_bar$new(
      format = "Evaluating the variable :what [:bar]:percent",
      clear = FALSE, total = nvar, width = 90)
  }
  for (var in 1:nvar) {
    data <- factors %>%
      mutate(Y = vars[[var]])
    if(has_na(data)){
      data <- remove_rows_na(data)
      has_text_in_num(data)
    }
    temp <- data %>%
      group_by(ENV) %>%
      mutate(grank = rank(-Y)) %>%
      group_by(GEN) %>%
      summarise(Y = mean(Y),
                TOP = sum(grank <= 3)) %>%
      as_tibble()
    if (verbose == TRUE) {
      pb$tick(tokens = list(what = names(vars[var])))
    }
    listres[[paste(names(vars[var]))]] <- temp
  }
  return(structure(listres, class = "Fox"))
}








#' Print an object of class Fox
#'
#' Print the \code{Fox} object in two ways. By default, the results
#' are shown in the R console. The results can also be exported to the directory
#' into a *.txt file.
#'
#'
#' @param x The \code{Fox} x
#' @param export A logical argument. If \code{TRUE}, a *.txt file is exported to
#'   the working directory.
#' @param file.name The name of the file if \code{export = TRUE}
#' @param digits The significant digits to be shown.
#' @param ... Options used by the tibble package to format the output. See
#'   \code{\link[tibble:formatting]{tibble::print()}} for more details.
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @method print Fox
#' @export
#' @examples
#' \donttest{
#' library(metan)
#' library(metan)
#' out <- Fox(data_ge2, ENV, GEN, PH)
#' print(out)
#' }
print.Fox <- function(x, export = FALSE, file.name = NULL, digits = 3, ...) {
  if (!class(x) == "Fox") {
    stop("The object must be of class 'Fox'")
  }
  opar <- options(pillar.sigfig = digits)
  on.exit(options(opar))
  if (export == TRUE) {
    file.name <- ifelse(is.null(file.name) == TRUE, "Fox print", file.name)
    sink(paste0(file.name, ".txt"))
  }
  for (i in 1:length(x)) {
    var <- x[[i]]
    cat("Variable", names(x)[i], "\n")
    cat("---------------------------------------------------------------------------\n")
    cat("Fox TOP third criteria\n")
    cat("---------------------------------------------------------------------------\n")
    print(var)
  }
  cat("\n\n\n")
  if (export == TRUE) {
    sink()
  }
}
