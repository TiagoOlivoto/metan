#' Shukla's stability variance parameter
#'
#' The function computes the Shukla's stability variance parameter (1972) and
#' uses the Kang's nonparametric stability (rank sum) to imcorporate the mean
#' performance and stability into a single selection criteria.
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
#'   single procedure use, for example, \code{resp = c(var1, var2, var3)}.
#' @param verbose Logical argument. If \code{verbose = FALSE} the code will run
#'   silently.
#' @return An object of class \code{Shukla}, which is a list containing the results for each
#'   variable used in the argument \code{resp}. For each variable, a tibble with the following
#'   columns is returned.
#' * \strong{GEN} the genotype's code.
#' * \strong{Y} the mean for the response variable.
#' * \strong{ShuklaVar} The Shukla's stability variance parameter.
#' * \strong{rMean} The rank for \strong{Y} (decreasing).
#' * \strong{rShukaVar} The rank for \strong{ShukaVar}.
#' * \strong{ssiShukaVar} The simultaneous selection index (\eqn{ssiShukaVar = rMean + rShukaVar}).
#' @md
#' @export
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @references Shukla, G.K. 1972. Some statistical aspects of partitioning
#'   genotype-environmental components of variability. Heredity. 29:238-245.
#'   \href{http://www.nature.com/articles/hdy197287}{doi:10.1038/hdy.1972.87}.
#' @references Kang, M.S., and H.N. Pham. 1991. Simultaneous Selection for High
#'   Yielding and Stable Crop Genotypes. Agron. J. 83:161.
#'   \href{https://acsess.onlinelibrary.wiley.com/doi/abs/10.2134/agronj1991.00021962008300010037x}{doi:10.2134/agronj1991.00021962008300010037x}.
#'
#' @examples
#' \donttest{
#' library(metan)
#'out <- Shukla(data_ge2,
#'              env = ENV,
#'              gen = GEN,
#'              rep = REP,
#'              resp = PH)
#'}
Shukla <- function(.data, env, gen, rep, resp, verbose = TRUE) {
  factors  <-
    .data %>%
    select({{env}}, {{gen}}, {{rep}}) %>%
    mutate(across(everything(), as.factor))
  vars <-
    .data %>%
    select({{resp}}, -names(factors)) %>%
    select_numeric_cols()
  factors %<>% set_names("ENV", "GEN", "REP")
  g <- nlevels(factors$GEN)
  e <- nlevels(factors$ENV)
  r <- nlevels(factors$REP)
  listres <- list()
  nvar <- ncol(vars)
  if (verbose == TRUE) {
    pb <- progress_bar$new(
      format = "Evaluating the variable :what [:bar]:percent",
      clear = FALSE, total = nvar, width = 90)
  }
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
    g_means <- means_by(data, GEN)
    ge_means <- means_by(data, GEN, ENV)
    ge_effect <- ge_means %>%
      mutate(ge = residuals(lm(Y ~ ENV + GEN, data = .))) %>%
      make_mat(GEN, ENV, ge) %>%
      as.matrix()
    Wi <- rowSums(ge_effect^2)
    ShuklaVar <- (g * (g - 1) * Wi - sum(Wi)) / ((e - 1) * (g - 1) * ( g - 2))
    temp <- as_tibble(cbind(g_means, ShuklaVar)) %>%
      mutate(rMean = rank(-Y),
             rShukaVar = rank(ShuklaVar),
             ssiShukaVar = rMean + rShukaVar)
    if (verbose == TRUE) {
      pb$tick(tokens = list(what = names(vars[var])))
    }
    listres[[paste(names(vars[var]))]] <- temp
  }
  return(structure(listres, class = "Shukla"))
}






#' Print an object of class Shukla
#'
#' Print the \code{Shukla} object in two ways. By default, the results
#' are shown in the R console. The results can also be exported to the directory
#' into a *.txt file.
#'
#'
#' @param x The \code{Shukla} x
#' @param export A logical argument. If \code{TRUE}, a *.txt file is exported to
#'   the working directory.
#' @param file.name The name of the file if \code{export = TRUE}
#' @param digits The significant digits to be shown.
#' @param ... Options used by the tibble package to format the output. See
#'   \code{\link[tibble:formatting]{tibble::print()}} for more details.
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @method print Shukla
#' @export
#' @examples
#' \donttest{
#' library(metan)
#' eco <- Shukla(data_ge2,
#'   env = ENV,
#'   gen = GEN,
#'   rep = REP,
#'   resp = PH
#' )
#' print(eco)
#' }
print.Shukla <- function(x, export = FALSE, file.name = NULL, digits = 3, ...) {
  if (!class(x) == "Shukla") {
    stop("The object must be of class 'Shukla'")
  }
  opar <- options(pillar.sigfig = digits)
  on.exit(options(opar))
  if (export == TRUE) {
    file.name <- ifelse(is.null(file.name) == TRUE, "Shukla print", file.name)
    sink(paste0(file.name, ".txt"))
  }
  for (i in 1:length(x)) {
    var <- x[[i]]
    cat("Variable", names(x)[i], "\n")
    cat("---------------------------------------------------------------------------\n")
    cat("Shukla stability variance\n")
    cat("---------------------------------------------------------------------------\n")
    print(var)
  }
  cat("\n\n\n")
  if (export == TRUE) {
    sink()
  }
}
