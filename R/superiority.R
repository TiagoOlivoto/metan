#' Lin e Binns' superiority index
#'
#' Nonparametric stability analysis using the superiority index proposed by Lin
#' & Binns (1988).
#'
#'
#' @param .data The dataset containing the columns related to Environments,
#'   Genotypes, replication/block and response variable(s)
#' @param env The name of the column that contains the levels of the
#'   environments.
#' @param gen The name of the column that contains the levels of the genotypes.
#' @param resp The response variable(s). To analyze multiple variables in a
#'   single procedure use, for example, \code{resp = c(var1, var2, var3)}.
#' @param verbose Logical argument. If \code{verbose = FALSE} the code will run
#'   silently.
#' @return An object of class \code{superiority} where each element is the
#'   result of one variable and contains the following items:
#'
#' * \strong{environments} The mean for each environment, the environment index
#' and classification as favorable and unfavorable environments.
#' * \strong{index} The superiority index computed for all (\code{Pi_a}),
#' favorable (\code{Pi_f}) and unfavorable (\code{Pi_u}) environments.
#'
#' @md
#' @author Tiago Olivoto, \email{tiagoolivoto@@gmail.com}
#' @seealso \code{\link{Annicchiarico}, \link{ecovalence}, \link{ge_stats}}
#' @references Lin, C.S., and M.R. Binns. 1988. A superiority measure of
#'   cultivar performance for cultivar x location data. Can. J. Plant Sci.
#'   68:193-198.
#'   \href{https://www.nrcresearchpress.com/doi/10.4141/cjps88-018#.XmEYwKhKgdU}{doi:10.4141/cjps88-018}
#'
#' @export
#' @examples
#' \donttest{
#' library(metan)
#' out <- superiority(data_ge2, ENV, GEN, PH)
#' print(out)
#'}
#'
superiority <- function(.data, env, gen, resp, verbose = TRUE) {
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
    environments <-
      data %>%
      means_by(ENV, na.rm = TRUE) %>%
      add_cols(index = Y - mean(Y),
               class = ifelse(index < 0, "unfavorable", "favorable")) %>%
      as_tibble()
    data <- left_join(data, environments %>% select(ENV, class), by = "ENV")
    lin_fun <- function(mat) {
      P <- apply(mat, 1, function(x) {
        sum((x - apply(mat, 2, max))^2)/(2 * length(x))
      })
      return(P)
    }
    mat_g <- make_mat(data, row = GEN, col = ENV, value = Y)
    ge_mf <- subset(data, class == "favorable")
    mat_f <- dplyr::select_if(make_mat(ge_mf, row = GEN,
                                       col = ENV, value = Y), function(x) !any(is.na(x)))
    ge_mu <- subset(data, class == "unfavorable")
    mat_u <- dplyr::select_if(make_mat(ge_mu, row = GEN,
                                       col = ENV, value = Y), function(x) !any(is.na(x)))
    temp <- list(environments = environments,
                 index = tibble(GEN = rownames(mat_g),
                                Y = apply(mat_g, 1, mean),
                                Pi_a = lin_fun(mat_g),
                                R_a = rank(lin_fun(mat_g)),
                                Pi_f = lin_fun(mat_f),
                                R_f = rank(lin_fun(mat_f)),
                                Pi_u = lin_fun(mat_u),
                                R_u = rank(lin_fun(mat_u))))
    rownames(temp) <- NULL
    if (verbose == TRUE) {
      pb$tick(tokens = list(what = names(vars[var])))
    }
    listres[[paste(names(vars[var]))]] <- temp
  }
  return(structure(listres, class = "superiority"))
}






#' Print an object ofclass \code{superiority}
#'
#' Print the \code{superiority} object in two ways. By default, the results are
#' shown in the R console. The results can also be exported to the directory
#' into a *.txt file.
#'
#'
#' @param x An object of class \code{superiority}.
#' @param export A logical argument. If \code{TRUE}, a *.txt file is exported to
#'   the working directory.
#' @param file.name The name of the file if \code{export = TRUE}
#' @param digits The significant digits to be shown.
#' @param ... Options used by the tibble package to format the output. See
#'   \code{\link[tibble:formatting]{tibble::print()}} for more details.
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @method print superiority
#' @export
#' @examples
#' \donttest{
#' library(metan)
#' model <- superiority(data_ge2, ENV, GEN, PH)
#' print(model)
#' }
print.superiority <- function(x, export = FALSE, file.name = NULL, digits = 3, ...) {
  if (!class(x) == "superiority") {
    stop("The object must be of class 'superiority'")
  }
  if (export == TRUE) {
    file.name <- ifelse(is.null(file.name) == TRUE, "superiority summary", file.name)
    sink(paste0(file.name, ".txt"))
  }
  opar <- options(pillar.sigfig = digits)
  on.exit(options(opar))
  for (i in 1:length(x)) {
    var <- x[[i]]
    cat("Variable", names(x)[i], "\n")
    cat("---------------------------------------------------------------------------\n")
    cat("Superiority index considering all, favorable and unfavorable environments\n")
    cat("---------------------------------------------------------------------------\n")
    print(var$index)
    cat("---------------------------------------------------------------------------\n")
    cat("\n\n\n")
  }
  if (export == TRUE) {
    sink()
  }
}
