#' Computes the coincidence index of genotype selection
#'
#' @description
#' Computes the coincidence index (Hamblin and Zimmermann, 1986) as follows:
#'
#' \loadmathjax
#' \mjsdeqn{CI = \frac{A-C}{M-C}\times 100}
#' where \emph{A} is the number of selected genotypes common to different
#' methods; \emph{C} is the number of expected genotypes selected by chance; and
#' \emph{M} is the number of genotypes selected according to the selection
#' intensity.
#'
#' @param ... A comma-separated list of objects of class \code{mgidi},
#'   \code{mtsi} \code{fai_blup}, or \code{sh}. When a model is informed, then
#'   the selected genotypes are extracted automatically.
#' @param total The total number of genotypes in the study.
#' @param sel1,sel2 The selected genotypes by the method 1 and 2, respectively. Defaults to \code{NULL}.
#' @return A list with the following elements:
#'
#' * \strong{coincidence}: A data frame with the coincidence index, number of
#' common genotypes and the list of common genotypes for each model combination.
#' * \strong{coincidence_mat}: A matrix-like containing the coincidence index.
#' * \strong{genotypes}: The number of common genotypes for all models, i.e.,
#' the insersection of the selected genotypes of all models
#'
#' @export
#' @md
#' @references
#' Hamblin, J., and M.J. de O. Zimmermann. 1986. Breeding Common Bean for Yield
#' in Mixtures. p. 245-272. In Plant Breeding Reviews. John Wiley & Sons, Inc.,
#' Hoboken, NJ, USA.\doi{10.1002/9781118061015.ch8}
#'
#' @examples
#' \donttest{
#' sel1 <- paste("G", 1:30, sep = "")
#' sel2 <- paste("G", 16:45, sep = "")
#' coincidence_index(sel1 = sel1, sel2 = sel2, total = 150)
#' }
coincidence_index <- function(..., total, sel1 = NULL, sel2 = NULL){
  comp_coinc <- function(sel1, sel2, total){
    common <- intersect(sel1, sel2)
    si <- round(length(sel1) / total, 2)
    rand <- length(sel1) * si
    coinc_idex <- ((length(common) - rand) / (length(sel1) - rand)) * 100
    return(coinc_idex)
  }
  if(!missing(...)){
    names <-
      as.character(
        sapply(quos(...), function(x){
          rlang::quo_get_expr(x)
        })
      )
    if(length(names)  < 2){
      stop("The coincidence index cannot be computed with only one model.", call. = FALSE)
    }
    models <- list(...)
    if(any(sapply(models, class) %in% c("mgidi", "fai_blup", "sh", "mtsi") == FALSE)){
      stop("Only objects of class 'mgidi', 'mtsi', 'fai_blup', and 'sh' are accepted.")
    }
    selected <- lapply(models, function(x){
      x[["sel_gen"]]
    })
    names(selected) <- names
    ngsel <- sapply(selected, length)
    if(length(unique(as.numeric(ngsel))) > 1){
      stop("The number of selected genotypes must be the same for all the models\n",
           paste(capture.output(print(data.frame(ngsel))), collapse = "\n"), call. = FALSE)
    }
    index <- combn(length(selected), 2)
    ncomb <- ncol(index)
    results <-
      data.frame(Model = combn(names(selected), 2, paste, collapse = "_::_")) %>%
      separate(Model, into = c("V1", "V2"), sep = "_::_")
    values <- NULL
    common <- NULL
    genotypes <- NULL
    common_gen_pairs <- NULL
    for (i in 1:ncomb) {
      V1 <- selected[[index[1, i]]]
      V2 <- selected[[index[2, i]]]
      values[i] <- comp_coinc(V1, V2, total)
      common[i] <- length(intersect(V1, V2))
      genotypes[i] <- list(intersect(sel1, sel2))
      common_gen_pairs[i] <- paste(intersect(V1, V2), collapse = ",")
    }
    results <-
      add_cols(results,
               index = values,
               common = common,
               genotypes = common_gen_pairs)
    final <-
      list(coincidence = results,
           coincidence_mat = make_mat(results, V1, V2, index),
           genotypes = Reduce(intersect, genotypes))
    return(final %>% set_class("coincidence"))
  } else{
    if(length(sel1) != length(sel2)){
      stop("The lenght of 'sel1' and 'sel2' must be equal")
    }
    return(comp_coinc(sel1, sel2, total))
  }
}




#' Print an object of class coincidence
#'
#' Print a \code{coincidence} object in two ways. By default, the results are shown in
#' the R console. The results can also be exported to the directory.
#'
#'
#' @param x An object of class \code{coincidence}.
#' @param export A logical argument. If \code{TRUE}, a *.txt file is exported
#'   to the working directory
#' @param file.name The name of the file if \code{export = TRUE}
#' @param digits The significant digits to be shown.
#' @param ... Options used by the tibble package to format the output. See
#'   \code{\link[tibble:formatting]{tibble::print()}} for more details.
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @method print coincidence
#' @export
#' @examples
#'\donttest{
#' library(metan)
#' sel1 <- paste("G", 1:30, sep = "")
#' sel2 <- paste("G", 16:45, sep = "")
#' coinc <- coincidence_index(sel1 = sel1, sel2 = sel2, total = 150)
#' print(coinc)
#' }
print.coincidence <- function(x, export = FALSE, file.name = NULL, digits = 4, ...) {
  if (!class(x) == "coincidence") {
    stop("The object must be of class 'coincidence'")
  }
  if (export == TRUE) {
    file.name <- ifelse(is.null(file.name) == TRUE, "waasb print", file.name)
    sink(paste0(file.name, ".txt"))
  }
  opar <- options(pillar.sigfig = digits)
  on.exit(options(opar))
  cat("---------------------------------------------------------------------------\n")
  cat("Coincidence index and common genotypes\n")
  cat("---------------------------------------------------------------------------\n")
  print(x$coincidence)
  if (export == TRUE) {
    sink()
  }
}
