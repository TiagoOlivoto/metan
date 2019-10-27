#' Print an object of class anova_joint
#'
#' Print the \code{anova_joint} object in two ways. By default, the results are shown
#' in the R console. The results can also be exported to the directory into a
#' *.txt file.
#'
#' @param x An object of class \code{anova_joint}.
#' @param export A logical argument. If \code{TRUE}, a *.txt file is exported to
#'   the working directory.
#' @param file.name The name of the file if \code{export = TRUE}
#' @param digits The significant digits to be shown.
#' @param ... Options used by the tibble package to format the output. See
#'   \code{\link[tibble]{trunc_mat}} for more details.
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @method print anova_joint
#' @export
#' @examples
#'
#' library(metan)
#' model <- data_ge %>% anova_joint(ENV, GEN, REP, c(GY, HM))
#' print(model)
print.anova_joint <- function(x, export = FALSE, file.name = NULL, digits = 3, ...) {
  if (!class(x) == "anova_joint") {
    stop("The object must be of class 'anova_joint'")
  }
  backup_options <- options()
  options(pillar.sigfig = digits, ...)
  if (export == TRUE) {
    file.name <- ifelse(is.null(file.name) == TRUE, "anova_joint print", file.name)
    sink(paste0(file.name, ".txt"))
  }
  for (i in 1:length(x)) {
    var <- x[[i]]
    cat("Variable", names(x)[i], "\n")
    cat("---------------------------------------------------------------------------\n")
    print(var)
    cat("---------------------------------------------------------------------------\n")
    cat("\n\n\n")
  }
  if (export == TRUE) {
    sink()
  }
  options(backup_options)
}
