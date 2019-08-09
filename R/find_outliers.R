#' Find possible outliers in a dataset
#'
#' Find possible outliers in the dataset.
#'
#'
#' @param .data The data to be analyzed. Must be a dataframe or an object of
#' class \code{split_factors}.
#' @param var The variable to be analyzed..
#' @param plots If \code{TRUE}, then histograms and boxplots are shown.
#' @param coef The multiplication coefficient, defaults to 1.5. For more details see
#' \code{?boxplot.stat}.
#' @param verbose If \code{verbose = TRUE} then some results are shown in the
#' console.
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @export
#' @examples
#'
#' library(metan)
#' library(dplyr)
#'
#' find_outliers(data_ge2, var = PH, plots = TRUE)
#'
#' data_ge2 %>%
#' split_factors(ENV) %>%
#' find_outliers(var = PH)
#'
#'
find_outliers <- function(.data, var, plots = FALSE, coef = 1.5,
                          verbose = TRUE) {
  if (!is.data.frame(.data) && !is.split_factors(.data)) {
    stop("The object 'x' must be a data.frame or an object of class split_factors")
  }
  if (missing(var) == TRUE) {
    stop("A variable must be declared.")
  }
  if (any(class(.data) == "split_factors")) {
    for (k in 1:length(.data)) {
      data <- .data[[k]]
      var <- dplyr::enquo(var)
      nam <- names(.data[k])
      var_name <- data %>% dplyr::select(!!var) %>% unlist() %>%
        as.numeric()
      tot <- sum(!is.na(var_name))
      na1 <- sum(is.na(var_name))
      m1 <- mean(var_name, na.rm = T)
      m11 <- (sd(var_name, na.rm = T)/m1) * 100
      if (any(plots == TRUE)) {
        op <- par(mfrow = c(2, 2), oma = c(0, 0, 2, 0),
                  mar = c(2, 2, 2, 2))
        boxplot(var_name, main = "With outliers")
        hist(var_name, main = "Without outliers", xlab = NA,
             ylab = NA)
      }
      outlier <- boxplot.stats(var_name, coef = coef)$out
      dd <- data.frame(data %>% select(!!var))
      names_out <- paste(which(dd[, 1] %in% outlier), sep = " ")
      if (length(outlier) >= 1) {
        mo <- mean(outlier)
        maxo <- max(outlier)
        mino <- min(outlier)
        names_out_max <- paste(which(dd[, 1] == maxo))
        names_out_min <- paste(which(dd[, 1] == mino))
      }
      var_name <- ifelse(var_name %in% outlier, NA, var_name)
      names <- rownames(var_name)
      if (any(plots == TRUE)) {
        boxplot(var_name, main = "With outliers")
        hist(var_name, main = "Without outliers", xlab = NA,
             ylab = NA)
        mtext(paste(nam), outer = TRUE, cex = 1.5)
        par(op)
      }
      na2 <- sum(is.na(var_name))
      if ((na2 - na1) > 0) {
        if (verbose == TRUE) {
          cat("\n----------------------------------------------------------------------------\n")
          cat("Level:", nam, "\n")
          cat("----------------------------------------------------------------------------\n")
          cat("Number of possible outliers:", na2 - na1,
              "\n")
          cat("Lines:", names_out, "\n")
          cat("Proportion: ", round((na2 - na1)/sum(!is.na(var_name)) *
                                      100, 1), "%\n", sep = "")
          cat("Mean of the outliers:", round(mo, 3),
              "\n")
          cat("Maximum of the outliers:", round(maxo,
                                                3), " | Line", names_out_max, "\n")
          cat("Minimum of the outliers:", round(mino,
                                                3), " | Line", names_out_min, "\n")
          m2 <- mean(var_name, na.rm = T)
          m22 <- (sd(var_name, na.rm = T)/m2) * 100
          cat("With outliers:    mean = ", round(m1,
                                                 3), " | CV = ", round(m11, 3), "%", sep = "",
              "\n")
          cat("Without outliers: mean = ", round(m2,
                                                 3), " | CV = ", round(m22, 3), "%", sep = "",
              "\n")
        }
      }
      if ((na2 - na1) == 0) {
        if (verbose == TRUE) {
          cat("\n----------------------------------------------------------------------------\n")
          cat("Level:", nam, "\n")
          cat("----------------------------------------------------------------------------\n")
          cat("No outlier identified. \n")
        }
      }
    }
  } else {
    if (sum(lapply(.data, is.factor) == TRUE) > 0) {
      if (verbose == TRUE) {
        message("The factors ", paste0(collapse = " ",
                                       names(.data[, unlist(lapply(.data, is.factor))])),
                " were ignored. Use 'split_factors()' before if you want to perform an analysis for each level of a factor.' ")
      }
    }
    var <- dplyr::enquo(var)
    var_name <- .data %>% dplyr::select(!!var) %>% unlist() %>%
      as.numeric()
    tot <- sum(!is.na(var_name))
    na1 <- sum(is.na(var_name))
    m1 <- mean(var_name, na.rm = T)
    m11 <- (sd(var_name, na.rm = T)/m1) * 100
    if (plots == TRUE) {
      op <- par(mfrow = c(2, 2), oma = c(0, 0, 2, 0), mar = c(2,
                                                              2, 2, 2))
      boxplot(var_name, main = "With outliers")
      hist(var_name, main = "With outliers", xlab = NA,
           ylab = NA)
    }
    outlier <- boxplot.stats(var_name, coef = coef)$out
    dd <- data.frame(.data %>% select(!!var))
    names_out <- paste(which(dd[, 1] %in% outlier), sep = " ")
    if (length(outlier) >= 1) {
      mo <- mean(outlier)
      maxo <- max(outlier)
      mino <- min(outlier)
      names_out_max <- paste(which(dd[, 1] == maxo))
      names_out_min <- paste(which(dd[, 1] == mino))
    }
    var_name <- ifelse(var_name %in% outlier, NA, var_name)
    names <- rownames(var_name)
    if (any(plots == TRUE)) {
      boxplot(var_name, main = "Without outliers")
      hist(var_name, main = "Without outliers", xlab = NA,
           ylab = NA)
      par(op)
    }
    na2 <- sum(is.na(var_name))
    if ((na2 - na1) > 0) {
      cat("Number of possible outliers:", na2 - na1, "\n")
      cat("Lines:", names_out, "\n")
      cat("Proportion: ", round((na2 - na1)/sum(!is.na(var_name)) *
                                  100, 1), "%\n", sep = "")
      cat("Mean of the outliers:", round(mo, 3), "\n")
      cat("Maximum of the outliers:", round(maxo, 3), " | Line",
          names_out_max, "\n")
      cat("Minimum of the outliers:", round(mino, 3), " | Line",
          names_out_min, "\n")
      m2 <- mean(var_name, na.rm = T)
      m22 <- (sd(var_name, na.rm = T)/m2) * 100
      cat("With outliers:    mean = ", round(m1, 3), " | CV = ",
          round(m11, 3), "%", sep = "", "\n")
      cat("Without outliers: mean = ", round(m2, 3), " | CV = ",
          round(m22, 3), "%", sep = "", "\n")
    }
    if ((na2 - na1) == 0) {
      cat("No possible outlier identified. \n")
    }
  }
}
