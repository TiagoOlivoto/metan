#' Check for common erros in multi-environment trial data
#'
#' \code{inspect()} scans a data.frame object for errors that may affect the use
#' of functions in \code{metan}. By default, all variables are checked regarding
#' the class (numeric or factor), missing values, and presence of possible
#' outliers. The function will return a warning if the data looks like
#' unbalanced, has missing values or possible outliers.
#'
#' @param .data The data to be analyzed
#' @param ... The variables in \code{.data} to check. Set to
#'   \code{NULL}, i.e., all the variables in \code{.data} are used.
#' @param plot Create a plot to show the check? Defaults to \code{FALSE}.
#' @param verbose Logical argument. If \code{TRUE} (default) then the results
#'   for checks are shown in the console.
#'
#' @return A tibble with the following variables:
#' * \strong{Variable} The name of variable
#' * \strong{Class} The class of the variable
#' * \strong{Missing} Contains missing values?
#' * \strong{Levels} The number of levels of a factor variable
#' * \strong{Valid_n} Number of valid n (omit NAs)
#' * \strong{Outlier} Contains possible outliers?
#' @md
#' @export
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#'
#' @examples
#' library(metan)
#' inspect(data_ge)
#'
#' # Create a toy example with messy data
#' df <- data_ge2[-c(2, 30, 45, 134), c(1:6)]
#' df[c(1, 20, 50), c(4, 6)] <- NA
#' df[40, 6] <- df[40, 6] * 2
#'
#' inspect(df, plot = TRUE)
inspect <- function (.data,
                     ...,
                     plot = FALSE,
                     verbose = TRUE) {
  if(!missing(...)){
    .data <- dplyr::select(.data, ...)
  } else{
    .data <- .data
  }
  df <-
    data.frame(
      Class = sapply(.data, class),
      Missing= sapply(.data, function(x) ifelse(any(is.na(x)), "Yes", "No")),
      Levels = sapply(.data, function(x) ifelse(!is.numeric(x), nlevels(x), "-")),
      Valid_n = sapply(.data, function(x) length(which(!is.na(x)))),
      Min = sapply(.data, function(x) ifelse(is.numeric(x), round(min(x),2), NA)),
      Median = sapply(.data, function(x) ifelse(is.numeric(x), round(median(x),2), NA)),
      Max = sapply(.data, function(x) ifelse(is.numeric(x), round(max(x),2), NA)),
      Outlier = sapply(.data, function(x) ifelse(is.numeric(x), find_outliers(values = x, verbose = F), NA))
    ) %>%
    rownames_to_column("Variable") %>%
    as_tibble()
  esp_nrows <- prod(as.numeric(levels(droplevels(df[which(df$Levels != "-"),][4])$Levels)))
  if(verbose == TRUE){
    print(df)
    nfactors <- sum(lapply(.data, is.factor) == TRUE)
    if(esp_nrows != nrow(.data)){
      warning("Considering the levels of factors, .data should have ", esp_nrows, " rows, but it has ", nrow(.data), call. = F)
    }
    if (nfactors < 3){
      warning("Expected three or more factor variables. The data has only ", nfactors, ": ", paste(collapse = " ", names(.data[, unlist(lapply(.data, is.factor))])), call. = F)
    }
    if(any(df$Missing == "Yes")){
      warning("Missing values in variable(s) ", paste(df$Variable[c(which(df$Missing == "Yes"))], collapse = " "), call. = F)
    }
    if(any(df$Outlier[!is.na(df$Outlier)] != 0)){
      warning("Possible outliers in variable(s) ", paste(df$Variable[c(which(df$Outlier != 0))], collapse = " "),". Use the function find_outliers() for more details.", call. = F)
    }
    if(nfactors >= 3 && all(df$Missing == "No") && all(df$Outlier[!is.na(df$Outlier)] == 0) == TRUE){
      message("No issues detected while inspecting data.")
    }
  }
  if(plot == TRUE){
    ggpair <-
      .data %>%
      ggpairs(lower = NULL,
              upper = list(continuous ="smooth"))+
      theme(panel.spacing = unit(0.1, "cm"))

    suppressMessages(suppressWarnings(print(ggpair, progress = FALSE)))
  }
  invisible(df)
}

