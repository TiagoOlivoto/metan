#' Transpose a data frame
#'
#' Is an alternative to \code{\link{t}()} to transpose a data frame. The first
#' column of \code{df} will become column names in the transposed data.
#' @param df A data frame to be transposed.
#'
#' @return A tibble containing the transposed data.
#' @export
#'
#' @examples
#' \donttest{
#' library(metan)
#' df <-
#'data.frame(
#'  GEN = c("G1", "G2", "G3","G4"),
#'  E1 = rnorm(4, 100, 20),
#'  E2 = rnorm(4, 10, 2),
#'  E3 = rnorm(4, 50, 5),
#'  E4 = rnorm(4, 1000, 150)
#')
#'df
#'t(df)
#'transpose_df(df)
#' }
transpose_df <- function(df){
  df %>%
    pivot_longer(-1) %>%
    pivot_wider(names_from = 1, values_from = 3)
}
