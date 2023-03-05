#' Random Sampling
#'
#' @description
#'
#' * [sample_random()] performs Simple Random Sampling or Stratified Random
#' Sampling
#' * [sample_systematic()] performs systematic sampling. In this case, a regular
#' interval of size k (`k = floor(N/n)`) is generated considering the population
#' size (N) and desired sample size (n). Then, the starting member (`r`) is
#' randomly chosen between `1-k`. The second element is `r` + `k`, and so on.
#'
#' @param data A data frame. If `data` is a `grouped_df`, the operation will be
#'   performed on each group (stratified).
#' @param n,prop Provide either `n`, the number of rows, or `prop`, the
#'   proportion of rows to select. If neither are supplied, `n = 1` will be
#'   used.
#' @param r The starting element. By default, `r` is randomly selected between
#'   `1:k`
#' @param by A categorical variable to compute the sample by. It is a
#'   shortcut to [dplyr::group_by()] that allows to group the data by one
#'   categorical variable. If more than one grouping variable needs to be used,
#'   use [dplyr::group_by()] to pass the data grouped.
#' @param weight Sampling weights. This must evaluate to a vector of
#'   non-negative numbers the same length as the input. Weights are
#'   automatically standardised to sum to 1.
#'
#' @return An object of the same type as `data`.
#' @export
#' @name utils_samples
#'
#' @examples
#' library(metan)
#' sample_random(data_ge, n = 5)
#' sample_random(data_ge,
#'               n = 3,
#'               by = ENV)
#'
#' sample_systematic(data_g, n = 6)
sample_random <- function(data,
                          n,
                          prop,
                          by = NULL,
                          weight = NULL){
  if(!missing(by)){
    data <- data |> group_by({{by}})
  }
  dplyr::slice_sample(data, n = n, prop = prop, weight_by = weight) |>
    ungroup()
}

#' @export
#' @name utils_samples
sample_systematic <- function(data,
                              n,
                              r = NULL,
                              by =  NULL){
  aux <- function(data, n, r = NULL){
    k <- floor(nrow(data) / n)
    message("k = ", k)
    if(is.null(r)){
      r <- sample(1:k, 1)
    }
    if(r == 1){
      rows <- sample(1:nrow(data), n)
    } else{
      rows <- seq(r, r + k*(n-1), k)
    }
    slice(data, rows) |>
      mutate(.id = rows, .before = 1)
  }
  if(!missing(by)){
    data <- data |> group_by({{by}})
  }
  if(is_grouped_df(data)){
    groups <- group_vars(data)
    data |>
      ungroup() |>
      nest(data = -c(!!!syms(groups))) |>
      mutate(sample = map(data, ~.x |>  aux(n = n, r = r))) |>
      select(-data) |>
      unnest(sample)
  } else{
    aux(data, n = n, r = r)
  }
}


