#' @title Utilities for text progress bar in the terminal
#'
#' @description
#' `r badge('experimental')`
#'
#' Progress bars are configurable, may include percentage, elapsed time, and
#' custom text.
#' * [progress()]: Initiate a custom progress bar of class `pb_metan`.
#' * [run_progress()]: Run the progress bar and should be called within a 'for
#' loop' statement, a [lapply()] family or [purrr::map()] family of functional
#' programming tools.
#' @name utils_progress
#' @param min,max Numeric values for the extremes of the progress bar. Must have
#'   `min < max`.
#' @param leftd,rightd The left and right delimiters for the progress bar.
#'   Defaults to `"|"`.
#' @param char The character (or character string) to form the progress bar.
#' @param style The 'style' of the progress bar. Elapsed time is counted from
#'   calling [progress()] up to each call of [run_progress()].
#'    - `type = 1`: Shows a progress bar without percentage or elapsed time.
#'    - `type = 2`: The default, shows the progress bar and its percentage.
#'    - `type = 3`: Shows the progress bar and elapsed time.
#'    - `type = 4`: Shows the progress bar, percentage, and elapsed time.
#' @param width The the width of the progress bar. Defaults to the number of
#'   characters is that which fits into `getOption("width")`.
#' @param time The system time used to compute the elapsed time from calling
#'   [progress()] to each call of [run_progress()]. Defaults to `Sys.time()`.
#' @return [progress()] returns a list of class `pb_metan` that contains the set
#'   parameters that will called by [run_progress()].
#' @param pb An object created with [progress()]
#' @param actual The actual value, for example, a loop variable that define the
#'   loop index value.
#' @param text An optional character string to be shown at the begining of the
#'   progress bar.
#' @param digits The number of significant figures in percentage value. Defaults
#'   to `0`.
#' @param sleep Suspend execution for a time interval with [Sys.sleep()] within
#'   [run_progress()]. Defaults to `0`.
#' @md
#' @export
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @examples
#' \donttest{
#'
#' ################### A for looping approach ################
#' pb <-  progress()
#' for (i in 1:100) {
#'  run_progress(pb, actual = i, sleep = 0.01)
#' }
#'
#' ################### An apply family approach ##############
#' pb <- progress(max = 10)
#' foo <- function(...){
#'    run_progress(pb, ...)
#'    rnorm(100) %>%  mean()
#'  }
#' (a <- sapply(1:10, FUN = foo, sleep = 0.05))
#'
#' ######## A purrr functional programming approach ##########
#' foo2 <- function(...){
#'      run_progress(pb2, ...)
#'      rnorm(100) %>%  mean()
#' }
#' pb2 <- progress(max = 10000,
#'                 style = 4,
#'                 leftd = "",
#'                 char = ".",
#'                 rightd = "!")
#'
#' b <- purrr::map_dbl(1:10000, foo2, text = "Progress bar for sampling")
#' hist(b)
#' }
#
progress <- function(min = 0,
                     max = 100,
                     leftd = "|",
                     rightd = "|",
                     char = "=",
                     style = 2,
                     width = getOption("width"),
                     time = Sys.time()){
  # Adapted from https://stackoverflow.com/a/26920123/15245107
  return(list(min = min,
              max = max,
              leftd = leftd,
              rightd = rightd,
              char = char,
              style = style,
              width = width,
              time = time) %>%
           set_class("pb_metan"))
}
#' @name utils_progress
#' @export
run_progress <- function(pb,
                         actual,
                         text = "",
                         digits = 0,
                         sleep = 0){
  Sys.sleep(sleep)
  elapsed <-
    as.numeric(difftime(Sys.time(), pb$time, units = "secs")) %>%
    sec_to_hms()
  temp <- switch(
    pb$style,
    list(extra = nchar(text) + nchar(pb$leftd) + nchar(pb$rightd),
         text = paste(text, paste(pb$leftd, '%s%s', pb$right, sep = ""))),
    list(extra = nchar(text) + nchar(pb$leftd) + nchar(pb$rightd) + 6,
         text =  paste(text, paste(pb$leftd, '%s%s', pb$right, sep = ""), '% s%%')),
    list(extra = nchar(text) + nchar(pb$leftd) + nchar(pb$rightd) + 9,
         text = paste(text, paste(pb$leftd, '%s%s', pb$rightd, sep = ""), elapsed)),
    list(extra = nchar(text) + nchar(pb$leftd) + nchar(pb$rightd) + 15,
         text = paste(text, paste(pb$leftd, '%s%s', pb$rightd, sep = ""), '% s%%', elapsed))
  )
  step <- round(actual / pb$max * (pb$width - temp$extra))
  temp$text %>%
    sprintf(strrep(pb$char, step),
            strrep(' ', pb$width - step - temp$extra),
            round(actual / pb$max * 100, digits = digits)) %>%
    cat("\r")
  if(actual == pb$max){
    cat("\n")
  }
}
