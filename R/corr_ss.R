#' Sample size planning for a desired Pearson's correlation confidence interval
#' @description
#' `r badge('stable')`
#'
#' Find the required (sufficient) sample size for computing a Pearson
#' correlation coefficient with a desired confidence interval (Olivoto et al.,
#' 2018) as follows
#'\loadmathjax
#'\mjsdeqn{n = {\left[ {\frac{{C{I_w}}}{{{{0.45304}^r} \times 2.25152}}} \right]^{{\rm{ - 0}}{\rm{.50089}}}}}
#'
#' where \mjseqn{CI_w} is desired confidence interval and \mjseqn{r} is the
#' correlation coefficient.
#'
#' @param r The magnitude of the correlation coefficient.
#' @param CI The half-width for confidence interval at p < 0.05.
#' @param verbose Logical argument. If `verbose = FALSE` the code is run
#'   silently.
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @references Olivoto, T., A.D.C. Lucio, V.Q. Souza, M. Nardino, M.I. Diel,
#'   B.G. Sari, D.. K. Krysczun, D. Meira, and C. Meier. 2018. Confidence
#'   interval width for Pearson's correlation coefficient: a
#'   Gaussian-independent estimator based on sample size and strength of
#'   association. Agron. J. 110:1-8.
#'   \doi{10.2134/agronj2016.04.0196}
#'
#' @export
#' @examples
#'
#' \donttest{
#' corr_ss(r = 0.60, CI = 0.1)
#' }
#'
#'
corr_ss <- function(r, CI, verbose = TRUE) {
    n <- round((CI/(0.45304^r * 2.25152))^(1/-0.50089), 0)
    if(verbose == TRUE){
    cat("-------------------------------------------------", "\n")
    cat("Sample size planning for correlation coefficient", "\n")
    cat("-------------------------------------------------", "\n")
    cat(paste0("Level of significance: 5%", "\nCorrelation coefficient: ", r, "\n95% half-width CI: ",
        CI, "\nRequired sample size: ", n, "\n"))
    cat("-------------------------------------------------", "\n")
    }
    invisible(tibble(`Description` = c("Significance level (%)", "Correlation", "95% half-width CI", "Sample size"),
                  `Value` = c(95, r, CI, n))
    )
}
