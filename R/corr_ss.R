corr_ss <- function(r, CI) {
    n <- (CI/(0.45304^r * 2.25152))^(1/-0.50089)
    cat("-------------------------------------------------", "\n")
    cat("Sample size planning for correlation coefficient", "\n")
    cat("-------------------------------------------------", "\n")
    cat(paste0("Level of significance: 5%", "\nCorrelation coefficient: ", r, "\n95% half-width CI: ",
        CI, "\nRequired sample size: ", round(n, 0)), "\n")
    cat("-------------------------------------------------", "\n")
}
