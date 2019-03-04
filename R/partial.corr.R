partial.corr <- function(data, n = NULL, method = "pearson") {
    
    if (is.matrix(data) == TRUE & is.null(n) == TRUE) {
        stop("You have a matrix as entire but do not declared the sample size used to estimate the matrix")
    }
    
    if (is.matrix(data) != TRUE & is.null(n) != TRUE) {
        stop("You have a dataframe as entire but declared the sample size used. When a dataframe is used, the sample size is automatically caltulated.")
    }
    if (is.null(n) == TRUE) {
        n <- nrow(data)
    } else {
        n <- n
    }
    nvar <- ncol(data)
    df <- n - nvar
    
    cl <- match.call()
    if (!is.matrix(data)) {
        data <- cor(data, use = "complete", method = "pearson")
    }
    m <- as.matrix(data)
    X.resid <- -(solve(m))
    diag(X.resid) <- 1/(1 - (1 - 1/diag(solve(m))))
    X.resid <- cov2cor(X.resid)
    
    results <- data.frame(linear = as.vector(t(m)[lower.tri(m, diag = F)]))
    results <- dplyr::mutate(results, partial = as.vector(t(X.resid)[lower.tri(X.resid, 
        diag = F)]), t = partial/(sqrt(1 - partial^2)) * sqrt(n - nvar), prob = 2 * 
        (1 - pt(abs(t), df = df)))
    names <- colnames(data)
    combnam <- combn(names, 2, paste, collapse = " x ")
    rownames(results) <- names(sapply(combnam, names))
    
    return(list(linear.mat = data.frame(m), partial.mat = data.frame(X.resid), results = results, 
        call = cl))
    
}
