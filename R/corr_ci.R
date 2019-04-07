corr_ci <- function(.data = NA, r = NULL, n = NULL, verbose = TRUE) {
  if (any(!is.na(.data))) {
    if(!is.matrix(.data) && !is.data.frame(.data) && !is.group_factors(.data)){
      stop("The object 'x' must be a correlation matrix, a data.frame or an object of class group_factors")
    }

    if(is.matrix(.data) && is.null(n)){
      stop("You have a matrix but the sample size used to compute the correlations (n) was not declared.")
    }


    if(is.data.frame(.data) && !is.null(n)){
      stop("You cannot informe the sample size because a data frame was used as input.")
    }

    if(is.group_factors(.data) && !is.null(n)){
      stop("You cannot informe the sample size because a data frame was used as input.")
    }


      internal = function(x){
        if(is.matrix(x)){
          cor.x = x
        }
        if(is.data.frame(x)){
          cor.x = cor(x)
          n = nrow(x)
        }
        m <- as.matrix(cor.x)

        results <- data.frame(Corr = as.vector(t(m)[lower.tri(m, diag = F)]))
        CI <- dplyr::mutate(results, CI = (0.45304^abs(Corr)) * 2.25152 * (n^-0.50089),
            LL = Corr - CI, UL = Corr + CI)
        combnam <- combn(colnames(x), 2, paste, collapse = " x ")
        rownames(CI) <- names(sapply(combnam, names))

        return(CI)
    }

      if (is.matrix(.data)){
        out = internal(.data)
      }

      if (any(class(.data) == "group_factors")) {
        out = lapply(seq_along(.data), function(x){
          internal(.data[[x]])
        })
        names(out) = names(.data)
      }

      if (is.data.frame(.data)) {
        if (sum(lapply(.data, is.factor)==TRUE)>0){
        }
        data = data.frame(.data[ , unlist(lapply(.data, is.numeric))])
        out = internal(data)
        if (verbose == TRUE){
          message("The factors ", paste0(collapse = " ", names(.data[ , unlist(lapply(.data, is.factor)) ])),
                  " where excluded to perform the analysis. If you want to perform an analysis for each level of a factor, use the function 'group_factors() before.' ")
        }
      }

invisible(out)

    } else {
        CI <- (0.45304^r) * 2.25152 * (n^-0.50089)
        UP <- r + CI
        LP <- r - CI
        cat("-------------------------------------------------", "\n")
        cat("Nonparametric 95% half-width confidence interval", "\n")
        cat("-------------------------------------------------", "\n")
        cat(paste0("Level of significance: 5%", "\n", "Correlation coefficient: ",
            r, "\nSample size: ", n, "\nConfidence interval: ", round(CI, 4), "\nTrue parameter range from: ",
            round(LP, 4), " to ", round(UP, 4)), "\n")
        cat("-------------------------------------------------", "\n")
    }
}

