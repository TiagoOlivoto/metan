lpcor <- function(.data, group_var = NULL, n = NULL, verbose = TRUE) {

  if(!is.matrix(.data) && !is.data.frame(.data)){
    stop("The object 'x' must be a correlation matrix or a data.frame")
  }
  if(is.matrix(.data) && is.null(n)){
    stop("You have a matrix but the sample size used to compute the correlations (n) was not declared.")
  }
  if(is.matrix(.data) && !missing(group_var)){
    stop("You cannot use a grouping variables because a correlation matrix was used as input.")
  }
  if(is.data.frame(.data) && !is.null(n)){
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
      nvar <- ncol(cor.x)
      df <- n - nvar
      if(df < 0){
        warning("The number of variables is higher than the number of individuals. Hypothesis testing will not be made.",
                call. = FALSE)
      }
    m <- as.matrix(cor.x)
    X.resid <- -(solve(m))
    diag(X.resid) <- 1/(1 - (1 - 1/diag(solve(m))))
    X.resid <- cov2cor(X.resid)
    results <- data.frame(linear = as.vector(t(m)[lower.tri(m, diag = F)]))
    results <- suppressWarnings(dplyr::mutate(results,
                             partial = as.vector(t(X.resid)[lower.tri(X.resid, diag = F)]),
                             t = partial/(sqrt(1 - partial^2)) * sqrt(n - nvar), prob = 2 *
                               (1 - pt(abs(t), df = df))))
    names <- colnames(x)
    combnam <- combn(names, 2, paste, collapse = " x ")
    rownames(results) <- names(sapply(combnam, names))
    if(verbose == TRUE){
      print.data.frame(results)
    }
    return(list(linear.mat = data.frame(m),
                partial.mat = data.frame(X.resid),
                results = results))
    }

    if (is.matrix(.data)) {
       out = internal(.data)
    }

    if(is.data.frame(.data)){
      if(!missing(group_var)){
        group_var <- dplyr::enquo(group_var)
        re =  .data  %>%
          split(dplyr::pull(., !!group_var))
        data = lapply(re, function(x){
          if(verbose == TRUE){
            message("The factors ", paste0(collapse = " ", names(x[ , unlist(lapply(x, is.factor)) ])),
                    " where excluded from the data")
          }
          x[ , unlist(lapply(x, is.numeric))]
        }
        )
        out = lapply(seq_along(data), function(x){
          if(verbose == TRUE){
            cat("\n----------------------------------------------------------------------------\n")
            cat("Level:", names(data)[[x]], "\n")
            cat("----------------------------------------------------------------------------\n")
          }
          internal(data[[x]])
        })
        names(out) = names(data)
      } else{
        if(sum(lapply(.data, is.factor)==TRUE)>0){
          if(verbose == TRUE){
            message("The factors ", paste0(collapse = " ", names(.data[ , unlist(lapply(.data, is.factor)) ])),
                    " where excluded from the data")
          }
        }
        data = .data[ , unlist(lapply(.data, is.numeric))]
        out = internal(data)
      }
    }
    if(!missing(group_var)){
    invisible(structure(out, class = "lpcor_group"))
    } else {
    invisible(structure(out, class = "lpcor"))
    }
}
