find_outliers = function(.data,
                         var,
                         group_var = NULL,
                         plots = FALSE,
                         coef = 1.5){
internal = function(x){
    var_name = x %>% select(!!enquo(var)) %>% unlist() %>% as.numeric()
    tot <- sum(!is.na(var_name))
    na1 <- sum(is.na(var_name))
    m1 <- mean(var_name, na.rm = T)
    m11 <- (sd(var_name, na.rm = T)/m1)*100
    if(plots == TRUE){
    par(mfrow=c(2, 2), oma=c(0,0,0,0))
    boxplot(var_name, main = "With outliers")
    hist(var_name, main = "Without outliers", xlab=NA, ylab=NA)
    }
    outlier <- boxplot.stats(var_name, coef =  coef)$out
    dd = .data %>% select(!!enquo(var))
    names_out = paste(which(dd[,1] %in% outlier), sep = " ")
    mo <- mean(outlier)
    var_name <- ifelse(var_name %in% outlier, NA, var_name)
    names = rownames(var_name)
    if(plots == TRUE){
    boxplot(var_name, main = "With outliers")
    hist(var_name, main = "Without outliers", xlab=NA, ylab=NA)
    }
    na2 <- sum(is.na(var_name))
    if((na2 - na1)>0){
      cat("Number of possible outliers:", na2 - na1, "\n")
      cat("Lines:", names_out, "\n")
      cat("Proportion: ", round((na2 - na1) / sum(!is.na(var_name))*100, 1), "%\n", sep = "")
      cat("Mean of the outliers:", round(mo, 3), "\n")
      m2 <- mean(var_name, na.rm = T)
      m22 <- (sd(var_name, na.rm = T)/m2)*100
      cat("With outliers:    mean = ", round(m1, 3), " | CV = ",round(m11, 3),"%", sep = "","\n")
      cat("Without outliers: mean = ", round(m2, 3), " | CV = ",round(m22, 3),"%", sep = "","\n")
    }
    if((na2 - na1)==0){
      cat("No outlier identified. \n")
    }
}

if(!missing(group_var)){
    group_var <- dplyr::enquo(group_var)
    re =  .data  %>%
      split(dplyr::pull(., !!group_var))
    data = lapply(re, function(x){
        message("The factors ", paste0(collapse = " ", names(x[ , unlist(lapply(x, is.factor)) ])),
                " where excluded from the data")
      x[ , unlist(lapply(x, is.numeric))]
    }
    )
    out = lapply(seq_along(data), function(x){
        cat("\n----------------------------------------------------------------------------\n")
        cat("Level", names(data)[[x]], "\n")
        cat("----------------------------------------------------------------------------\n")
      internal(data[[x]])
    })
    names(out) = names(data)
  } else{
    if(sum(lapply(.data, is.factor)==TRUE)>0){
        message("The factors ", paste0(collapse = " ", names(.data[ , unlist(lapply(.data, is.factor)) ])),
                " where excluded from the data")
    }
    data = .data[ , unlist(lapply(.data, is.numeric))]
    out = internal(data)
  }
}
