find_outliers = function(.data,
                         var,
                         plots = FALSE,
                         coef = 1.5,
                         verbose = TRUE){
  if(!is.data.frame(.data) && !is.group_factors(.data)){
    stop("The object 'x' must be a data.frame or an object of class group_factors")
  }

  if (missing(var) == TRUE){
    stop("A variable must be declared.")
  }

  if(any(class(.data) == "group_factors")){
    for (k in 1:length(.data)){
      data = .data[[k]]
      var = dplyr::enquo(var)
      nam  = names(.data[k])

    var_name = data %>% dplyr::select(!!var) %>% unlist() %>% as.numeric()
    tot <- sum(!is.na(var_name))
    na1 <- sum(is.na(var_name))
    m1 <- mean(var_name, na.rm = T)
    m11 <- (sd(var_name, na.rm = T)/m1)*100
    if(any(plots == TRUE)){
      op <- par(mfrow = c(2,2),
                oma = c(0, 0, 2, 0) ,
                mar = c(2,2,2,2))
    boxplot(var_name, main = "With outliers")
    hist(var_name, main = "Without outliers", xlab=NA, ylab=NA)
    }
    outlier <- boxplot.stats(var_name, coef =  coef)$out
    dd = data %>% select(!!var)
    names_out = paste(which(dd[,1] %in% outlier), sep = " ")
    mo <- mean(outlier)
    var_name <- ifelse(var_name %in% outlier, NA, var_name)
    names = rownames(var_name)
    if(any(plots == TRUE)){
    boxplot(var_name, main = "With outliers")
    hist(var_name, main = "Without outliers", xlab=NA, ylab=NA)
    mtext(paste(nam), outer = TRUE, cex = 1.5)
    par(op)
    }
    na2 <- sum(is.na(var_name))
    if((na2 - na1)>0){
      if(verbose == TRUE){
        cat("\n----------------------------------------------------------------------------\n")
        cat("Level:", nam, "\n")
        cat("----------------------------------------------------------------------------\n")
      cat("Number of possible outliers:", na2 - na1, "\n")
      cat("Lines:", names_out, "\n")
      cat("Proportion: ", round((na2 - na1) / sum(!is.na(var_name))*100, 1), "%\n", sep = "")
      cat("Mean of the outliers:", round(mo, 3), "\n")
      m2 <- mean(var_name, na.rm = T)
      m22 <- (sd(var_name, na.rm = T)/m2)*100
      cat("With outliers:    mean = ", round(m1, 3), " | CV = ",round(m11, 3),"%", sep = "","\n")
      cat("Without outliers: mean = ", round(m2, 3), " | CV = ",round(m22, 3),"%", sep = "","\n")
      }
    }
    if((na2 - na1)==0){
      if(verbose == TRUE){
        cat("\n----------------------------------------------------------------------------\n")
        cat("Level:", nam, "\n")
        cat("----------------------------------------------------------------------------\n")
      cat("No outlier identified. \n")
      }
    }
    }
  } else {
    if(sum(lapply(.data, is.factor)==TRUE)>0){
    if(verbose == TRUE){
      message("The factors ", paste0(collapse = " ", names(.data[ , unlist(lapply(.data, is.factor)) ])),
              " were ignored. Use 'group_factors()' if you want to perform an analysis for each level of a factor.' ")
    }
    }
    var = dplyr::enquo(var)
    var_name = .data %>% dplyr::select(!!var) %>% unlist() %>% as.numeric()
    tot <- sum(!is.na(var_name))
    na1 <- sum(is.na(var_name))
    m1 <- mean(var_name, na.rm = T)
    m11 <- (sd(var_name, na.rm = T)/m1)*100
    if(plots == TRUE){
      op <- par(mfrow = c(2,2),
                oma = c(0, 0, 2, 0) ,
                mar = c(2,2,2,2))
      boxplot(var_name, main = "With outliers")
      hist(var_name, main = "With outliers", xlab=NA, ylab=NA)
    }
    outlier <- boxplot.stats(var_name, coef =  coef)$out
    dd = .data %>% select(!!var)
    names_out = paste(which(dd[,1] %in% outlier), sep = " ")
    mo <- mean(outlier)
    var_name <- ifelse(var_name %in% outlier, NA, var_name)
    names = rownames(var_name)
    if(any(plots == TRUE)){
      boxplot(var_name, main = "Without outliers")
      hist(var_name, main = "Without outliers", xlab=NA, ylab=NA)
      par(op)
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
}