can_corr = function(.data = NULL, FG = NULL, SG = NULL, use = "cor",
  test = "Bartlett", prob = 0.05, center = TRUE, verbose = TRUE, collinearity = TRUE) {
  if (missing(.data) & missing(FG) || missing(SG)) {
    stop("No valid data imput for analysis.")
  }
  if (!missing(.data) & missing(FG) || missing(SG)) {
    stop("If a dataset is used as impute then 'FG' and 'SG' must be declared.")
  }
  if (!use %in% c("cov", "cor")) {
    stop("The argument  'use' is incorrect, it should be 'cov' or 'cor'.")
  }
  if (!missing(.data)) {
    FG = as.data.frame(dplyr::select(.data, !!!dplyr::quos(!!dplyr::enquo(FG))))
    SG = as.data.frame(dplyr::select(.data, !!!dplyr::quos(!!dplyr::enquo(SG))))
  }

  if (missing(.data) & !missing(FG) & !missing(SG)) {
    if (!is.data.frame(FG)) {
      stop("'FG' should be data frame.")
    }
    if (!is.data.frame(SG)) {
      stop("'SG' should be data frame.")
    }
  }

  if (nrow(FG) != nrow(SG)) {
    stop("The number of observations of 'FG', should be equal to 'SG'.")
  }
  if (ncol(FG) > ncol(SG)) {
    stop("The number of variables in 'FG' should be lesser than or equal to the number of variables in 'SG'.")
  }
  if (!test %in% c("Bartlett", "Rao")) {
    stop("The argument 'test' is incorrect, it should be 'Bartlett' or 'Rao'.")
  }
  if (!is.numeric(prob) | prob <= 0 || prob > 1) {
    stop("The argument 'prob' is incorrect. It should be numeric with values between 0 and 1.")
  }

  if (use == "cov") {
    MC = cov(cbind(FG, SG))
    S11 = cov(FG)
    S22 = cov(SG)
    S12 = cov(FG, SG)
    S21 = cov(SG, FG)
  }
  if (use == "cor") {
    MC = cor(cbind(FG, SG))
    S11 = cor(FG)
    S22 = cor(SG)
    S12 = cor(FG, SG)
    S21 = cor(SG, FG)
  }
  M1 <- eigen(S11)
  megval1 <- M1$values
  megvec1 <- M1$vectors
  S11_12 = megvec1 %*% diag(1/sqrt(megval1)) %*% t(megvec1)
  S22_Inv = solve(S22)
  M2 <- eigen(S11_12 %*% S12 %*% S22_Inv %*% S21 %*% S11_12)
  megval2 <- M2$values
  megvec2 <- M2$vectors
  mtr = megval2
  varuv <- as.data.frame(matrix(NA, length(mtr), 3))
  rownames(varuv) <- paste("U", 1:length(mtr), "V", 1:length(mtr),
    sep = "")
  colnames(varuv) <- c("Variance", "Proportion", "Cum_proportion")
  varuv[, "Variance"] <- mtr
  varuv[, "Proportion"] <- (mtr/sum(mtr)) * 100
  varuv[, "Cum_proportion"] <- cumsum(varuv[, "Proportion"])
  coruv <- as.matrix(sqrt(mtr), ncol = length(coruv), nrow = 1)
  rownames(coruv) <- paste("U", 1:length(coruv), "V", 1:length(coruv),
    sep = "")
  colnames(coruv) <- c("Correlation")
  Coef_FG = S11_12 %*% megvec2
  rownames(Coef_FG) <- colnames(FG)
  colnames(Coef_FG) <- paste("U", 1:ncol(Coef_FG), sep = "")
  Coef_SG = S22_Inv %*% S21 %*% Coef_FG %*% solve(diag(sqrt(megval2)))
  colnames(Coef_SG) <- paste("V", 1:ncol(Coef_SG), sep = "")
  M3 <- eigen(diag(diag(S11)))
  megval3 <- M3$values
  megvec3 <- M3$vectors
  D11_12 = megvec3 %*% diag(1/sqrt(megval3)) %*% t(megvec3)
  M4 <- eigen(diag(diag(S22)))
  megval4 <- M4$values
  megvec4 <- M4$vectors
  D22_12 = megvec4 %*% diag(1/sqrt(megval4)) %*% t(megvec4)
  Rux <- t(t(Coef_FG) %*% S11 %*% D11_12)
  rownames(Rux) <- colnames(FG)
  Rvy <- t(t(Coef_SG) %*% S22 %*% D22_12)
  rownames(Rvy) <- colnames(SG)

  if (center == TRUE) {
    FG_A = scale(FG, center = TRUE, scale = FALSE)
    SG_A = scale(SG, center = TRUE, scale = FALSE)
  } else {
    FG_A = FG
    SG_A = SG
  }
  FG_A[is.na(FG_A)] = 0
  SG_A[is.na(SG_A)] = 0
  FG_SC = FG_A %*% Coef_FG
  SG_SC = SG_A %*% Coef_SG
  if (test == "Bartlett") {
    n <- nrow(FG)
    p <- ncol(FG)
    q <- ncol(SG)
    QtdF <- length(coruv)
    Bartlett <- as.data.frame(matrix(NA, QtdF, 5))
    colnames(Bartlett) <- c("Canonical_pairs", "Lambda_Wilks",
      "Chi_square", "DF", "p_value")
    Bartlett[, 1] <- paste("U", 1:QtdF, "V", 1:QtdF, sep = "")
    i = 1
    for (i in 1:QtdF) {
      Lambda <- prod(1 - coruv[i:QtdF]^2)
      chisq <- -((n - 1) - (p + q + 1)/2) * log(Lambda)
      gl <- (p - i + 1) * (q - i + 1)
      pValor <- pchisq(chisq, gl, ncp = 0, lower.tail = F)
      Bartlett[i, 2] <- round(Lambda, 5)
      Bartlett[i, 3] <- round(chisq, 5)
      Bartlett[i, 4] <- gl
      Bartlett[i, 5] <- round(pValor, 5)
    }
    teste <- Bartlett
  }
  if (test == "Rao") {
    n <- nrow(FG)
    p1 <- ncol(FG)
    q1 <- ncol(SG)
    QtdF <- length(coruv)
    Rao <- as.data.frame(matrix(NA, QtdF, 6))
    colnames(Rao) <- c("Canonical pairs", "Lambda_Wilks",
      "F_value", "DF1", "DF2", "p_value")
    Rao[, 1] <- paste("U", 1:QtdF, "V", 1:QtdF, sep = "")
    for (i in 1:QtdF) {
      p <- p1 - i + 1
      q <- q1 - i + 1
      t <- (n - 1) - (p + q + 1)/2
      s <- ifelse((p^2 + q^2) <= 5, 1, sqrt((p^2 * q^2 -
        4)/(p^2 + q^2 - 5)))
      Lambda <- prod(1 - coruv[i:QtdF]^2)
      gl1 <- p * q
      gl2 <- (1 + t * s - p * q/2)
      FVAL <- ((1 - Lambda^(1/s))/Lambda^(1/s)) *
        gl2/gl1
      pValor <- pf(FVAL, gl1, gl2, ncp = 0, lower.tail = FALSE)
      Rao[i, 2] <- round(Lambda, 5)
      Rao[i, 3] <- round(FVAL, 5)
      Rao[i, 4] <- gl1
      Rao[i, 5] <- round(gl2, 5)
      Rao[i, 6] <- round(pValor, 5)
    }
    teste <- Rao
  }
  results = data.frame(cbind(cbind(varuv, coruv), teste[-1]))
if (collinearity == TRUE){
  colin = list(FG = colindiag(FG, verbose = FALSE),
               SG = colindiag(SG, verbose = FALSE))
} else {
  colin = NULL
}
if (verbose == TRUE){
  cat("---------------------------------------------------------------------------\n")
  cat("Matrix (correlation/covariance) between variables of first group (FG)\n")
  cat("---------------------------------------------------------------------------\n")
  print(S11)
    if (collinearity == TRUE){
  cat("\n---------------------------------------------------------------------------\n")
  cat("Collinearity within first group \n")
  cat("---------------------------------------------------------------------------\n")
  colindiag(FG)
    }
  cat("\n---------------------------------------------------------------------------\n")
  cat("Matrix (correlation/covariance) between variables of second group (SG)\n")
  cat("---------------------------------------------------------------------------\n")
  print(S22)
  if (collinearity == TRUE){
  cat("\n---------------------------------------------------------------------------\n")
  cat("Collinearity within second group \n")
  cat("---------------------------------------------------------------------------\n")
  colindiag(SG)
  }
  cat("\n---------------------------------------------------------------------------\n")
  cat("Matrix (correlation/covariance) between FG and SG)\n")
  cat("---------------------------------------------------------------------------\n")
  print(S12)
  cat("\n---------------------------------------------------------------------------\n")
  cat("Correlation of the canonical pairs and hipotesis testing \n")
  cat("---------------------------------------------------------------------------\n")
  print(results)
  cat("\n---------------------------------------------------------------------------\n")
  cat("Canonical coefficients of the first group \n")
  cat("---------------------------------------------------------------------------\n")
  print(Coef_FG)
  cat("\n---------------------------------------------------------------------------\n")
  cat("Canonical coefficients of the second group \n")
  cat("---------------------------------------------------------------------------\n")
  print(Coef_SG)
  cat("\n---------------------------------------------------------------------------\n")
  cat("Canonical loads of the first group \n")
  cat("---------------------------------------------------------------------------\n")
  print(Rux)
  cat("\n---------------------------------------------------------------------------\n")
  cat("Canonical loads of the second group \n")
  cat("---------------------------------------------------------------------------\n")
  print(Rvy)
}
  invisible(structure(list(Matrix = MC, MFG = S11, MSG = S22,
    MFG_SG = S12, Coef_FG = Coef_FG, Coef_SG = Coef_SG,
    Loads_FG = Rux, Loads_SG = Rvy, Score_FG = FG_SC, Score_SG = SG_SC,
    Sigtest = results, collinearity = colin), class = "can_cor"))
}
