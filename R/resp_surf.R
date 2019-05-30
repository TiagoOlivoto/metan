#' Response surface model
#'
#' Compute a surface model and find the best combination of factor1 and factor2
#' to obtain the stacionary point.
#'
#'
#' @param .data The dataset containing the columns related to Environments,
#' factor1, factor2, replication/block and response variable(s).
#' @param factor1 The first factor, for example, dose of Nitrogen.
#' @param factor2 The second factor, for example, dose of potassium.
#' @param rep The name of the column that contains the levels of the
#' replications/blocks.
#' @param resp The response variable(s).
#' @param prob The probability error.
#' @param verbose If \code{verbose = TRUE} then some results are shown in the
#' console.
#' @export
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#'
resp_surf = function(.data,
                     factor1,
                     factor2,
                     rep,
                     resp,
                     prob = 0.05,
                     verbose = TRUE) {

data = .data  %>%
        dplyr::select(!!dplyr::enquo(factor1),
                      !!dplyr::enquo(factor2),
                      !!dplyr::enquo(rep),
                      !!dplyr::enquo(resp)) %>%
  as.data.frame()
  names = colnames(data)
  A = names[1]
  D = names[2]
  Bloco = names[3]
  Resp = names[4]
  F1 = as.formula(paste0(Resp, "~", paste(Bloco),"+", paste(A), "+", paste(D), "+", paste(A), "*", paste(D)))
  ANOVA = aov(F1, data = data)
  sum_test = unlist(summary(ANOVA))
  pValue_Bloco = sum_test["Pr(>F)1"]
  pValue_A = sum_test["Pr(>F)2"]
  pValue_D = sum_test["Pr(>F)3"]
  pValue_AD = sum_test["Pr(>F)4"]
  F2 = as.formula(paste0(Resp, " ~", paste(Bloco),"+", paste(A),"*", paste(D), "+", "I(",paste(A),"^2)", "+","I(",paste(D), "^2)"))
  SurfMod = lm(F2, data=data)
  B0 = SurfMod$coef[2]
  B1 = SurfMod$coef[3]
  B2 = SurfMod$coef[4]
  B11 = SurfMod$coef[5]
  B22 = SurfMod$coef[6]
  B12 = SurfMod$coef[7]
  B12c = B12/2
  P = cbind(c(B11,B12c), c(B12c,B22))
  P = as.matrix(P)
  invA = solve(P)
  X = cbind(c(B1, B2))
  Pontos = -0.5*(invA %*% X)
  dA = Pontos[1]
  dD = Pontos[2]
  AV1 = eigen(P)$values[1]
  AV2 = eigen(P)$values[2]
  results = dplyr::mutate(data,
                          predicted = SurfMod$fitted.values,
                          residuals = SurfMod$residuals)

if (verbose == TRUE) {
  cat("-----------------------------------------------------------------\n")
  cat("Result for the analysis of variance", "\n")
  cat("Model: Y = m + bk + Ai + Dj + (AD)ij + eijk", "\n")
  cat("-----------------------------------------------------------------\n")
  cat("\n")
  print(summary(ANOVA))
  cat("\n")
  Norm = shapiro.test(ANOVA$residuals)
  cat("-----------------------------------------------------------------\n")
  cat("Shapiro-Wilk's test for normality of residuals:", "\n")
  cat("W = ",Norm$statistic, "p-valor = ", Norm$p.value, "\n")
  if (Norm$p.value>0.05){
    cat("According to S-W test, residuals may be considered normal.","\n")
  }
  cat("-----------------------------------------------------------------\n")
  cat("Anova table for the response surface model", "\n")
  print(F2)
  cat("-----------------------------------------------------------------\n")

  print(anova(SurfMod))
  cat("-----------------------------------------------------------------\n")

  cat("Parameter estimates for response surface model", "\n")
  print(summary(SurfMod))
  cat("-----------------------------------------------------------------\n")

  cat("Model equation for response surface model", "\n")
  cat("Y = B0 + B1*A + B2*D + B11*A^2 + B22*D^2 + B12*A*D", "\n")
  cat("-----------------------------------------------------------------\n")
  cat("Estimated parameters", "\n")
  cat(paste0("B0: ", format(round(B0,7), nsmall = 7),"\n"))
  cat(paste0("B1: ", format(round(B1,7), nsmall = 7),"\n"))
  cat(paste0("B2: ", format(round(B2,7), nsmall = 7),"\n"))
  cat(paste0("B11: ", format(round(B11,7), nsmall = 7),"\n"))
  cat(paste0("B22: ", format(round(B22,7), nsmall = 7),"\n"))
  cat(paste0("B12: ", format(round(B12,7), nsmall = 7),"\n"))
  cat("-----------------------------------------------------------------\n")

  cat("Matrix of parameters (A)","\n")
  cat(paste0(format(round(P[1,1],7), nsmall = 7)),"  ", format(round(P[1,2],7), nsmall = 7),"\n")
  cat(paste0(format(round(P[2,1],7), nsmall = 7)),"  ", format(round(P[2,2],7), nsmall = 7),"\n")
  cat("-----------------------------------------------------------------\n")

  cat("Inverse of the matrix A (invA)","\n")
  cat(paste0(format(round(invA[1,1],7), nsmall = 7)),"  ", format(round(invA[1,2],7), nsmall = 7),"\n")
  cat(paste0(format(round(invA[2,1],7), nsmall = 7)),"  ", format(round(invA[2,2],7), nsmall = 7),"\n")
  cat("-----------------------------------------------------------------\n")

  cat("Vetor of parameters B1 e B2 (X)","\n")
  cat(paste0("B1: ", format(round(B1,7), nsmall = 7),"\n"))
  cat(paste0("B2: ", format(round(B2,7), nsmall = 7),"\n"))
  cat("-----------------------------------------------------------------\n")

  cat("Equation for the optimal points (A and D)","\n")
  cat("-0.5*(invA*X)")
  cat(paste0("\nEigenvalue 1: " ,round(AV1,6),
             "\nEigenvalue 2: " ,round(AV2,6)))
  cat("\n")
  if (AV1 > 0 && AV2 > 0) {
    cat(paste0("Stacionary point is minimum!"))
  } else if (AV1 < 0 && AV2 < 0) {
    cat(paste0("Stacionary point is maximum!"))
  } else
    cat(paste0("The stacionary point is outside the intervals of the trataments"))
  cat("\n")
  cat("-----------------------------------------------------------------\n")
  cat("Stacionary point obtained with the following original units:", "\n")
  cat(paste0("Optimal dose (",A,"): ", round(dA,4),"\n"))
  cat(paste0("Optimal dose (",D,"): ", round(dD,4),"\n"))
  cat("-----------------------------------------------------------------\n")
  cat("Fitted model", "\n")
  cat(paste0("A = ", A,"\n"))
  cat(paste0("D = ", D,"\n"))
  cat(paste0("y = " , round(B0,5), "+",round(B1,5),"A+",round(B2,5),"D+" , round(B11,5),"A^2+" , round(B22,5),"D^2+" , round(B12,5),"A*D","\n"))
  cat("-----------------------------------------------------------------\n")

  cat("Observed, predicted and residuals","\n")
  cat("-----------------------------------------------------------------\n")
  print(results)
  cat("-----------------------------------------------------------------\n")

  pvalor.shapiro = shapiro.test(results$residuals)$p.value
  cat("Shapiro-Wilk normality test\n")
  cat("p-value: ", pvalor.shapiro, "\n")
  if (pvalor.shapiro < 0.05) {
    cat("WARNING: at 5% of significance, residuals can not be considered normal!", "\n")
    cat("------------------------------------------------------------------")
  } else {
    cat("According to Shapiro-Wilk normality test at 5% of significance, residuals can be considered normal.", "\n")
    cat("------------------------------------------------------------------\n")
  }
}

  return(structure(list(results = results,
                        model = SurfMod),
                   class = "resp_sup"))
}
