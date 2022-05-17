#' Thennarasu's stability statistics
#' @description
#' `r badge('stable')`
#'
#' Performs a stability analysis based on Thennarasu (1995) statistics.
#'
#' @param .data The dataset containing the columns related to Environments,
#'   Genotypes, replication/block and response variable(s).
#' @param env The name of the column that contains the levels of the
#'   environments.
#' @param gen The name of the column that contains the levels of the genotypes.
#' @param resp The response variable(s). To analyze multiple variables in a
#'   single procedure use, for example, `resp = c(var1, var2, var3)`.
#' @param verbose Logical argument. If `verbose = FALSE` the code will run
#'   silently.
#' @return An object of class `Thennarasu`, which is a list containing the results
#'   for each variable used in the argument `resp`. For each variable, a
#'   tibble with the columns GEN, N1, N2, N3 and N4 is returned.
#' @md
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @references Thennarasu, K. 1995. On certain nonparametric procedures for
#'   studying genotype x environment interactions and yield stability. Ph.D.
#'   thesis. P.J. School, IARI, New Delhi, India.
#' @export
#' @examples
#' \donttest{
#' library(metan)
#' out <- Thennarasu(data_ge, ENV, GEN, GY)
#' print(out)
#' }
#'
Thennarasu <- function(.data, env, gen, resp, verbose = TRUE) {
  factors  <-
    .data %>%
    select({{env}}, {{gen}}) %>%
    mutate(across(everything(), as.factor))
  vars <-
    .data %>%
    select({{resp}}, -names(factors)) %>%
    select_numeric_cols()
  factors %<>% set_names("ENV", "GEN")
  listres <- list()
  nvar <- ncol(vars)
  if (verbose == TRUE) {
    pb <- progress(max = nvar, style = 4)
  }
  for (var in 1:nvar) {
    data <- factors %>%
      mutate(Y = vars[[var]])
    if(has_na(data)){
      data <- remove_rows_na(data)
      has_text_in_num(data)
    }
    data %<>% make_mat(GEN, ENV, Y)
    nr <- nrow(data)
    nc <- ncol(data)
    data_m <- as.matrix(data)
    N4.1 <- matrix(NA, nr, nc)
    N4 <- numeric()
    k <- 2
    data_r <- sweep(data_m, 1, rowMeans(data_m, na.rm = TRUE), FUN = "-")
    ranks <- apply(-data_r, 2, rank)
    ranks_y <- apply(-data_m, 2, rank)
    r_means <- rowMeans(ranks, na.rm = TRUE)
    r_means_y <- rowMeans(ranks_y, na.rm = TRUE)
    N1 <- round(rowMeans(abs(sweep(ranks, 1, apply(ranks, 1, median, na.rm = TRUE))), na.rm = TRUE),  4)
    N2 <- rowMeans(sweep(abs(sweep(ranks, 1, apply(ranks, 1, median, na.rm = TRUE))), 1, apply(ranks_y, 1, median, na.rm = TRUE), FUN = "/"))
    N3 <- round(sqrt(rowSums((sweep(ranks, 1, apply(ranks, 1, mean, na.rm = TRUE))^2) / nc)) / (rowMeans(ranks_y, na.rm = TRUE)), digits = 4)
    for (i in 1:nrow(data)) {
      for (j in 1:(nc - 1)) {
        N4.1[i, j] <- abs(ranks[i, j] - ranks[i, k])
        while (k < nc) k <- k + 1
      }
      N4[i] <- round((2/(nc * (nc - 1))) * (sum((N4.1[i, j])/(mean(ranks_y[i, ])))), digits = 4)
    }
    temp <- tibble(GEN = rownames(data),
                   Y = apply(data, 1, mean, na.rm = TRUE),
                   Y_R = rank(-Y),
                   N1 = N1,
                   N1_R = rank(N1),
                   N2 = N2,
                   N2_R = rank(N2),
                   N3 = N3,
                   N3_R = rank(N3),
                   N4 = N4,
                   N4_R = rank(N4))
    if (verbose == TRUE) {
      run_progress(pb,
                   actual = var,
                   text = paste("Evaluating trait", names(vars[var])))
    }
    listres[[paste(names(vars[var]))]] <- temp
  }
  return(structure(listres, class = "Thennarasu"))
}



#' Print an object ofclass `Thennarasu`
#'
#' Print the `Thennarasu` object in two ways. By default, the results are
#' shown in the R console. The results can also be exported to the directory
#' into a *.txt file.
#'
#'
#' @param x An object of class `Thennarasu`.
#' @param export A logical argument. If `TRUE`, a *.txt file is exported to
#'   the working directory.
#' @param file.name The name of the file if `export = TRUE`
#' @param digits The significant digits to be shown.
#' @param ... Options used by the tibble package to format the output. See
#'   [`tibble::print()`][tibble::formatting] for more details.
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @method print Thennarasu
#' @export
#' @examples
#' \donttest{
#' library(metan)
#' model <- Thennarasu(data_ge2, ENV, GEN, PH)
#' print(model)
#' }
print.Thennarasu <- function(x, export = FALSE, file.name = NULL, digits = 3, ...) {
  if (export == TRUE) {
    file.name <- ifelse(is.null(file.name) == TRUE, "Thennarasu summary", file.name)
    sink(paste0(file.name, ".txt"))
  }
  opar <- options(pillar.sigfig = digits)
  on.exit(options(opar))
  for (i in 1:length(x)) {
    var <- x[[i]]
    cat("Variable", names(x)[i], "\n")
    cat("---------------------------------------------------------------------------\n")
    cat("Thennarasu's stability indexes\n")
    cat("---------------------------------------------------------------------------\n")
    print(var)
    cat("---------------------------------------------------------------------------\n")
    cat("\n\n")
  }
  if (export == TRUE) {
    sink()
  }
}
