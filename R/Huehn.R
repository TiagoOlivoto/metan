#' Huehn's stability statistics
#' @description
#' `r badge('stable')`
#'
#' Performs a stability analysis based on Huehn (1979) statistics. The four
#' nonparametric measures of phenotypic stability are: S1 (mean of the absolute
#' rank differences of a genotype over the n environments), S2 (variance among
#' the ranks over the k environments), S3 (sum of the absolute deviations), and
#' S6 (relative sum of squares of rank for each genotype).
#' @param .data The dataset containing the columns related to Environments,
#'   Genotypes, replication/block and response variable(s).
#' @param env The name of the column that contains the levels of the
#'   environments.
#' @param gen The name of the column that contains the levels of the genotypes.
#' @param resp The response variable(s). To analyze multiple variables in a
#'   single procedure use, for example, `resp = c(var1, var2, var3)`.
#' @param verbose Logical argument. If `verbose = FALSE` the code will run
#'   silently.
#' @return An object of class `Huehn`, which is a list containing the results
#'   for each variable used in the argument `resp`. For each variable, a
#'   tibble with the following columns is returned.
#' * **GEN** The genotype's code.
#' * **Y** The mean for the response variable.
#' * **S1** Mean of the absolute rank differences of a genotype over the n environments.
#' * **S2** variance among the ranks over the k environments.
#' * **S3** Sum of the absolute deviations.
#' * **S6** Relative sum of squares of rank for each genotype.
#' @md
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @references Huehn, V.M. 1979. Beitrage zur erfassung der phanotypischen
#'   stabilitat. EDV Med. Biol. 10:112.
#'
#' @export
#' @examples
#' \donttest{
#' library(metan)
#' out <- Huehn(data_ge2, ENV, GEN, PH)
#' print(out)
#' }
#'
Huehn <- function(.data, env, gen, resp, verbose = TRUE) {
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
    data <- make_mat(data, GEN, ENV, Y)
    nr <- nrow(data)
    nc <- ncol(data)
    data_m <- as.matrix(data)
    S.1 <- matrix(NA, nr, nc)
    S1 <- numeric()
    k <- 2
    data_r <- sweep(data_m, 1, rowMeans(data_m, na.rm = TRUE), FUN = "-") + mean(data_m, na.rm = TRUE)
    ranks <- apply(data_r, 2, rank)
    ranks_y <- apply(data_m, 2, rank)
    r_means <- rowMeans(ranks, na.rm = TRUE)
    r_means_y <- rowMeans(ranks_y, na.rm = TRUE)
    S2 <- round(rowSums(sweep(ranks, 1, r_means)^2, na.rm = TRUE) / (nc-1), 4)
    S3 <- round(rowSums(sweep(ranks_y, 1, r_means_y)^2, na.rm = TRUE) / r_means_y, 4)
    S6 <- round(rowSums(abs(sweep(ranks_y, 1, r_means_y)), na.rm = TRUE) / r_means_y, 4)
    for (i in 1:nrow(data)) {
      for (j in 1:(nc - 1)) {
        S.1[i, j] <- abs(ranks[i, j] - ranks[i, k])
        while (k < nc) k <- k + 1
      }
      S1[i] <- round((2 * (sum(S.1[i, j])))/(nc * (nc - 1)), digits = 4)
    }
    temp <- tibble(GEN = rownames(data),
                   Y = apply(data, 1, mean, na.rm = TRUE),
                   Y_R = rank(-Y),
                   S1 = S1,
                   S1_R = rank(S1),
                   S2 = S2,
                   S2_R = rank(S2),
                   S3 = S3,
                   S3_R = rank(S3),
                   S6 = S6,
                   S6_R = rank(S6))
    if (verbose == TRUE) {
      run_progress(pb,
                   actual = var,
                   text = paste("Evaluating trait", names(vars[var])))
    }
    listres[[paste(names(vars[var]))]] <- temp
  }
  return(structure(listres, class = "Huehn"))
}





#' Print an object ofclass `Huehn`
#'
#' Print the `Huehn` object in two ways. By default, the results are
#' shown in the R console. The results can also be exported to the directory
#' into a *.txt file.
#'
#'
#' @param x An object of class `Huehn`.
#' @param export A logical argument. If `TRUE`, a *.txt file is exported to
#'   the working directory.
#' @param file.name The name of the file if `export = TRUE`
#' @param digits The significant digits to be shown.
#' @param ... Options used by the tibble package to format the output. See
#'   [`tibble::print()`][tibble::formatting] for more details.
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @method print Huehn
#' @export
#' @examples
#' \donttest{
#' library(metan)
#' model <- Huehn(data_ge2, ENV, GEN, PH)
#' print(model)
#' }
print.Huehn <- function(x, export = FALSE, file.name = NULL, digits = 3, ...) {
  if (!class(x) == "Huehn") {
    stop("The object must be of class 'Huehn'")
  }
  if (export == TRUE) {
    file.name <- ifelse(is.null(file.name) == TRUE, "Huehn summary", file.name)
    sink(paste0(file.name, ".txt"))
  }
  opar <- options(pillar.sigfig = digits)
  on.exit(options(opar))
  for (i in 1:length(x)) {
    var <- x[[i]]
    cat("Variable", names(x)[i], "\n")
    cat("---------------------------------------------------------------------------\n")
    cat("Huehn's stability indexes\n")
    cat("---------------------------------------------------------------------------\n")
    print(var)
    cat("---------------------------------------------------------------------------\n")
    cat("\n\n")
  }
  if (export == TRUE) {
    sink()
  }
}
