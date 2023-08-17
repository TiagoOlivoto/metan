#' Stability analysis based on Plaisted and Peterson (1959)
#' @description
#' `r badge('stable')`
#'
#' The function computes the stability as the arithmetic mean of the variance
#' component of the genotype-environment interaction between environment pairs
#' that includes a given genotype
#'
#'
#' @param .data The dataset containing the columns related to Environments,
#'   Genotypes, replication/block and response variable(s).
#' @param env The name of the column that contains the levels of the
#'   environments.
#' @param gen The name of the column that contains the levels of the genotypes.
#' @param rep The name of the column that contains the levels of the
#'   replications/blocks.
#' @param resp The response variable(s). To analyze multiple variables in a
#'   single procedure use, for example, `resp = c(var1, var2, var3)`.
#' @param verbose Logical argument. If `verbose = FALSE` the code will run
#'   silently.
#' @return An object of class `plaisted_peterson` containing the results for each
#'   variable used in the argument `resp`.
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @references Plaisted, R.L., and L.C. Peterson. 1959. A technique for
#'   evaluating the ability of selections to yield consistently in different
#'   locations or seasons. American Potato Journal 36(11): 381–385.
#'   \doi{10.1007/BF02852735}

#' @export
#' @examples
#' \donttest{
#' library(metan)
#'  plaisted_peterson(data_ge, ENV, GEN, REP, GY)
#'}
#'
plaisted_peterson <- function(.data, env, gen, rep, resp, verbose = TRUE) {
  factors  <-
    .data %>%
    select({{env}}, {{gen}}, {{rep}}) %>%
    mutate(across(everything(), as.factor))
  vars <-
    .data %>%
    select({{resp}}, -names(factors)) %>%
    select_numeric_cols()
  factors %<>% set_names(c("ENV", "GEN", "REP"))
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
    joint_an <- anova_joint(data, ENV, GEN, REP, Y, verbose = FALSE)
    r <- nlevels(data$REP)
    nenv <- nlevels(data$ENV)
    ngen <- nlevels(data$GEN)
    qmr <- joint_an$Y$anova[nrow(joint_an$Y$anova) - 3, 4]
    combs <- combn(ngen, 2)

    df2 <- make_mat(data, GEN, ENV, Y)
    gens <- rownames(df2)
    dists <- NULL
    for (j in 1:ncol(combs)) {
      gens_sel <- df2[c(combs[1, j], combs[2, j]), ]
      tmp <- dist(gens_sel) ^ 2
      dif <- diff(apply(gens_sel, 1, sum))
      tmp2 <- (r / 2) * (tmp - 1 / nenv * (dif ^ 2))
      dists[j] <- (tmp2 / (nenv - 1) - qmr) / r
    }
    mat <- matrix(NA, nrow = ngen, ncol = ngen)
    mat[lower.tri(mat)] <- dists
    mat <- make_sym(mat)
    mat_res <-
      as.data.frame(mat) |>
      mutate(theta = apply(mat, 1, mean, na.rm = TRUE),
             theta_perc = theta / sum(theta) * 100)
    rownames(mat_res) <- c(gens)
    colnames(mat_res) <- c(gens, "theta", "theta_prop")
    if (verbose == TRUE) {
      run_progress(pb,
                   actual = var,
                   text = paste("Evaluating trait", names(vars[var])))
    }
    listres[[paste(names(vars[var]))]] <- mat_res
  }
  return(structure(listres, class = "plaisted_peterson"))
}



#' Print an object of class plaisted_peterson
#'
#' Print the `plaisted_peterson` object in two ways. By default, the results
#' are shown in the R console. The results can also be exported to the directory
#' into a *.txt file.
#'
#'
#' @param x The `plaisted_peterson` x
#' @param export A logical argument. If `TRUE`, a *.txt file is exported to
#'   the working directory.
#' @param file.name The name of the file if `export = TRUE`
#' @param digits The significant digits to be shown.
#' @param ... Options used by the tibble package to format the output. See
#'   [`tibble::print()`][tibble::formatting] for more details.
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @method print plaisted_peterson
#' @export
#' @examples
#' \donttest{
#' library(metan)
#' eco <- ecovalence(data_ge2,
#'                   env = ENV,
#'                   gen = GEN,
#'                   rep = REP,
#'                   resp = PH)
#' print(eco)
#' }
print.plaisted_peterson <- function(x, export = FALSE, file.name = NULL, digits = 3, ...) {
  opar <- options(pillar.sigfig = digits)
  on.exit(options(opar))
  if (export == TRUE) {
    file.name <- ifelse(is.null(file.name) == TRUE, "plaisted_peterson print", file.name)
    sink(paste0(file.name, ".txt"))
  }
  for (i in 1:length(x)) {
    var <- x[[i]]
    cat("Variable", names(x)[i], "\n")
    print(var)
  }
  cat("\n\n\n")
  if (export == TRUE) {
    sink()
  }
}
