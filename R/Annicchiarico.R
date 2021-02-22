#' Annicchiarico's genotypic confidence index
#' @description
#' `r badge('stable')`
#'
#' Stability analysis using the known genotypic confidence index (Annicchiarico,
#' 1992).
#'
#' @param .data The dataset containing the columns related to Environments,
#'   Genotypes, replication/block and response variable(s)
#' @param env The name of the column that contains the levels of the
#'   environments.
#' @param gen The name of the column that contains the levels of the genotypes.
#' @param rep The name of the column that contains the levels of the
#'   replications/blocks
#' @param resp The response variable(s). To analyze multiple variables in a
#'   single procedure use, for example, `resp = c(var1, var2, var3)`.
#' @param prob The probability of error assumed.
#' @param verbose Logical argument. If `verbose = FALSE` the code will run
#'   silently.
#' @author Tiago Olivoto, \email{tiagoolivoto@@gmail.com}
#' @seealso [superiority()], [ecovalence()], [ge_stats()]
#' @references Annicchiarico, P. 1992. Cultivar adaptation and recommendation
#'   from alfalfa trials in Northern Italy. J. Genet. Breed. 46:269-278.
#' @return
#' A list where each element is the result for one variable and contains the
#' following data frames:
#' * **environments** Contains the mean, environmental index and
#' classification as favorable and unfavorable environments.
#' * **general** Contains the genotypic confidence index considering all
#' environments.
#' * **favorable** Contains the genotypic confidence index considering
#' favorable environments.
#' * **unfavorable** Contains the genotypic confidence index considering
#' unfavorable environments.
#' @md
#' @export
#' @examples
#' \donttest{
#' library(metan)
#'Ann <- Annicchiarico(data_ge2,
#'                     env = ENV,
#'                     gen = GEN,
#'                     rep = REP,
#'                     resp = PH)
#' print(Ann)
#'}
#'
Annicchiarico <- function(.data,
                          env,
                          gen,
                          rep,
                          resp,
                          prob = 0.25,
                          verbose = TRUE) {
  factors  <-
    .data %>%
    select({{env}}, {{gen}}, {{rep}}) %>%
    mutate(across(everything(), as.factor))
  vars <-
    .data %>%
    select({{resp}}, -names(factors)) %>%
    select_numeric_cols()
  factors %<>% set_names("ENV", "GEN", "REP")
  nvar <- ncol(vars)
  listres <- list()
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
    environments <-
      data %>%
      means_by(ENV, na.rm = TRUE) %>%
      add_cols(index = Y - mean(Y),
               class = ifelse(index < 0, "unfavorable", "favorable")) %>%
      as_tibble()
    data <- left_join(data, environments %>% select(ENV, class), by = "ENV")
    mat_g <- make_mat(data, row = GEN, col = ENV, value = Y)
    rp_g <- sweep(mat_g, 2, colMeans(mat_g, na.rm = TRUE), "/") * 100
    Wi_g <- rowMeans(rp_g, na.rm = TRUE) - qnorm(1 - prob) * apply(rp_g, 1, sd, na.rm = TRUE)
    general <- tibble(GEN = rownames(mat_g),
                      Y = rowMeans(mat_g, na.rm = TRUE),
                      Mean_rp = rowMeans(rp_g, na.rm = TRUE),
                      Sd_rp = apply(rp_g, 1, sd, na.rm = TRUE),
                      Wi = Wi_g,
                      rank = rank(-Wi_g))
    ge_mf <- subset(data, class == "favorable")
    mat_f <- dplyr::select_if(make_mat(ge_mf, row = GEN, col = ENV, value = Y), function(x) !any(is.na(x)))
    rp_f <- sweep(mat_f, 2, colMeans(mat_f, na.rm = TRUE), "/") * 100
    Wi_f <- rowMeans(rp_f, na.rm = TRUE) - qnorm(1 - prob) * apply(rp_f, 1, sd)
    favorable <- tibble(GEN = rownames(mat_f),
                        Y = rowMeans(mat_f, na.rm = TRUE),
                        Mean_rp = rowMeans(rp_f, na.rm = TRUE),
                        Sd_rp = apply(rp_f, 1, sd, na.rm = TRUE),
                        Wi = Wi_f,
                        rank = rank(-Wi_f))
    ge_mu <- subset(data, class == "unfavorable")
    mat_u <- dplyr::select_if(make_mat(ge_mu, row = GEN, col = ENV, value = Y), function(x) !any(is.na(x)))
    rp_u <- sweep(mat_u, 2, colMeans(mat_u, na.rm = TRUE), "/") * 100
    Wi_u <- rowMeans(rp_u, na.rm = TRUE) - qnorm(1 - prob) * apply(rp_u, 1, sd, na.rm = TRUE)
    unfavorable <- tibble(GEN = rownames(mat_u),
                          Y = rowMeans(mat_u, na.rm = TRUE),
                          Mean_rp = rowMeans(rp_u, na.rm = TRUE),
                          Sd_rp = apply(rp_u, 1, sd, na.rm = TRUE),
                          Wi = Wi_u,
                          rank = rank(-Wi_u))
    temp <- list(environments = environments,
                 general = general,
                 favorable = favorable,
                 unfavorable = unfavorable)
    if (verbose == TRUE) {
      run_progress(pb,
                   actual = var,
                   text = paste("Evaluating trait", names(vars[var])))
    }
    listres[[paste(names(vars[var]))]] <- temp
  }
  return(structure(listres, class = "Annicchiarico"))
}








#' Print an object of class Annicchiarico
#'
#' Print the `Annicchiarico` object in two ways. By default, the results
#' are shown in the R console. The results can also be exported to the directory
#' into a *.txt file.
#'
#'
#' @param x The `Annicchiarico` x
#' @param export A logical argument. If `TRUE`, a *.txt file is exported to
#'   the working directory.
#' @param file.name The name of the file if `export = TRUE`
#' @param digits The significant digits to be shown.
#' @param ... Options used by the tibble package to format the output. See
#'   [`tibble::print()`][tibble::formatting] for more details.
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @method print Annicchiarico
#' @export
#' @examples
#' \donttest{
#' library(metan)
#' Ann <- Annicchiarico(data_ge2,
#'   env = ENV,
#'   gen = GEN,
#'   rep = REP,
#'   resp = PH
#' )
#' print(Ann)
#' }
print.Annicchiarico <- function(x, export = FALSE, file.name = NULL, digits = 3, ...) {
  if (!class(x) == "Annicchiarico") {
    stop("The object must be of class 'Annicchiarico'")
  }
  opar <- options(pillar.sigfig = digits)
  on.exit(options(opar))
  if (export == TRUE) {
    file.name <- ifelse(is.null(file.name) == TRUE, "Annicchiarico print", file.name)
    sink(paste0(file.name, ".txt"))
  }
  for (i in 1:length(x)) {
    var <- x[[i]]
    cat("Variable", names(x)[i], "\n")
    cat("---------------------------------------------------------------------------\n")
    cat("Environmental index\n")
    cat("---------------------------------------------------------------------------\n")
    print(var$environments)
    cat("---------------------------------------------------------------------------\n")
    cat("Analysis for all environments\n")
    cat("---------------------------------------------------------------------------\n")
    print(var$general)
    cat("---------------------------------------------------------------------------\n")
    cat("Analysis for favorable environments\n")
    cat("---------------------------------------------------------------------------\n")
    print(var$favorable)
    cat("---------------------------------------------------------------------------\n")
    cat("Analysis for unfavorable environments\n")
    cat("---------------------------------------------------------------------------\n")
    print(var$unfavorable)
    cat("\n\n\n")
  }
  if (export == TRUE) {
    sink()
  }
}
