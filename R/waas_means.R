#' Weighted Average of Absolute Scores
#' @description
#' `r badge('stable')`
#'
#' Compute the Weighted Average of Absolute Scores (Olivoto et al., 2019) based
#' on means for genotype-environment data as follows:
#'\loadmathjax
#'\mjsdeqn{WAAS_i = \sum_{k = 1}^{p} |IPCA_{ik} \times EP_k|/ \sum_{k =
#' 1}^{p}EP_k}
#'
#' where \mjseqn{WAAS_i} is the weighted average of absolute scores of the
#' *i*th genotype; \mjseqn{PCA_{ik}} is the score of the *i*th genotype
#' in the *k*th IPCA; and \mjseqn{EP_k} is the explained variance of the *k*th
#' IPCA for *k = 1,2,..,p*, where *p* is the number of IPCAs that
#' explain at least an amount of the genotype-interaction variance declared in
#' the argument `min_expl_var`.
#'
#' @param .data The dataset containing the columns related to Environments,
#'   Genotypes, replication/block and response variable(s).
#' @param env The name of the column that contains the levels of the
#'   environments.
#' @param gen The name of the column that contains the levels of the genotypes.
#' @param resp The response variable(s). To analyze multiple variables in a
#'   single procedure a vector of variables may be used. For example `resp
#'   = c(var1, var2, var3)`. Select helpers are also allowed.
#' @param mresp  The new maximum value after rescaling the response variable. By
#'   default, all variables in `resp` are rescaled so that de maximum value
#'   is 100 and the minimum value is 0 (i.e., `mresp = NULL`). It must be a
#'   character vector of the same length of `resp` if rescaling is assumed
#'   to be different across variables, e.g., if for the first variable smaller
#'   values are better and for the second one, higher values are better, then
#'   `mresp = c("l, h")` must be used. Character value of length 1 will be
#'   recycled with a warning message.
#' @param wresp The weight for the response variable(s) for computing the WAASBY
#'   index. Must be a numeric vector of the same length of `resp`. Defaults
#'   to 50, i.e., equal weights for stability and mean performance.
#' @param min_expl_var The minimum explained variance. Defaults to 85.
#'   Interaction Principal Compoment Axis are iteractively retained up to the
#'   explained variance (eigenvalues in the singular value decomposition of the
#'   matrix with the interaction effects) be greather than or equal to
#'   `min_expl_var`. For example, if the explained variance (in percentage)
#'   in seven possible IPCAs are `56, 21, 9, 6, 4, 3, 1`, resulting in a
#'   cumulative proportion of `56,  77,  86,  92,  96, 99, 100`, then
#'   `p = 3`, i.e., three IPCAs will be used to compute the index WAAS.
#' @param verbose Logical argument. If `verbose = FALSE` the code is run
#'   silently.
#' @param ... Arguments passed to the function
#'   [impute_missing_val()] for imputation of missing values in case
#'   of unbalanced data.
#' @references
#' Olivoto, T., A.D.C. L{\'{u}}cio, J.A.G. da silva, V.S. Marchioro, V.Q. de
#' Souza, and E. Jost. 2019a. Mean performance and stability in
#' multi-environment trials I: Combining features of AMMI and BLUP techniques.
#' Agron. J. 111:2949-2960. \doi{10.2134/agronj2019.03.0220}
#'
#' @return An object of class `waas_means` with the following items for each
#'   variable:
#'
#' * **model** A data frame with the response variable, the scores of all
#' Principal Components, the estimates of Weighted Average of Absolute Scores,
#' and WAASY (the index that consider the weights for stability and productivity
#' in the genotype ranking.
#' * **ge_means** A tbl_df containing the genotype-environment means.
#' * **ge_eff** A *gxe* matrix containing the genotype-environment effects.
#' * **eigenvalues** The eigenvalues from the singular value decomposition
#' of the matrix withe the genotype-environment interaction effects.
#' * **proportion** The proportion of the variance explained by each IPCA.
#' * **cum_proportion** The cumulative proportion of the variance explained.
#' @md
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @seealso [waas()] [waasb()]
#' @export
#' @examples
#'\donttest{
#' library(metan)
#' # Data with replicates
#' model <- waas(data_ge,
#'               env = ENV,
#'               gen = GEN,
#'               rep = REP,
#'               resp = everything())
#'
#' # Based on means of genotype-environment data
#' data_means <- means_by(data_ge, ENV, GEN)
#' model2 <- waas_means(data_ge,
#'                      env = ENV,
#'                      gen = GEN,
#'                      resp = everything())
#' # The index WAAS
#' get_model_data(model, what = "OrWAAS")
#' get_model_data(model2, what = "OrWAAS")
#'
#'}
#'
waas_means <- function(.data,
                       env,
                       gen,
                       resp,
                       mresp = NULL,
                       wresp = NULL,
                       min_expl_var = 85,
                       verbose = TRUE,
                       ...){
  if(is.numeric(mresp)){
    stop("Using a numeric vector in 'mresp' is deprecated as of metan 1.9.0. use 'h' or 'l' instead.\nOld code: 'mresp = c(100, 100, 0)'.\nNew code: 'mresp = c(\"h, h, l\")", call. = FALSE)
  }
  factors  <-
    .data %>%
    select({{env}}, {{gen}}) %>%
    mutate(across(everything(), as.factor))
  vars <-
    .data %>%
    select({{resp}}, -names(factors)) %>%
    select_numeric_cols()
  factors %<>% set_names("ENV", "GEN")
  nvar <- ncol(vars)
  if (verbose == TRUE) {
    pb <- progress(max = nvar, style = 4)
  }
  if (is.null(mresp)) {
    mresp <- replicate(nvar, 100)
    minresp <- 100 - mresp
  } else {
    mresp <- unlist(strsplit(mresp, split="\\s*(\\s|,)\\s*")) %>% all_lower_case()
    if(!any(mresp %in% c("h", "l", "H", "L"))){
      if(!mresp[[1]] %in% c("h", "l")){
        stop("Argument 'mresp' must have only h or l.", call. = FALSE)
      } else{
        warning("Argument 'mresp' must have only h or l. Setting mresp = ", mresp[[1]],
                " to all the ", nvar, " variables.", call. = FALSE)
        mresp <- replicate(nvar, mresp[[1]])
      }
    }
    if (length(mresp) != nvar) {
      warning("Invalid length in 'mresp'. Setting mresp = ", mresp[[1]],
              " to all the ", nvar, " variables.", call. = FALSE)
      mresp <- replicate(nvar, mresp[[1]])
    }

    mresp <- ifelse(mresp == "h", 100, 0)
    minresp <- 100 - mresp
  }
  if (is.null(wresp)) {
    PesoResp <- replicate(nvar, 50)
    PesoWAASB <- 100 - PesoResp
  } else {
    if (length(wresp) != nvar) {
      stop("The length of the numeric vector 'wresp' must be ", nvar, ", the number of variables in argument 'resp'")
    }
    if (min(wresp) < 0 | max(wresp) > 100) {
      stop("The range of the numeric vector 'wresp' must be equal between 0 and 100.")
    }
    PesoResp <- wresp
    PesoWAASB <- 100 - PesoResp
  }
  listres <- list()
  vin <- 0
  for (var in 1:nvar) {
    data <- factors %>%
      mutate(Y = vars[[var]])
    check_labels(data)
    if(has_na(data)){
      data <- remove_rows_na(data)
      has_text_in_num(data)
    }
    data <-
      means_by(data, GEN, ENV) %>%
      make_mat(GEN, ENV, Y)
    if(has_na(data)){
      data <- impute_missing_val(data, verbose = verbose, ...)$.data
      warning("Data imputation used to fill the GxE matrix", call. = FALSE)
    }
    data %<>% make_long()

    vin <- vin + 1
    MGEN <-
      means_by(data, GEN) %>%
      rename(Code = GEN) %>%
      add_cols(type = "GEN")%>%
      select_cols(type, Code, Y)
    MENV <-
      means_by(data, ENV) %>%
      rename(Code = ENV) %>%
      add_cols(type = "ENV") %>%
      select_cols(type, Code, Y)
    ngen <- nrow(MGEN)
    nenv <- nrow(MENV)
    minimo <- min(ngen - 1, nenv - 1)
    res <- residuals(lm(Y ~ GEN + ENV, data = data))
    res_mat <- t(matrix(res, nenv, byrow = T))
    colnames(res_mat) <- get_levels(MENV, Code)
    rownames(res_mat) <- get_levels(MGEN, Code)
    s <- svd(res_mat)
    U <- s$u[, 1:minimo]
    LL <- diag(s$d[1:minimo])
    V <- s$v[, 1:minimo]
    SS <- (s$d[1:minimo]^2)
    SUMA <- sum(SS)
    weights <- (1/SUMA) * SS * 100
    naxis <- 1
    while (cumsum(weights)[naxis] < min_expl_var)  {
      naxis <- naxis + 1
    }
    SCOREG <- U %*% LL^0.5 %>%
      as.data.frame() %>%
      set_names(paste("PC", 1:minimo, sep = ""))
    SCOREE <- V %*% LL^0.5 %>%
      as.data.frame() %>%
      set_names(paste("PC", 1:minimo, sep = ""))
    Escores <-
      rbind(cbind(MGEN, SCOREG),
            cbind(MENV, SCOREE))
    WAAS <-
      Escores %>%
      select(contains("PC")) %>%
      abs() %>%
      t() %>%
      as.data.frame() %>%
      slice(1:naxis) %>%
      mutate(Percent = weights[1:naxis])
    WAASAbs <- mutate(Escores, WAAS = sapply(WAAS[, -ncol(WAAS)], weighted.mean, w = WAAS$Percent))
    if (nvar > 1) {
      WAASAbs %<>% group_by(type) %>%
        mutate(PctResp = (mresp[vin] - minresp[vin])/(max(Y) - min(Y)) * (Y - max(Y)) +
                 mresp[vin], PctWAAS = (0 - 100)/(max(WAAS) -  min(WAAS)) * (WAAS - max(WAAS)) + 0,
               wRes = PesoResp[vin],
               wWAAS = PesoWAASB[vin],
               OrResp = rank(-Y),
               OrWAAS = rank(WAAS),
               OrPC1 = rank(abs(PC1)),
               WAASY = ((PctResp * wRes) + (PctWAAS * wWAAS))/(wRes +  wWAAS),
               OrWAASY = rank(-WAASY)) %>%
        ungroup()
    }
    else {
      WAASAbs %<>% group_by(type) %>%
        mutate(PctResp = (mresp - minresp)/(max(Y) - min(Y)) * (Y - max(Y)) +
                 mresp, PctWAAS = (0 - 100)/(max(WAAS) - min(WAAS)) * (WAAS - max(WAAS)) + 0,
               wRes = PesoResp,
               wWAAS = PesoWAASB,
               OrResp = rank(-Y),
               OrWAAS = rank(WAAS),
               OrPC1 = rank(abs(PC1)),
               WAASY = ((PctResp * wRes) + (PctWAAS * wWAAS))/(wRes + wWAAS),
               OrWAASY = rank(-WAASY)) %>%
        ungroup()
    }
    temp <- structure(list(model = WAASAbs,
                           ge_means = data,
                           ge_eff = res_mat,
                           eigenvalues = SS,
                           proportion = weights,
                           cum_proportion = cumsum(weights)),
                      class = "waas_means")
    if (verbose == TRUE) {
      run_progress(pb,
                   actual = var,
                   text = paste("Evaluating trait", names(vars[var])))
    }
    listres[[paste(names(vars[var]))]] <- temp
  }
  return(structure(listres, class = "waas_means"))
}


#' Print an object of class waas_means
#'
#' Print the `waas_means` object in two ways. By default, the results are shown
#' in the R console. The results can also be exported to the directory.
#'
#'
#' @param x An object of class `waas_means`.
#' @param export A logical argument. If `TRUE`, a *.txt file is exported to
#'   the working directory
#' @param file.name The name of the file if `export = TRUE`
#' @param digits The significant digits to be shown. See
#'   [`tibble::print()`][tibble::formatting] for more details.
#' @param ... Currently not used.
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @method print waas_means
#' @export
#' @examples
#'\donttest{
#' library(metan)
#' data_means <- means_by(data_ge, ENV, GEN)
#' model <- waas_means(data_ge,
#'                     env = ENV,
#'                     gen = GEN,
#'                     resp = everything())
#' print(model)
#' }
print.waas_means <- function(x, export = FALSE, file.name = NULL, digits = 4, ...) {
  if (!class(x) == "waas_means") {
    stop("The object must be of class 'waas_means'")
  }
  if (export == TRUE) {
    file.name <- ifelse(is.null(file.name) == TRUE, "waas_means print", file.name)
    sink(paste0(file.name, ".txt"))
  }
  opar <- options(pillar.sigfig = digits)
  on.exit(options(opar))
  for (i in 1:length(x)) {
    var <- x[[i]]
    cat("Variable", names(x)[i], "\n")
    cat("---------------------------------------------------------------------------\n")
    cat("Weighted average of the absolute scores\n")
    cat("---------------------------------------------------------------------------\n")
    print(var$model)
    cat("---------------------------------------------------------------------------\n")
    cat("Genotype-environment interaction effects\n")
    cat("---------------------------------------------------------------------------\n")
    print(var$ge_eff, digits = digits)
    cat("---------------------------------------------------------------------------\n")
    cat("Proportion of the variance explained\n")
    cat("---------------------------------------------------------------------------\n")
    print(var$proportion, digits = digits)
    cat("---------------------------------------------------------------------------\n")
    cat("Cumulative proportion of the variance explained\n")
    cat("---------------------------------------------------------------------------\n")
    print(var$cum_proportion, digits = digits)
    cat("\n\n\n")
  }
  if (export == TRUE) {
    sink()
  }
}

