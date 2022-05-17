#' Smith-Hazel index
#'
#' @description
#' `r badge('stable')`
#'
#' Computes the Smith (1936) and Hazel (1943) index given economic weights and
#' phenotypic and genotypic variance-covariance matrices. The Smith-Hazel index
#' is computed as follows:
#'\loadmathjax
#'\mjsdeqn{\bf{b = P^{-1}Aw}}
#'
#' where \mjseqn{\bf{P}} and \mjseqn{\bf{G}} are phenotypic and genetic
#' covariance matrices, respectively, and \mjseqn{\bf{b}} and \mjseqn{\bf{w}}
#' are vectors of index coefficients and economic weightings, respectively.
#'
#' The genetic worth \mjseqn{I} of an individual
#' genotype based on traits *x*, *y*, ..., *n*, is calculated as:
#'
#'\mjsdeqn{I = b_xG_x + b_yG_y + ... + b_nG_n}
#'
#' where *b* the index coefficient for the traits *x*, *y*, ...,
#' *n*, respectively, and *G* is the individual genotype BLUPs for the
#' traits *x*, *y*, ..., *n*, respectively.
#'
#' @details
#' When using the phenotypic means in `.data`, be sure the genotype's code
#' are in rownames. If `.data` is an object of class `gamem` them the
#' BLUPs for each genotype are used to compute the index. In this case, the
#' genetic covariance components are estimated by mean cross products.
#'
#' @param .data The input data. It can be either a two-way table with genotypes
#'   in rows and traits in columns, or an object fitted with the function
#'   [gamem()]. Please, see **Details** for more details.
#' @param use_data Define which data to use If `.data` is an object of
#'   class `gamem`. Defaults to `"blup"` (the BLUPs for genotypes).
#'   Use `"pheno"` to use phenotypic means instead BLUPs for computing the
#'   index.
#' @param pcov,gcov The phenotypic and genotypic variance-covariance matrix,
#'   respectively. Defaults to `NULL`. If a two-way table is informed in
#'   `.data` these matrices are mandatory.
#' @param SI The selection intensity (percentage). Defaults to `20`
#' @param weights The vector of economic weights. Defaults to a vector of 1s
#'   with the same length of the number of traits.
#'
#' @references
#' Smith, H.F. 1936. A discriminant function for plant selection. Ann. Eugen.
#' 7:240-250.
#' \doi{10.1111/j.1469-1809.1936.tb02143.x}
#'
#' Hazel, L.N. 1943. The genetic basis for constructing selection indexes.
#' Genetics 28:476-90. https://www.genetics.org/content/28/6/476.short
#'
#' @return An object of class `hz` containing:
#'  * **b**: the vector of index coefficient.
#'  * **index**: The genetic worth.
#'  * **sel_dif_trait**: The selection differencial.
#'  * **sel_gen**: The selected genotypes.
#'  * **gcov**: The genotypic variance-covariance matrix
#'  * **pcov**: The phenotypic variance-covariance matrix
#'
#' @export
#' @md
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @seealso [mtsi()], [mgidi()], [fai_blup()]
#' @examples
#' \donttest{
#' vcov <- covcor_design(data_g, GEN, REP, everything())
#' means <- as.matrix(vcov$means)
#' pcov <- vcov$phen_cov
#' gcov <- vcov$geno_cov
#'
#' index <- Smith_Hazel(means, pcov = pcov, gcov = gcov, weights = rep(1, 15))
#' }
Smith_Hazel <- function(.data,
                        use_data = "blup",
                        pcov = NULL,
                        gcov = NULL,
                        SI = 15,
                        weights = NULL){
  if(!use_data %in% c("blup", "pheno")){
    stop("Argument 'use_data = ", match.call()["use_data"], "'", "invalid. It must be either 'blup' or 'pheno'.")
  }
  if(has_class(.data, "gamem")){
    ifelse(is.null(weights),
           weights <- rep(1, length(.data)),
           weights <- weights)
    mat <-
      gmd(.data, ifelse(use_data == "blup", "blupg", "data"), verbose = FALSE) %>%
      means_by(GEN) %>%
      column_to_rownames("GEN") %>%
      as.matrix()
    pcov <-
      gmd(.data, "data", verbose = FALSE) %>%
      means_by(GEN) %>%
      column_to_rownames("GEN") %>%
      as.matrix() %>%
      cov()
    gcov <- gmd(.data, "gcov", verbose = FALSE)
    gcov <- gcov[rownames(pcov), rownames(pcov)]
    h2 <- gmd(.data, "h2", verbose = FALSE)
  } else{
    ifelse(missing(weights), weights <- rep(1, ncol(.data)), weights <- weights)
    if(missing(pcov) || missing(gcov)){
      stop("The genotypic and phenotypic covariance matrices must be informed.", call. = FALSE)
    } else{
      if(!isSymmetric(pcov) || !isSymmetric(gcov)){
        stop("The genotypic and phenotypic covariance matrices must be symmetric.", call. = FALSE)
      } else{
        pcov <- as.matrix(pcov)
        gcov <- as.matrix(gcov)
      }
    }
    mat <- .data

  }
  ngs <- round(nrow(mat) * (SI/100), 0)
  b <- solve_svd(pcov) %*% gcov %*% weights
  index <-
    mat %*% b %>%
    as.data.frame() %>%
    rownames_to_column("GEN") %>%
    arrange(-V1)
  sel_gen <- head(index, ngs)[[1]]
  if(length(sel_gen)==1){
    xsel <- t(as.matrix(mat[sel_gen, ]))
  } else{
    xsel <- mat[sel_gen, ]
  }
  sel_dif_trait <-
    tibble(VAR = colnames(mat),
           Xo = colMeans(mat),
           Xs = colMeans(xsel),
           SD = Xs - Xo,
           SDperc = (Xs - Xo) / Xo * 100)
  vars <- tibble(VAR = colnames(mat),
                 sense = weights) %>%
    mutate(sense = ifelse(sense < 0, "decrease", "increase"))
  if(has_class(.data, "gamem")){
    sel_dif_trait <-
      left_join(sel_dif_trait, h2, by = "VAR") %>%
      add_cols(SG = SD * h2,
               SGperc = SG / Xo * 100)
  }
  sel_dif_trait <-
    sel_dif_trait %>%
    left_join(vars, by = "VAR") %>%
    mutate(goal = case_when(
      sense == "decrease" & SDperc < 0  |  sense == "increase" & SDperc > 0 ~ 100,
      TRUE ~ 0
    ))
  total_gain <-
    desc_stat(sel_dif_trait,
              by = sense,
              any_of(c("SDperc", "SGperc")),
              stats = c("min, mean, max, sum"))
  b <-
    data.frame(b) %>%
    rownames_to_column("VAR") %>%
    add_cols(gen_weights = weights)



  results <-
    list(b = b,
         index = index,
         sel_dif_trait = sel_dif_trait,
         total_gain = total_gain,
         sel_gen = sel_gen,
         gcov = gcov,
         pcov = pcov)
  return(results %>% set_class("sh"))
}


#' Plot the Smith-Hazel index
#'
#' Makes a radar plot showing the individual genetic worth for the Smith-Hazel index
#'
#'
#' @param x An object of class `sh`
#' @param SI An integer (0-100). The selection intensity in percentage of the
#'   total number of genotypes.
#' @param radar Logical argument. If true (default) a radar plot is generated
#'   after using `coord_polar()`.
#' @param arrange.label Logical argument. If `TRUE`, the labels are
#'   arranged to avoid text overlapping. This becomes useful when the number of
#'   genotypes is large, say, more than 30.
#' @param size.point The size of the point in graphic. Defaults to 2.5.
#' @param size.line The size of the line in graphic. Defaults to 0.7.
#' @param size.text The size for the text in the plot. Defaults to 10.
#' @param col.sel The colour for selected genotypes. Defaults to `"red"`.
#' @param col.nonsel The colour for nonselected genotypes. Defaults to `"black"`.
#' @param ... Other arguments to be passed from ggplot2::theme().
#' @return An object of class `gg, ggplot`.
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @method plot sh
#' @export
#' @examples
#' \donttest{
#' library(metan)
#' vcov <- covcor_design(data_g, GEN, REP, everything())
#' means <- as.matrix(vcov$means)
#' pcov <- vcov$phen_cov
#' gcov <- vcov$geno_cov
#'
#' index <- Smith_Hazel(means, pcov = pcov, gcov = gcov, weights = rep(1, 15))
#' plot(index)
#'}
#'
plot.sh <- function(x,
                    SI = 15,
                    radar = TRUE,
                    arrange.label = FALSE,
                    size.point = 2.5,
                    size.line = 0.7,
                    size.text = 10,
                    col.sel = "red",
                    col.nonsel = "black",
                    ...) {
    data <- x$index %>% add_cols(sel = "Selected")
    data[["sel"]][(round(nrow(data) * (SI/100), 0) + 1):nrow(data)] <- "Nonselected"
    cutpoint <- min(subset(data, sel == "Selected")$V1)
    p <-
      ggplot(data = data, aes(x = reorder(GEN, V1), y = V1)) +
      geom_hline(yintercept = cutpoint, col = col.sel, size = size.line) +
      geom_path(colour = "black", group = 1, size = size.line) +
      geom_point(size = size.point, aes(fill = sel), shape = 21, colour = "black", stroke  = size.point / 10) +
      scale_x_discrete() +
      theme_minimal() +
      theme(legend.position = "bottom",
            legend.title = element_blank(),
            axis.title.x = element_blank(),
            panel.border = element_blank(),
            axis.text = element_text(colour = "black"),
            text = element_text(size = size.text)) +
      labs(y = "Individual genetic worth") +
      scale_fill_manual(values = c(col.nonsel, col.sel))
    if (radar == TRUE) {
      if(arrange.label == TRUE){
        tot_gen <- length(unique(data$GEN))
        fseq <- c(1:(tot_gen/2))
        sseq <- c((tot_gen/2 + 1):tot_gen)
        fang <- c(90 - 180/length(fseq) * fseq)
        sang <- c(-90 - 180/length(sseq) * sseq)
        p <- p + coord_polar() +
          theme(axis.text.x = element_text(angle = c(fang, sang)), legend.margin = margin(-120, 0, 0, 0), ...)
      } else{
        p <- p + coord_polar()
      }
    }
  return(p)
}





#' Print an object of class sh
#'
#' Print a `sh` object in two ways. By default, the results are shown in
#' the R console. The results can also be exported to the directory.
#'
#' @param x An object of class `sh`.
#' @param export A logical argument. If `TRUE`, a *.txt file is exported
#'   to the working directory
#' @param file.name The name of the file if `export = TRUE`
#' @param digits The significant digits to be shown.
#' @param ... Options used by the tibble package to format the output. See
#'   [`tibble::print()`][tibble::formatting] for more details.
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @method print sh
#' @export
#' @examples
#' \donttest{
#' vcov <- covcor_design(data_g, GEN, REP, everything())
#' means <- as.matrix(vcov$means)
#' pcov <- vcov$phen_cov
#' gcov <- vcov$geno_cov
#'
#' index <- Smith_Hazel(means, pcov = pcov, gcov = gcov, weights = rep(1, 15))
#' print(index)
#' }
print.sh <- function(x, export = FALSE, file.name = NULL, digits = 4, ...) {
  if (export == TRUE) {
    file.name <- ifelse(is.null(file.name) == TRUE, "Smith-Hazel print", file.name)
    sink(paste0(file.name, ".txt"))
  }
  opar <- options(pillar.sigfig = digits)
  on.exit(options(opar))
  cat("\n-----------------------------------------------------------------------------------\n")
  cat("Index coefficients\n")
  cat("-----------------------------------------------------------------------------------\n")
  print(x$b)
  cat("\n-----------------------------------------------------------------------------------\n")
  cat("Genetic worth\n")
  cat("-----------------------------------------------------------------------------------\n")
  print(x$index)
  cat("\n-----------------------------------------------------------------------------------\n")
  cat("Selection gain\n")
  cat("-----------------------------------------------------------------------------------\n")
  print(x$sel_dif_trait)
  cat("\n-----------------------------------------------------------------------------------\n")
  cat("Phenotypic variance-covariance matrix\n")
  cat("-----------------------------------------------------------------------------------\n")
  print(x$pcov, digits = 2)
  cat("\n-----------------------------------------------------------------------------------\n")
  cat("Genotypic variance-covariance matrix\n")
  cat("-----------------------------------------------------------------------------------\n")
  print(x$gcov, digits = 2)
  if (export == TRUE) {
    sink()
  }
}
