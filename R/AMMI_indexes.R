#' AMMI-based stability indexes
#' @description
#' * `ammi_indexes()` `r badge('stable')` computes several AMMI-based stability statistics. See
#' **Details** for a detailed overview.
#' * `AMMI_indexes()` `r badge('deprecated')` use `ammi_indexes()` instead.
#'
#' @details
#'First, let's define some symbols: \mjseqn{N'} is the number of significant
#'interation principal component axis (IPCs) that were retained in the AMMI
#'model via F tests); \mjseqn{\lambda_{n}} is the singular value for th IPC and
#'correspondingly \mjseqn{\lambda_{n}^{2}} its eigen value; \mjseqn{\gamma_{in}}
#'is the eigenvector value for ith genotype; \mjseqn{\delta_{jn}} is the
#'eigenvector value for the th environment. \mjseqn{PC_{1}}, \mjseqn{PC_{2}},
#'and \mjseqn{PC_{n}} are the scores of 1st, 2nd, and nth IPC; respectively;
#'\mjseqn{\theta_{1}}, \mjseqn{\theta_{2}}, and \mjseqn{\theta_{n}}  are
#'percentage sum of squares explained by the 1st, 2nd, and nth IPC,
#'respectively.
#'
#' \loadmathjax
#'
#' *  AMMI Based Stability Parameter (ASTAB) (Rao and Prabhakaran 2005).
#' \mjsdeqn{ASTAB = \sum_{n=1}^{N'}\lambda_{n}\gamma_{in}^{2}}
#'
#' * AMMI Stability Index (ASI) (Jambhulkar et al. 2017)
#' \mjsdeqn{ASI = \sqrt{\left \[ PC_{1}^{2} \times \theta_{1}^{2} \right \]+\left\[ PC_{2}^{2} \times \theta_{2}^{2} \right \]}}
#'
#'
#' * AMMI-stability value (ASV) (Purchase et al., 2000).
#'  \mjsdeqn{ASV_{i}=\sqrt{\frac{SS_{IPCA1}}{SS_{IPCA2}}(\mathrm{IPC} \mathrm{A} 1)^{2}+(\mathrm{IPCA} 2)^{2}}}
#'
#' * Sum Across Environments of Absolute Value of GEI Modelled by AMMI (AVAMGE) (Zali et al. 2012)
#' \mjsdeqn{AV_{(AMGE)} = \sum_{j=1}^{E} \sum_{n=1}^{N'} \left |\lambda_{n}\gamma_{in} \delta_{jn} \right |}
#'
#' * Annicchiarico's D Parameter values (Da) (Annicchiarico 1997)
#' \mjsdeqn{D_{a} = \sqrt{\sum_{n=1}^{N'}(\lambda_{n}\gamma_{in})^2}}
#'
#' * Zhang's D Parameter (Dz) (Zhang et al. 1998)
#' \mjsdeqn{D_{z} = \sqrt{\sum_{n=1}^{N'}\gamma_{in}^{2}}}
#'
#' * Sums of the Averages of the Squared Eigenvector Values (EV) (Zobel 1994)
#' \mjsdeqn{EV = \sum_{n=1}^{N'}\frac{\gamma_{in}^2}{N'}}
#'
#' * Stability Measure Based on Fitted AMMI Model (FA) (Raju 2002)
#' \mjsdeqn{FA = \sum_{n=1}^{N'}\lambda_{n}^{2}\gamma_{in}^{2}}
#'
#' * Modified AMMI Stability Index (MASI) (Ajay et al. 2018)
#' \mjsdeqn{MASI = \sqrt{ \sum_{n=1}^{N'} PC_{n}^{2} \times \theta_{n}^{2}}}
#'
#' * Modified AMMI Stability Value (MASV) (Ajay et al. 2019)
#' \mjsdeqn{MASV = \sqrt{\sum_{n=1}^{N'-1}\left (\frac{SSIPC_{n}}{SSIPC_{n+1}} \right ) \times (PC_{n})^2   + \left (PC_{N'}\right )^2}}
#'
#' * Sums of the Absolute Value of the IPC Scores (SIPC) (Sneller et al. 1997)
#' \mjsdeqn{SIPC = \sum_{n=1}^{N'} | \lambda_{n}^{0.5}\gamma_{in}|}
#'
#'* Absolute Value of the Relative Contribution of IPCs to the Interaction (Za) (Zali et al. 2012)
#' \mjsdeqn{Za = \sum_{i=1}^{N'} | \theta_{n}\gamma_{in} |}
#'
#'* Weighted average of absolute scores (WAAS) (Olivoto et al. 2019)
#' \mjsdeqn{WAAS_i = \sum_{k = 1}^{p} |IPCA_{ik} \times \theta_{k}/ \sum_{k = 1}^{p}\theta_{k}}
#'
#'
#' For all the statistics, simultaneous selection indexes (SSI) are also
#' computed by summation of the ranks of the stability and mean performance,
#' Y_R, (Farshadfar, 2008).
#'
#' @param .data An object of class `waas` or `performs_ammi`
#' @param order.y A vector of the same length of `x` used to order the
#' response variable. Each element of the vector must be one of the `'h'`
#' or `'l'`. If `'h'` is used, the response variable will be ordered
#' from maximum to minimum. If `'l'` is used then the response variable
#' will be ordered from minimum to maximum. Use a comma-separated vector of
#' names. For example, `order.y = c("h, h, l, h, l")`.
#' @param level The confidence level. Defaults to 0.95.
#' @return
#'
#' A list where each element contains the result AMMI-based stability indexes
#' for one variable.
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @references
#' Ajay BC, Aravind J, Abdul Fiyaz R, Bera SK, Kumar N, Gangadhar K, Kona P
#' (2018). “Modified AMMI Stability Index (MASI) for stability analysis.”
#' ICAR-DGR Newsletter, 18, 4–5.
#'
#' Ajay BC, Aravind J, Fiyaz RA, Kumar N, Lal C, Gangadhar K, Kona P, Dagla MC,
#' Bera SK (2019). “Rectification of modified AMMI stability value (MASV).”
#' Indian Journal of Genetics and Plant Breeding (The), 79, 726–731.
#' https://www.isgpb.org/article/rectification-of-modified-ammi-stability-value-masv.
#'
#' Annicchiarico P (1997). “Joint regression vs AMMI analysis of
#' genotype-environment interactions for cereals in Italy.” Euphytica, 94(1),
#' 53–62. \doi{10.1023/A:1002954824178}
#'
#' Farshadfar E (2008) Incorporation of AMMI stability value and grain yield in
#' a single non-parametric index (GSI) in bread wheat. Pakistan J Biol Sci
#' 11:1791–1796. \doi{10.3923/pjbs.2008.1791.1796}
#'
#' Jambhulkar NN, Rath NC, Bose LK, Subudhi HN, Biswajit M, Lipi D, Meher J
#' (2017). “Stability analysis for grain yield in rice in demonstrations
#' conducted during rabi season in India.” Oryza, 54(2), 236–240.
#' \doi{10.5958/2249-5266.2017.00030.3}
#'
#' Olivoto T, LUcio ADC, Silva JAG, et al (2019) Mean Performance and Stability
#' in Multi-Environment Trials I: Combining Features of AMMI and BLUP
#' Techniques. Agron J 111:2949–2960. \doi{10.2134/agronj2019.03.0220}
#'
#' Raju BMK (2002). “A study on AMMI model and its biplots.” Journal of the
#' Indian Society of Agricultural Statistics, 55(3), 297–322.
#'
#' Rao AR, Prabhakaran VT (2005). “Use of AMMI in simultaneous selection of
#' genotypes for yield and stability.” Journal of the Indian Society of
#' Agricultural Statistics, 59, 76–82.
#'
#' Sneller CH, Kilgore-Norquest L, Dombek D (1997). “Repeatability of yield
#' stability statistics in soybean.” Crop Science, 37(2), 383–390.
#' \doi{10.2135/cropsci1997.0011183X003700020013x}
#'
#' Zali H, Farshadfar E, Sabaghpour SH, Karimizadeh R (2012). “Evaluation of
#' genotype × environment interaction in chickpea using measures of stability
#' from AMMI model.” Annals of Biological Research, 3(7), 3126–3136.
#'
#' Zhang Z, Lu C, Xiang Z (1998). “Analysis of variety stability based on AMMI
#' model.” Acta Agronomica Sinica, 24(3), 304–309.
#' http://zwxb.chinacrops.org/EN/Y1998/V24/I03/304.
#'
#' Zobel RW (1994). “Stress resistance and root systems.” In Proceedings of the
#' Workshop on Adaptation of Plants to Soil Stress. 1-4 August, 1993. INTSORMIL
#' Publication 94-2, 80–99. Institute of Agriculture and Natural Resources,
#' University of Nebraska-Lincoln.
#' @name ammi_indexes
#' @export
#'
#' @examples
#'\donttest{
#' library(metan)
#' model <-
#'   performs_ammi(data_ge,
#'                 env = ENV,
#'                 gen = GEN,
#'                 rep = REP,
#'                 resp = c(GY, HM))
#' model_indexes <- ammi_indexes(model)
#'
#'
#' # Alternatively (and more intuitively) using %>%
#' # If resp is not declared, all traits are analyzed
#' res_ind <- data_ge %>%
#'            performs_ammi(ENV, GEN, REP, verbose = FALSE) %>%
#'            ammi_indexes()
#'
#' rbind_fill_id(res_ind, .id = "TRAIT")
#'}
#'
ammi_indexes <- function(.data, order.y = NULL, level = 0.95) {
    if(!is.null(order.y)){
        order.y <- unlist(strsplit(order.y, split = ", "))
    } else {
        order.y <- rep("h", length(.data))
    }
    if (!class(.data)  %in% c("waas", "performs_ammi")) {
        stop("The object 'x' must be an object of class \"waas\" or \"performs_ammi\"")
    }
    if (any(!order.y %in% c("h", "l")) == TRUE) {
        stop("The argument 'order.y' must be a comma-separated vector with 'h' or 'l'. Did you accidentally omit the space between the comma and the following word?")
    }
    if (length(order.y) != length(.data)) {
        stop("The lenght of argument 'order.y' must be ", length(.data), ", the length of '.data'")
    }
    listres <- list()
    varin <- 1
    for (var in 1:length(.data)) {
        model <- .data[[var]]
        n <- sum(model$PCA$`Pr(>F)` <= (1 - level), na.rm = TRUE)
        n <- ifelse(n == 0, 1, n)
        meange <- model$MeansGxE
        effects <- residuals(lm(Y ~ ENV + GEN, data = meange))
        meange$residual <- effects
        ge <- by(meange[, 7], meange[, c(2, 1)], function(x) sum(x, na.rm = TRUE))
        ge <- array(ge, dim(ge), dimnames(ge))
        ngen <- nrow(ge)
        nenv <- ncol(ge)
        svdge <- svd(ge)
        gamma.n <- svdge$u[, 1:n]
        delta.n <- svdge$v[, 1:n]
        lambda  <- svdge$d[1:n]
        theta.n <- model$PCA$Proportion[1:n]/100
        PCA <- data.frame(model$model)
        SCOR <- as.matrix(PCA[PCA[, 1] == "GEN", c(seq(4, n + 3))])
        mean <- PCA[PCA[, 1] == "GEN", c(2:3)]

        if (length(order.y) == 1) {
            if (order.y == "h") {
                rY <- rank(-mean[2])
            }
            if (order.y == "l") {
                rY <- rank(mean[2])
            }
        } else {
            if (order.y[[varin]] == "h") {
                rY <- rank(-mean[2])
            }
            if (order.y[[varin]] == "l") {
                rY <- rank(mean[2])
            }
        }
        varin <- varin + 1
        # sum of absolute scores
        if (n == 1) {
            SIPC <- abs(SCOR)
        } else {
            SIPC <- unname(rowSums(apply(SCOR, 2, FUN = abs)))
        }
        rS <- rank(SIPC)
        ssiSIPC <- rS + rY

        # Za
        if (n == 1) {
            Za <- abs(gamma.n * theta.n)
        } else {
            Za <- rowSums(abs(gamma.n %*% diag(theta.n)))
        }
        rZA <- rank(Za)
        ssiZA <- rZA + rY

        # averages of the squared eigenvector
        if (n == 1) {
            EV <- gamma.n^2/n
        } else {
            EV <- rowSums(gamma.n^2/n)
        }
        rEV <- rank(EV)
        ssiEV <- rEV + rY

        # AMMI stability values
        pc <- model$PCA$`Sum Sq`[1]/model$PCA$`Sum Sq`[2]
        SCOR2 <- PCA[PCA[, 1] == "GEN", c(seq(4, 2 + 3))]
        ASV <- sqrt((pc * SCOR2[, 1])^2 + SCOR2[, 2]^2)
        rASV <- rank(ASV)
        ssiASV <- rASV + rY

        # modified AMMI stability values
        ssquares <- model$ANOVA[c(5:(5 + n)), 3]
        MASV <- rep(0, nrow(SCOR))
        for (i in 1:ncol(SCOR)) {
            pc <- ssquares[i]/ssquares[i + 1]
            MASV <- MASV + (SCOR[, i] * pc)^2
            if ((i + 1) == ncol(SCOR))
                (break)()
        }
        MASV <- sqrt(MASV + (SCOR[, ncol(SCOR)]^2))
        rMASV <- rank(MASV)
        ssiMASV <- rMASV + rY


        # DI
        DZ <- sqrt(rowSums((gamma.n)^2))
        rDZ <- rank(DZ)
        ssiDZ<- rDZ + rY


        # FA
        FA <- rowSums((gamma.n^2) %*% diag(lambda^2))
        rFA <- rank(FA)
        ssiFA <- rFA + rY

        # Annicchiarico's Da
        DA <- sqrt(rowSums((gamma.n %*% diag(lambda))^2))
        rDA <- rank(DA)
        ssiDA <- rDA + rY


        # ASTAB
        ASTAB <- rowSums((gamma.n^2) %*% diag(lambda))
        rASTAB <- rank(ASTAB)
        ssiASTAB <- rASTAB + rY


        # ASI
        ASI <- sqrt((SCOR[,1]^2 * theta.n[1]^2) + (SCOR[,2]^2 * theta.n[2]^2))
        rASI <- rank(ASI)
        ssiASI <- rASI + rY


        # MASI
        MASI <- sqrt(rowSums(SCOR^2 %*% diag(theta.n^2)))
        rMASI <- rank(MASI)
        ssiMASI <- rMASI + rY

        # AVAMGE
        ge.n <- gamma.n %*% diag(lambda) %*% t(delta.n)
        AVAMGE <- rowSums(apply(ge.n, 2, FUN = abs))
        rAVAMGE <- rank(AVAMGE)
        ssiAVAMGE <- rAVAMGE + rY

        # Weighted average of absolute scores
        if(n == 1){
            explan <- model$PCA[1, ][7]
        } else{
            explan <- model$PCA[which(model$PCA[6] < 1 - level),][7]
        }
        WAAS <-
            SCOR %>%
            abs() %>%
            t() %>%
            as.data.frame()
        WAAS <- sapply(WAAS, weighted.mean, w = explan$Proportion)
        rWAAS <- rank(WAAS)
        ssiWAAS <- rWAAS + rY

        temp <- tibble(
            GEN = mean[1] %>% pull(),
            Y = mean[2] %>% pull(),
            Y_R = rY,
            ASTAB = ASTAB,
            ASTAB_R = rASTAB,
            ssiASTAB = ssiASTAB,
            ASI = ASI,
            ASI_R = rASI,
            ASI_SSI = ssiASI,
            ASV = ASV,
            ASV_R = rASV,
            ASV_SSI = ssiASV,
            AVAMGE = AVAMGE,
            AVAMGE_R = rAVAMGE,
            AVAMGE_SSI = ssiAVAMGE,
            DA = DA,
            DA_R = rDA,
            DA_SSI = ssiDA,
            DZ = DZ,
            DZ_R = rDZ,
            DZ_SSI = ssiDZ,
            EV = EV,
            EV_R = rEV,
            EV_SSI = ssiEV,
            FA = FA,
            FA_R = rFA,
            FA_SSI = ssiFA,
            MASI = MASI,
            MASI_R = rMASI,
            MASI_SSI = ssiMASI,
            MASV = MASV,
            MASV_R = rMASV,
            MASV_SSI = ssiMASV,
            SIPC = SIPC,
            SIPC_R = rS,
            SIPC_SSI = ssiSIPC,
            ZA = Za,
            ZA_R = rZA,
            ZA_SSI = ssiZA,
            WAAS = WAAS,
            WAAS_R = rWAAS,
            WAAS_SSI = ssiWAAS)
        listres[[paste(names(.data[var]))]] <- temp
    }
    invisible(structure(listres, class = "ammi_indexes"))
}

#' @name ammi_indexes
#' @export
AMMI_indexes <- function(.data, order.y = NULL, level = 0.95) {
    warning("`AMMI_indexes()` is deprecated as of metan 1.16.0. use `ammi_indexes()` instead.")
    ammi_indexes(.data, order.y, level)
}



#' Print an object of class ammi_indexes
#'
#' Print the `ammi_indexes` object in two ways. By default, the results are shown
#' in the R console. The results can also be exported to the directory into a
#' *.txt file.
#'
#' @param x An object of class `ammi_indexes`.
#' @param which Which should be printed. Defaults to `"stats"`. Other
#'   possible values are `"ranks"` for genotype ranking and `"ssi"`
#'   for the simultaneous selection index.
#' @param export A logical argument. If `TRUE`, a *.txt file is exported to
#'   the working directory.
#' @param file.name The name of the file if `export = TRUE`
#' @param digits The significant digits to be shown.
#' @param ... Options used by the tibble package to format the output. See
#'   [`tibble::print()`][tibble::formatting] for more details.
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @method print ammi_indexes
#' @export
#' @examples
#' \donttest{
#' library(metan)
#' model <- performs_ammi(data_ge, ENV, GEN, REP, GY) %>%
#'          ammi_indexes()
#' print(model)
#' }
print.ammi_indexes <- function(x, which = "stats", export = FALSE, file.name = NULL, digits = 3, ...) {
    if (!class(x) == "ammi_indexes") {
        stop("The object must be of class 'ammi_indexes'")
    }
    opar <- options(pillar.sigfig = digits)
    on.exit(options(opar))
    if (export == TRUE) {
        file.name <- ifelse(is.null(file.name) == TRUE, "ammi_indexes print", file.name)
        sink(paste0(file.name, ".txt"))
    }
    if(!which %in% c("stats", "ssi", "ranks")){
        stop("Argument 'which' must be one of 'stats', 'ranks', or 'ssi'", call. = FALSE)
    }
    for (i in 1:length(x)) {
        if(which == "stats"){
            var <- x[[i]] %>%
                select_cols(-contains("_R"), -contains("_SSI"))
        }
        if(which == "ranks"){
            var <- x[[i]] %>%
                select_cols(GEN, contains("_R"))
        }
        if(which == "ssi"){
            var <- x[[i]] %>%
                select_cols(GEN, contains("SSI"))
        }
        cat("Variable", names(x)[i], "\n")
        cat("---------------------------------------------------------------------------\n")
        cat("AMMI-based stability indexes\n")
        cat("---------------------------------------------------------------------------\n")
        print(var)
    }
    if (export == TRUE) {
        sink()
    }
}

