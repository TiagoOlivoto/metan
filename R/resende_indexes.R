Resende_indexes <- function(.data) {
    if (!is(.data, "WAASB")) {
        stop("The object '.data' must be an object of class \"WAASB\"")
    }

    # Helper functions
    hmean_fun <- function(x) {
        hmean <- length(x)/sum(1/x)
        return(hmean)
    }

    listres <- list()

    for (var in 1:length(.data)) {

        gge <- .data[[var]][["BLUPgge"]]
        # Harmonic mean
        GEPRED = make_mat(gge, GEN, ENV, Predicted)
        HMGV <- data.frame(HMGV = apply(GEPRED, 1, FUN = hmean_fun), HMGV_order = rank(-apply(GEPRED,
            1, FUN = hmean_fun)))
        ## Relative performance
        y <- .data[[var]][["MeansGxE"]]

        GEMEAN <- make_mat(y, GEN, ENV, Y)
        ovmean <- mean(y$Y)
        mean_env <- apply(GEMEAN, 2, FUN = mean)
        RPGV <- data.frame(RPGV = apply(t(t(GEPRED)/mean_env), 1, mean))
        RPGV_data <- dplyr::mutate(RPGV, GY_RPGV = RPGV * ovmean, RPGV_order = rank(-GY_RPGV))
        ## Harmonic mean of Relative performance
        HMRPGV <- data.frame(HMRPGV = apply(t(t(GEPRED)/mean_env), 1, hmean_fun))
        HMRPGV_data <- dplyr::mutate(HMRPGV, GY_HMRPGV = HMRPGV * ovmean, HMRPGV_order = rank(-GY_HMRPGV))

        temp <- data.frame(cbind(HMGV, RPGV_data, HMRPGV_data))
        rownames(temp) <- rownames(RPGV)
        listres[[paste(names(.data[var]))]] <- temp
    }
    invisible(listres)
}
