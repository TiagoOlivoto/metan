rbind_fill <- function(..., fill = "-") {
    df <- list(...)
    ns <- unique(unlist(sapply(df, names)))
    do.call(rbind, lapply(df, function(x) {
        for (n in ns[!ns %in% names(x)]) {
            x[[n]] <- fill
        }
        x
    }))
}
