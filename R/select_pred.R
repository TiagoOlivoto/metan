#' Selects a best subset of predictor variables.
#'
#' Selects among a set of `covariates` the best set of `npred`
#' predictors for a given response trait `resp` based on AIC values.
#'
#' @param .data A data frame with the response variable and covariates.
#' @param resp The response variable.
#' @param covariates The covariates. Defaults to *NULL*. In this case, all
#'   numeric traits in `.data`, except that in `resp` are selected. To
#'   select specific covariates from `.data`, use a list of unquoted
#'   comma-separated variable names (e.g. *traits = c(var1, var2, var3)*),
#'   an specific range of variables, (e.g. *traits = c(var1:var3)*), or
#'   even a select helper like `starts_with("N")`.
#' @param npred An integer specifying the size of the subset of predictors to be
#'   selected
#'
#' @return A list with the following elements:
#' * **sel_mod** An object of class `lm` that is the selected model.
#' * **predictors** The name of the selected predictors.
#' * **AIC** The Akaike's Information Criterion for the selected model.
#' * **pred_models** The Akaike's Information Criterion and the predictors
#' selected in each step.
#' * **predicted** The predicted values considering the model in
#' `sel_mod`.
#' @md
#' @export
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @examples
#' \donttest{
#' library(metan)
#' mod <- select_pred(data_ge2, resp = PH, npred = 10)
#' mod$predictors
#' mod$AIC
#' }
select_pred <- function (.data,
                         resp,
                         covariates = NULL,
                         npred){
  y <- select(.data, {{resp}})
  if(ncol(y) > 1){
    stop("Only one response trait is allowed.")
  }
  if(missing(covariates)){
    x <-
      select(.data %>% select_numeric_cols(), everything(), -names(y)) %>%
      as.data.frame()
  } else{
    x <-
      select(.data %>% select_numeric_cols(), {{covariates}}, -names(y)) %>%
      as.data.frame()
  }
  y <- pull(y)
  npreds <- ncol(x)
  varin <- integer(npred)
  n <- length(y)
  if (npred >= npreds) {
    stop("The number of predictors ('npred') must be lesser than the number of covariates.")
  }
  mod <- lm(y ~ NULL)
  x_nam <- names(x)
  x_yes <- names(x)
  f2 <-
    update(
      as.formula(mod),
      paste(". ~ ", paste(x_yes, collapse = "+"))
    )
  step_help <- function(j) {
    f1 <- as.formula(mod, env = environment(fun = NULL))
    f2 <- update(f1, . ~ . + x[, j])
    mods <- lm(f2)
    return(AIC(mods))
  }
  out <- 1:npreds
  x_yes <- NULL
  aic_step <- list()
  pred_step <- list()
  for (k in 1:npred) {
    ic <- sapply(out, step_help)
    var_out <- which.min(ic)
    varin[k] <- out[var_out]
    out <- out[-var_out]
    x_nam <- paste("x[,", varin[[k]], "]", sep = "")
    x_yes[k] <- x_nam
    f2 <-
      update(
        as.formula(mod),
        paste(". ~ ", paste(x_yes, collapse = "+"))
      )
    mod <- lm(f2)
    best_ic <- AIC(mod)
    pred_step[[paste("Model", k)]] <- names(x[varin])
    aic_step[[paste("Model", k)]] <- best_ic
  }
  best_models <-
    sapply(aic_step, function(x){x}) %>%
    as.data.frame() %>%
    rownames_to_column("Model") %>%
    set_names("Model", "AIC") %>%
    add_cols(Predictors = pred_step) %>%
    arrange(AIC)
  res <- list(sel_mod = mod,
              predictors = names(x[varin]),
              AIC = AIC(mod),
              pred_models = best_models,
              predicted = predict(mod, type = "response")) %>%
    set_class("sel_pred")
  return(res)
}
