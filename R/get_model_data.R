#' Easily get data from a model
#'
#' Easily get data from a model of class \code{waas} or \code{waasb} (Olivoto et al., 2019a, 2019b)
#'
#'
#' @param x An object created with the functions \code{\link{waas}} or \code{\link{waasb}}.
#' @param what What should be captured from the model. The following options are allowed:
#' \code{"blupg"} for genotype's predicted mean and \code{"blupge"} for genotype-vs-environment's
#' predicted mean (if the object is of class \code{waasb} only); \code{"Y"} for raw means (Default);
#' \code{"PCn"} for the scores of the nth Interaction Principal Component Axis (IPCA),
#'  where \code{n} is the desired IPCA (e.g. \code{"PC3"}); \code{"WAASB"} for the weighted
#'  average of the absolute scores; \code{"PctResp"} and \code{"PctWAASB"} for the reescaled
#'  values of the response variable and WAASB, respectively; \code{"wResp"} and \code{"wWAASB"}
#'  for weight of the response variable and stability, respectively; \code{"OrResp"} and
#'  \code{"OrWAASB"} for the ranking regarding the response variable and WAASB, respectively;
#'  \code{"OrPC1"} for the ranking regarding the first principal component axix; \code{"WAASBY"}
#'  for the superiority index, \code{"OrWAASBY"} for the ranking regarding the superiorirty
#'  index.
#' @param type Chose if the statistics must be show by genotype (\code{type = "GEN"}, default)
#' or environment (\code{type = "ENV"}), when possible.
#' @importFrom dplyr starts_with matches case_when
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @export
#' @references Olivoto, T., A.D.C. L{\'{u}}cio, J.A.G. da silva, V.S. Marchioro,
#'  V.Q. de Souza, and E. Jost. 2019a. Mean performance and stability in multi-environment
#'   trials I: Combining features of AMMI and BLUP techniques. Agron. J. doi:10.2134/agronj2019.03.0220.
#' @references Olivoto, T., A.D.C. L{\'{u}}cio, J.A.G. da silva, B.G. Sari, and M.I. Diel. 2019b.
#'  Mean performance and stability in multi-environment trials II: Selection based on multiple
#'   traits. Agron. J.doi:10.2134/agronj2019.03.0221.
#' @examples
#'
#' library(metan)
#'
#' model = waas(data_ge2,
#'              env = ENV,
#'              gen = GEN,
#'              rep = REP,
#'              resp = c(PH, ED, TKW, NKR))
#' get_model_data(model, what = "OrWAASBY")
#'
#'
#' waasb(data_ge2, ENV, GEN, REP,
#'       resp = c(PH, ED, TKW, NKR)) %>%
#' get_model_data(what = "blupg")
#'
#'

get_model_data <- function(x, what = "Y", type = "GEN") {
if(!class(x) %in% c("waasb", "waas")){
  stop("Invalid input in 'x' argument. It must be one object of class 'waas' or 'waasb'")
}
  if(substr(what, 1, 2) == "PC"){
    npc = ncol(x[[1]][["model"]] %>%
                 select(starts_with("PC")) %>%
                 select(matches("PC\\d+")))
  npcwhat = as.numeric(substr(what, 3, nchar(what)))
  if(npcwhat > npc){
    stop("The number of principal components informed seems to be larger than those in the model informed in 'x'.")
  }
}
check = c("blupg", "blupge", "Y", "WAASB", "PctResp",
          "PctWAASB", "wRes", "wWAASB", "OrResp", "OrWAASB",
          "OrPC1", "WAASBY", "OrWAASBY")
check2 = paste("PC", 1:200, sep = "")
if(what %in% c("blupg", "blupge") && class(x) == "waas"){
  stop("Blups can only be captured from a model fitted with the function 'waasb'.")
}
if(!what %in% c(check, check2)){
  stop("The argument 'what' is invalid. Please, check the function help (?get_model_data) for more details.")
}

if(!type %in% c("GEN", "ENV")){
  stop("Argument 'type' invalid. It must be either 'GEN' or 'ENV'.")
}
  if(what %in% c("WAASB", "PctWAASB", "wWAASB", "OrWAASB", "WAASBY", "OrWAASBY") && class(x) == "waas"){
  what =  case_when(what == "WAASB" ~ "WAAS",
                    what == "PctWAASB" ~ "PctWAAS",
                    what == "wWAASB" ~ "wWAAS",
                    what == "OrWAASB" ~ "OrWAAS",
                    what == "WAASBY" ~ "WAASY",
                    what == "OrWAASBY" ~ "OrWAASY")
  }
if(what %in% c("blupg", "blupge")){
  if(what == "blupg"){
  datt <- x[[1]][["blupGEN"]]
  bind <- do.call(cbind, lapply(x, function(x) {
    val <- arrange(x[["blupGEN"]], GEN)$Predicted
  }))  %>%
    as_tibble() %>%
    mutate(gen = datt %>% arrange(GEN) %>% pull(GEN)) %>%
    select(gen, everything())
  }
  if(what == "blupge"){
    bind <- do.call(cbind, lapply(x, function(x) {
      val <- x[["BLUPgge"]][["Predicted"]]
    })) %>%
      as_tibble() %>%
      mutate(ENV = x[[1]][["BLUPgge"]][["ENV"]],
             GEN = x[[1]][["BLUPgge"]][["GEN"]]) %>%
      select(ENV, GEN, everything())
  }
} else{
  bind = do.call(cbind,
                 lapply(x, function(x){
                   x$model[[what]]
                 })) %>%
    as_tibble() %>%
    mutate(gen = x[[1]][["model"]][["Code"]],
           type = x[[1]][["model"]][["type"]]) %>%
    dplyr::filter(type == {{type}}) %>%
    select(-type) %>%
    select(gen, everything())
}
  return(bind)
}

