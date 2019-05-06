make_sym <- function(.matrix, make = "upper") {
  if(make == "upper"){
    .matrix[upper.tri(.matrix)] <- t(.matrix)[upper.tri(.matrix)]
  }
  if (make == "lower"){
    .matrix[lower.tri(.matrix)] <- t(.matrix)[lower.tri(.matrix)]
  }
  return(.matrix)
}
