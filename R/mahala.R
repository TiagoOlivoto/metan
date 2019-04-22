mahala = function(.means, covar, inverted = FALSE){
  cb = data.frame(t(combn(nrow(.means), 2)))
  cb = cb[with(cb, order(X2)), ]
  nvars = data.frame(matrix(nrow = nrow(cb), ncol = ncol(.means)))
  for (i in 1:nrow(cb)){
    nvars[i ,] = .means[cb[i, 1], ] - .means[cb[i, 2],]
  }
  mmean = as.matrix(nvars)
  maha = matrix(ncol = nrow(.means), nrow = nrow(.means))
  diag(maha) = 0
  if (inverted == TRUE){
    invmat = as.matrix(covar)
  } else {
    invmat = solve(covar)
  }
  dists = diag(mmean %*% invmat %*% t(mmean))
  maha[upper.tri(maha, diag = F)] <- maha[lower.tri(maha, diag = F)] <- dists
  rownames(maha) <- colnames(maha) <- rownames(.means)
  return(maha)
}
