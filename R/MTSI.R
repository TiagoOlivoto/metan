MTSI = function(x, index = "WAASB", show = TRUE, SI = 15, mineval = 1){
  if (length(x) == 1){
    stop("The multitrait stability index cannot be computed with one single variable.")
  }
  if(index == "WAASBY"){
    ideotype.D = rep(100, length(x))
  }

  if (class(x) == "WAAS.AMMI"){
    if (index == "WAASB"){
      bind = data.frame(do.call(cbind,
                                lapply(x, function(x) {
                                  val = x[["model"]][["WAAS"]]
                                })
      )
      )
    }
    if (index == "WAASBY"){
      bind = data.frame(do.call(cbind,
                                lapply(x, function(x) {
                                  val = x[["model"]][["WAASY"]]
                                })
      )
      )
    }
    bind$gen = x[[1]][["model"]][["Code"]]
    bind$type = x[[1]][["model"]][["type"]]
    data = data.frame(subset(bind, type == "GEN") %>% select(- type) %>% select(gen, everything()))
}


  if(class(x) == "WAASB"){
    if (index == "WAASB"){
      bind = data.frame(do.call(cbind,
                                lapply(x, function(x) {
                                  val = x[["model"]][["WAASB"]]
                                })
      )
      )
    }
    if (index == "WAASBY"){
      bind = data.frame(do.call(cbind,
                                lapply(x, function(x) {
                                  val = x[["model"]][["WAASBY"]]
                                })
      )
      )
    }
    bind$gen = x[[1]][["model"]][["Code"]]
    bind$type = x[[1]][["model"]][["type"]]
    data = data.frame(subset(bind, type == "GEN") %>% select(-type) %>% select(gen, everything()))
  }

if (any(apply(data, 2, is.na)) == TRUE) {
  stop("The MTSI index cannot be computed because the genotype-vs-environment interaction effect was not significant for some variable.")
}

if(is.null(SI)){
  ngs = NULL
} else {
  ngs  = round(nrow(data)*(SI/100),0)
}
  means <- data[,2:ncol(data)]
  rownames(means) <- data[,1]
  cor.means <- cor(means)
  eigen.decomposition <- eigen(cor.means)
  eigen.values <- eigen.decomposition$values
  eigen.vectors <- eigen.decomposition$vectors
  colnames(eigen.vectors) <- paste("PC",1:ncol(cor.means),sep="")
  rownames(eigen.vectors) <- colnames(means)
  if(length(eigen.values[eigen.values >= mineval]) == 1){
    eigen.values.factors <- as.vector(c(as.matrix(sqrt(eigen.values[eigen.values >= mineval]))))
    initial.loadings <- cbind(eigen.vectors[, eigen.values >= mineval]*eigen.values.factors)
    A <- initial.loadings
  } else {
    eigen.values.factors <- t(replicate(ncol(cor.means), c(as.matrix(sqrt(eigen.values[eigen.values >= mineval])))))
    initial.loadings <- eigen.vectors[, eigen.values >= mineval]*eigen.values.factors
    A <- varimax(initial.loadings)[[1]][]
  }
  partial = solve(cor.means)
  k = ncol(means)
  seq_k <- seq_len(ncol(means))
  for (j in seq_k) {
    for (i in seq_k) {
      if (i == j) {
        next
      }
      else {
        partial[i, j] <- -partial[i, j]/sqrt(partial[i, i] * partial[j, j])
      }
    }
  }
  KMO <- sum((cor.means[!diag(k)])^2)/(sum((cor.means[!diag(k)])^2) +
                                      sum((partial[!diag(k)])^2))
  MSA <- unlist(lapply(seq_k, function(i) {
    sum((cor.means[i, -i])^2)/(sum((cor.means[i, -i])^2) + sum((partial[i, -i])^2))
  }))
  names(MSA) <- colnames(means)

  colnames(A) <- paste("FA",1:ncol(initial.loadings),sep="")
  variance = (eigen.values/sum(eigen.values))*100
  cumulative.var <- cumsum(eigen.values/sum(eigen.values))*100
  pca <- cbind(eigen.values, variance, cumulative.var)
  colnames(pca) = c("Eigenvalues", "Variance (%)", "Cum. variance (%)")
  rownames(pca) <- paste("PC",1:ncol(means),sep="")
  Communality = diag(A %*% t(A))
  Uniquenesses = 1 - Communality
  fa <- cbind(A, Communality, Uniquenesses)
  z = scale(means, center = F, scale = apply(means, 2, sd))
  canonical.loadings = t(t(A) %*% solve(cor.means))
  scores <- z %*% canonical.loadings
  colnames(scores) <- paste("FA",1:ncol(scores),sep="")
  rownames(scores) <- data[,1]
  pos.var.factor <- which(abs(A) == apply(abs(A),1,max) , arr.ind = T)
  var.factor <- lapply(1:ncol(A),function(i){rownames(pos.var.factor)[pos.var.factor[,2] == i]})
  names(var.factor) <- paste("FA",1:ncol(A),sep="")
  names.pos.var.factor <- rownames(pos.var.factor)
  if (index == "WAASB"){
    ideotype.D <- apply(means, 2, min)
  } else {
    names(ideotype.D) <- colnames(means)
  }
  ideotypes.matrix <- t(as.matrix(ideotype.D))/apply(means, 2, sd)
  rownames(ideotypes.matrix) <- "ID1"
  ideotypes.scores <- ideotypes.matrix %*% canonical.loadings
  gen_ide = sweep(scores, 2, ideotypes.scores, "-")
  MTSI = sort(apply(gen_ide, 1, function(x) sqrt(sum(x^2))), decreasing = F)
  contr.factor = (gen_ide^2 / apply(gen_ide, 1, function(x) sum(x^2)))* 100
  means.factor <- means[,names.pos.var.factor]
  if (!is.null(ngs)){
    selection.diferential <- data.frame(cbind(Factor = pos.var.factor[,2],
                                  Xo = colMeans(means.factor),
                                  Xs = colMeans(means.factor[names(MTSI)[1:ngs],]),
                                  SD = colMeans(means.factor[names(MTSI)[1:ngs],]) - colMeans(means.factor) ,
                                  SDperc = (colMeans(means.factor[names(MTSI)[1:ngs],]) - colMeans(means.factor)) /
                                           colMeans(means.factor) * 100))
    selection.diferential[,1] = paste("FA", selection.diferential[,1], sep = "")
    sd_mean = apply(selection.diferential[,2:5], 2, mean)
  }
  if (is.null(ngs)){
    selection.diferential <- NULL
  }
  if(show){
    cat("\n-------------------------------------------------------------------------------\n")
    cat("Principal Component Analysis\n")
    cat("-------------------------------------------------------------------------------\n")
    print(pca)
    cat("\n-------------------------------------------------------------------------------\n")
    cat("Factor Analysis - factorial loadings after rotation-\n")
    cat("-------------------------------------------------------------------------------\n")
    print(fa)
    cat("\n-------------------------------------------------------------------------------\n")
    cat("Comunalit Mean:",mean(Communality),"\n")
    cat("\n-------------------------------------------------------------------------------\n")
    cat("Multitrait stability index\n")
    cat("-------------------------------------------------------------------------------\n")
    print(round(MTSI, 4))
    cat("\n-------------------------------------------------------------------------------\n")
    if(!is.null(ngs)){
      cat("Selection differential\n")
      cat("-------------------------------------------------------------------------------\n")
      print(selection.diferential)
      cat("\n------------------------------------------------------------------------------\n")
      cat("Mean of selection differential\n")
      cat("-------------------------------------------------------------------------------\n")
      print(sd_mean)
      cat("\n------------------------------------------------------------------------------\n")
      cat("Selected genotypes\n")
      cat(names(MTSI)[1:ngs])
      cat("\n-------------------------------------------------------------------------------\n")
    }
  }

  return(structure(list(data = data,
                        cormat = as.matrix(cor.means),
                        PCA = data.frame(pca),
                        FA = data.frame(fa),
                        KMO = KMO,
                        MSA = MSA,
                        comunalits = Communality,
                        comunalits.mean = mean(Communality),
                        initial.loadings = initial.loadings,
                        finish.loadings = A,
                        canonical.loadings = canonical.loadings,
                        scores.gen = scores,
                        scores.ide = ideotypes.scores,
                        MTSI = MTSI,
                        contri.fac = data.frame(contr.factor),
                        selection.diferential = selection.diferential,
                        selec.dif.mean = sd_mean,
                        Selected = names(MTSI)[1:ngs]),
                   class = "MTSI"))
}
