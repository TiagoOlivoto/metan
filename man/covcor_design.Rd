\name{covcor_design}
\alias{covcor_design}
\title{Variance-covariance and correlation matrices for designed experiments}
\usage{
covcor_design(.data, gen, rep, resp, design = "RCBD", type = NULL)
}
\description{
Compute variance-covariance and correlation matrices using data from a designed (RCBD or CRD).

}

\arguments{
 \item{.data}{The dataset containing the columns related to Genotypes, replication/block and response variables. Alternatively, it is possible to use an object of class 'group_factors' to compute the results for each level of the grouping factor. See \code{?group_factors}.}

  \item{gen}{The name of the column that contains the levels of the genotypes.}

  \item{rep}{The name of the column that contains the levels of the replications/blocks.}

  \item{resp}{The response variables. For example \code{resp = c(var1, var2, var3)}.}

  \item{design}{The experimental design. Must be RCBD or CRD.}

  \item{type}{What the matrices should return? Set to \code{NULL}, i.e., a list of matrices is returned. The argument type allow the following values 'pcor', 'gcor', 'rcor', (which will return the phenotypic, genotypic and residual correlation matrices, respectively) or 'pcov', 'gcov', 'rcov' (which will return the phenotypic, genotypic and residual variance-covariance matrices, respectively).}

}


\author{Tiago Olivoto \email{tiagoolivoto@gmail.com}}




