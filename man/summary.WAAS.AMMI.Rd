\name{summary.WAAS.AMMI}
\alias{summary.WAAS.AMMI}
\title{Summary a WAAS.AMMI object}
\usage{
\method{summary}{WAAS.AMMI}(object, export = FALSE, file.name = NULL, digits = 4, ...)
}

\description{
Summary the  \code{WAAS.AMMI} object in two ways. By default, the results are shown in the R console. The results can also be exported to the directory.

}

\arguments{

  \item{object}{The \code{WAAS.AMMI} object}

  \item{export}{A logical argument. If \code{TRUE|T}, a *.txt file is exported to the working directory}

  \item{file.name}{The name of the file if \code{export = TRUE} }

  \item{digits}{The significant digits to be shown.}

  \item{...}{Other arguments of the function}

}

\author{Tiago Olivoto \email{tiagoolivoto@gmail.com}}


\examples{
library(METAAB)
model = WAAS.AMMI(data_ge,
                  resp = c(GY, HM),
                  gen = GEN,
                  env = ENV,
                  rep = REP)
summary(model)
summary(model,
        export = TRUE,
        file.name = "my results")
}


