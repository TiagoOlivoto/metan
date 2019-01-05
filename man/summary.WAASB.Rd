\name{summary.WAASB}
\alias{summary.WAASB}
\title{Summary a WAASB object}
\usage{
\method{summary}{WAASB}(object, export = FALSE, file.name = NULL, digits = 4, ...)
}

\description{
Summary the  \code{WAASB} object in two ways. By default, the results are shown in the R console. The results can also be exported to the directory.

}

\arguments{

  \item{object}{The \code{WAASB} object}

  \item{export}{A logical argument. If \code{TRUE|T}, a *.txt file is exported to the working directory}

  \item{file.name}{The name of the file if \code{export = TRUE} }

  \item{digits}{The significant digits to be shown.}

  \item{...}{Other arguments of the function}

}

\author{Tiago Olivoto \email{tiagoolivoto@gmail.com}}


\examples{
library(METAAB)
model = WAASB(WAASBdata, GY, GEN, ENV, REP)
summary(model, digits = 2)
summary(model, export = T, file.name = "my results")
}


