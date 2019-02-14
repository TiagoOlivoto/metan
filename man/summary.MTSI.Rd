\name{summary.MTSI}
\alias{summary.MTSI}
\title{Summary a MTSI object}
\usage{
\method{summary}{MTSI}(object, export = FALSE, file.name = NULL, digits = 4, ...)
}

\description{
Summary a \code{MTSI} object in two ways. By default, the results are shown in the R console. The results can also be exported to the directory.

}

\arguments{

  \item{object}{The \code{MTSI} object}

  \item{export}{A logical argument. If \code{TRUE|T}, a *.txt file is exported to the working directory}

  \item{file.name}{The name of the file if \code{export = TRUE} }

  \item{digits}{The significant digits to be shown.}

  \item{...}{Other arguments of the function}

}

\author{Tiago Olivoto \email{tiagoolivoto@gmail.com}}


\examples{
library(METAAB)
# Based on stability only
MTSI_MODEL = WAASB(data_ge,
                   resp = c(GY, HM),
                   gen = GEN,
                   env = ENV,
                   rep = REP)

MTSI_index = MTSI(MTSI_MODEL)
summary(MTSI_index)
summary(MTSI_index,
        export = TRUE,
        file.name = "my results")
}


