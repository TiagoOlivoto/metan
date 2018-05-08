\name{summary.WAASB}
\alias{summary.WAASB}
\title{Summary a WAASB object}
\usage{
\method{summary}{WAASB}(object, export = FALSE, file.name = "WAASB", ...)
}

\description{
Summary the  \code{WAASB} object in two ways. By default, the results are shown in the R console. The results can also be exported to the directory.

}

\arguments{

  \item{object}{The \code{WAASB} object}

  \item{export}{A logical argument. If \code{TRUE|T}, a *.txt file is exported to the working directory}

  \item{file.name}{The name of the file if \code{export = TRUE} }

  \item{...}{Other parameters of the plot.}

}

\author{Tiago Olivoto \email{tiagoolivoto@gmail.com}}


\examples{
library(WAASB)
model = WAASB(WAASBdata, resp = "GY", prob = 0.95, weight.response = 50, weight.WAAS = 50)
summary.WAASB(model, export = T, file.name = "my results")
}


