\name{summary.WAASB}
\alias{summary.WAASB}
\title{Summary a WAASB object}
\usage{
\method{summary}{WAASB}(data)
}

\description{
Summary the  \code{WAAS} object in two ways. By default, the results are shown in the R console. The results can also be exported to the directory.

}

\arguments{

  \item{data}{The \code{WAAS object} }

}

\author{Tiago Olivoto \email{tiagoolivoto@gmail.com}}


\seealso{\code{\link{write.WAASB}}}

\examples{
library(WAASB)
model = WAASB(WAASBdata, prob = 0.95, weight.response = 50, weight.WAAS = 50)
summary.WAASB(model)
}


