\name{write.WAASB}
\alias{write.WAASB}
\title{Write the output of a WAASB function}
\usage{
write.WAASB(data, file.type = "xlsx", file.name = "Results WAAS")
}

\description{
Write \code{*.xlsx} and \code{*.txt} files with the outputs of \code{WAASB} function

}
\arguments{

  \item{data}{The \code{WAAS object} }

  \item{file.type}{The type of file to be exported. Default is \code{xlsx}, but *txt files can be created by declaring \code{file.type = "txt"}.}

  \item{file.name}{The name of file to be exported. Default is \code{Results WAAS}.}
}

\examples{
library(WAASB)
model = WAASB(WAASBdata, prob = 0.95, weight.response = 50, weight.WAAS = 50)
write.WAASB(model$PCA, file.name = "results PC")
}

\author{Tiago Olivoto \email{tiagoolivoto@gmail.com}}


\seealso{\code{\link{plot.scores}}}
