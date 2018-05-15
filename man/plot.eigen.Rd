\name{plot.eigen}
\alias{plot.eigen}
\title{Plot the eigenvalues }
\usage{
\method{plot}{eigen}(x, export = FALSE, file.type = "pdf", file.name = NULL, width = 6, height = 6, size.lab = 12, size.tex = 12, size.shape = 3.5, size.line = 1, col.shape = "blue", col.line = "blue", y.lab = "Eigenvalue", x.lab = "Number of multiplicative terms", resolution = 300, \dots)
}

\description{
Plot the eigenvalues for from singular value decomposition of BLUP interaction effects matrix.


}
\arguments{

  \item{x}{The \code{WAASB object} }

  \item{export}{Export (or not) the plot. Default is \code{TRUE}.}

  \item{file.type}{If \code{export = TRUE}, define the type of file to be exported. Default is \code{pdf}, Graphic can also be exported in \code{*.tiff} format by declaring  \code{file.type = "tiff"}.}

  \item{file.name}{The name of the file for exportation, default is \code{NULL}, i.e. the files are automatically named.}

  \item{width}{The width "inch" of the plot. Default is \code{6}.}

  \item{height}{The height "inch" of the plot. Default is \code{6}. }

  \item{size.lab}{The size of the labels in the plot. Default is \code{12}.}

  \item{size.tex}{The size of the tick labels in the plot. Default is \code{12}.}

  \item{size.shape}{The size of the shape. Default is \code{3.5}.}

  \item{size.shape}{The size of the line. Default is \code{1}.}


  \item{col.shape}{The color for shape. Default is \code{c("blue")}. }

  \item{col.line}{The color for line. Default is \code{c("blue")}. }

  \item{y.lab}{The label of the y-axis in the plot. Default is \code{"Eigenvalue"}.}

  \item{x.lab}{The label of the x-axis in the plot. Default is \code{"Number of multiplicative terms"}.}

  \item{resolution}{The resolution of the plot. Parameter valid if \code{file.type = "tiff"} is used. Default is \code{300} (300 dpi)}

  \item{...}{Additional parameter for the function}


}
\examples{
library(WAASB)
model = WAASB(WAASBdata, resp = "GY", weight.response = 50, weight.WAAS = 50)
plot.eigen(model)


}

\author{Tiago Olivoto \email{tiagoolivoto@gmail.com}}


\seealso{\code{\link{plot.scores}}, \code{\link{plot.WAASBY}}}
