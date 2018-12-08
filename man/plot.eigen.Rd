\name{plot.eigen}
\alias{plot.eigen}
\title{Plot the eigenvalues }
\usage{
\method{plot}{eigen}(x, export = FALSE, theme = theme_waasb(),
     file.type = "pdf", file.name = NULL, width = 6, height = 6,
     size.shape = 3.5, size.line = 1, y.lab = "Eigenvalue",
     y2.lab = "Accumulated variance", x.lab = "Number of multiplicative terms",
     resolution = 300, ...)
}

\description{
Plot the eigenvalues for from singular value decomposition of BLUP interaction effects matrix.


}
\arguments{

  \item{x}{The \code{WAASB object} }

  \item{export}{Export (or not) the plot. Default is \code{TRUE}.}

  \item{theme}{The graphical theme of the plot. Default is `theme = theme_waasb()`. Please, see `?WAASB::theme_waasb`. An own theme can be applied using the arguments: `theme = theme_waasb() + theme(some stuff here)`. For more details, please, see `?ggplot2::theme`}

  \item{file.type}{If \code{export = TRUE}, define the type of file to be exported. Default is \code{pdf}, Graphic can also be exported in \code{*.tiff} format by declaring  \code{file.type = "tiff"}.}

  \item{file.name}{The name of the file for exportation, default is \code{NULL}, i.e. the files are automatically named.}

  \item{width}{The width "inch" of the plot. Default is \code{6}.}

  \item{height}{The height "inch" of the plot. Default is \code{6}. }

  \item{size.shape}{The size of the shape. Default is \code{3.5}.}

  \item{size.line}{The size of the line. Default is \code{1}.}

  \item{y.lab}{The label of the y-axis in the plot. Default is \code{"Eigenvalue"}.}

  \item{y2.lab}{The label of the second y-axis in the plot. Default is \code{"Accumulated variance"}.}

  \item{x.lab}{The label of the x-axis in the plot. Default is \code{"Number of multiplicative terms"}.}

  \item{resolution}{The resolution of the plot. Parameter valid if \code{file.type = "tiff"} is used. Default is \code{300} (300 dpi)}

  \item{...}{Other arguments of the function}

}
\examples{
library(WAASB)
model = WAASB(WAASBdata, GY, GEN, ENV, REP)
plot.eigen(model)


}

\author{Tiago Olivoto \email{tiagoolivoto@gmail.com}}


\seealso{\code{\link{plot.scores}}, \code{\link{plot.WAASBY}}}
