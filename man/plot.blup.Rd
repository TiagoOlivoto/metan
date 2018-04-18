\name{plot.blup}
\alias{plot.blup}
\title{Plot the estimated BLUPs of genotypes }
\usage{
\method{plot}{blup}(data, prob, export, width, height, \dots)
}

\description{
Plot the predicted BLUP of the genotypes.


}
\arguments{

  \item{data}{The \code{WAASB object} }

  \item{prob}{The probability to be considered in confidence interval prediction. }

  \item{export}{Export (or not) the plot. Default is \code{TRUE}.}

  \item{file.type}{If \code{export = TRUE}, define the type of file to be exported. Default is \code{pdf}, Graphic can also be exported in \code{*.tiff} format by declaring  \code{file.type = "tiff"}. Editable vector-based graphics can be saved into a *.pptx file using \code{file.type = "pptx"}}

  \item{width}{The width "inch" of the plot. Default is \code{6}.}

  \item{height}{The height "inch" of the plot. Default is \code{6}. }

  \item{size.lab}{The size of the labels in the plot. Default is \code{12}.}

  \item{size.tex}{The size of the tick labels in the plot. Default is \code{12}.}

  \item{size.leg}{The size of the legend text. Default is \code{12}.}

  \item{size.err.bar}{The size of the error bar for the plot. Default is \code{0.5}.}

  \item{size.shape}{The size of the shape (both for genotypes and environments). Default is \code{3.5}.}

  \item{height.err.bar}{The height for error bar. Default is \code{0.3}.}

   \item{x.lim}{The range of x-axis. Default is \code{NULL} (maximum and minimum values of the data set). New arguments can be inserted as \code{x.lim = c(x.min, x.max)}.}

  \item{x.breaks}{The breaks to be plotted in the x-axis. Default is \code{authomatic breaks}. New arguments can be inserted as \code{x.breaks = c(breaks)}}



  \item{leg.pos}{A vector of length 2 to configure, respectively, the horizontal and vertical,  position of the legend in the plot. Default is \code{c(0.9, 0.1)} (bottom-right. }

  \item{col.shape}{A vector of length 2 that contains the color of shapes for genotypes above and below of the mean, respectively. Default is \code{c("blue", "red")}. }

  \item{x.lab}{The label of the x-axis in the plot. Default is \code{"Predicted Grain Yield"}.}

  \item{y.lab}{The label of the y-axis in the plot. Default is \code{"Genotypes"}.}


  \item{resolution}{The resolution of the plot. Parameter valid if \code{file.type = "tiff"} is used. Default is \code{300} (300 dpi)}

  \item{...}{Additional parameter for the function}



}
\examples{
library(WAASB)
model = WAASB(WAASBdata, resp = "GY", weight.response = 50, weight.WAAS = 50)
plot.blup(model)


}

\author{Tiago Olivoto \email{tiagoolivoto@gmail.com}}


\seealso{\code{\link{plot.scores}}, \code{\link{plot.WAASBY}}}
