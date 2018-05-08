\name{plot.scores}
\alias{plot.scores}
\title{Plot scores in different graphical interpretations}

\description{

  Plot scores of genotypes and environments in different graphics. \code{1 = PC1 x PC2}, \code{2 =  GY x PC1}, and \code{3 = GY x WAASB}.)

}

\usage{
\method{plot}{WAASB, WAAS.AMMI}(x, type, file.type = "pdf", export = FALSE, width = 8, height = 7, x.lim = NULL, x.breaks = waiver(), x.lab = "Grain yield", y.lim = NULL, y.breaks = waiver(), shape.gen = 21, shape.env = 23, size.shape = 3.5, size.lab = 18, size.tex = 3.5, size.line = 0.5, size.leg = 14, size.segm.line = 0.5, leg.position = "tr", leg.lab = c("Gen", "Env"), line.type = "dashed", line.alpha = 0.9, col.line = "gray", col.gen = "gray75", col.env = "red", col.alpha = 0.9, col.segm.gen = "transparent", col.segm.env = "grey50", resolution = 300, \dots)

}

\arguments{

  \item{x}{The object \code{WAASB} or \code{WAAS.AMMI}}

  \item{type}{Three types of graphics can be generated: \code{1 = PC1 x PC2}, to make inferences related to the interaction effects; \code{2 = GY x PC1} to make inferences related to stability and productivity; and \code{3 = GY x WAASB}. In this plot, the scores of both genotypes and environements are plotted considering the response variable and the WAASB (stability index that consider all significant principal component axis of traditional AMMI models or all principal component axis estimated with BLUP-interaction effects.}

    \item{export}{Export (or not) the plot. Default is \code{TRUE}.}

  \item{file.type}{The type of file to be exported. Valid parameter if \code{export = T|TRUE}.  Default is \code{"pdf"}. The graphic can also be exported in \code{*.tiff} format by declaring  \code{file.type = "tiff"}. Editable vector-based graphics can be saved into a *.pptx file using \code{file.type = "pptx"}}


  \item{width}{The width "inch" of the plot. Default is \code{8}.}

  \item{height}{The height "inch" of the plot. Default is \code{7}.}

  \item{x.lim}{The range of x-axis. Default is \code{NULL} (maximum and minimum values of the data set). New arguments can be inserted as \code{x.lim = c(x.min, x.max)}.}

  \item{x.breaks}{The breaks to be plotted in the x-axis. Default is \code{authomatic breaks}. New arguments can be inserted as \code{x.breaks = c(breaks)}}

  \item{x.lab}{The label of x-axis for the plot \code{type = 2 and 3} . Default is \code{Grain yield}. New arguments can be inserted as \code{x.lab = c(my label)}}

  \item{y.lim}{The range of x-axis. Default is \code{NULL}. The same arguments than \code{x.lim} can be used.}

  \item{y.breaks}{The breaks to be plotted in the x-axis. Default is \code{authomatic breaks}. The same arguments than \code{x.breaks} can be used.}

  \item{shape.gen}{The shape for genotype indication in the biplot. Default is \code{21} (circle). Values must be between \code{21-25}: \code{21} (circle), \code{22} (square), \code{23} (diamond), \code{24} (up triangle), and \code{25} (low triangle).}

  \item{shape.env}{The shape for environment indication in the biplot. Default is \code{23} (diamond). The same arguments than \code{"shape.gen"}.}

  \item{sizes.hape}{The size of the shape (both for genotypes and environments). Default is \code{3.5}.}

  \item{size.lab}{The size of the plot label. Default is \code{18}.}

  \item{size.tex}{The size of the text in the background of the plot. Default is \code{3.5}.}

  \item{size.line}{The size of the line that indicate the means in the biplot. Default is \code{0.5}.}

  \item{size.leg}{The size of the legend text. Default is \code{16}.}

  \item{size.segm.line}{The size of the segment that start in the origin of the biplot and end in the scores values. Default is \code{0.5}.}


  \item{leg.position}{The position of the legend in the graphic. Default is \code{tr} (top right). Possible values: \code{tr} (top right), \code{tl} (top left), \code{br} (bottom right) and \code{bl} (bottom left). }

  \item{leg.lab}{The labs of legend. Default is \code{Gen} and \code{Env}. }

  \item{line.type}{The type of the line that indicate the means in the biplot. Default is \code{"solid"}. Other values that can be attributed are: \code{"blank"}, no lines in the biplot, \code{"dashed", "dotted", "dotdash", "longdash", and "twodash"}.}

  \item{line.alpha}{The alpha value that combine the line with the background to create the appearance of partial or full transparency. Default is \code{0.4}. Values must be between "0" (full transparency) to "1" (full color).}

  \item{col.line}{The color of the line that indicate the means in the biplot. Default is \code{"gray"}}

  \item{col.gen}{The shape color for genotypes. Default is \code{"gray75"}. Other values can be attributed. For example, \code{"transparent"}, will make a plot with only an outline around the shape area.}

  \item{col.env}{The shape color for environments. Default is \code{"red"}. The same usability than \code{"col.gen"}.}

  \item{col.alpha}{The alpha value for the color. Default is \code{0.9}. Values must be between \code{0} (full transparency) to \code{1} (full color).}

  \item{col.segm.gen}{The color of segment for genotypes.Default is \code{"transparent"}. Parameter valid for \code{type = 1} and \code{type = 2} graphics. This segment start in the origin of the biplot and end in the scores values.}

  \item{col.segm.env}{The color of segment for environments. Default is \code{"gray50"}. The same usability than \code{"col.segm.gen"} }

  \item{resolution}{The resolution of the plot. Parameter valid if \code{file.type = "tiff"} is used. Default is \code{300} (300 dpi)}

  \item{...}{Additional parameter for the function}

}

\examples{
  library(WAASB)
  model = WAASB(WAASBdata, resp = "GY", weight.response = 50, weight.WAAS = 50)
  plot.scores(model, 1, file.type = "tiff", width = 8, height = 7, col.gen = "green")
  # PC1 x PC2.tiff

  plot.scores(model, 2, width = 8, height = 7, col.gen = "green", shape.gen = 23)
  # GY x PC1.pdf

  plot.scores(model, 3, col.gen = "blue", col.env = "red", size.leg = 17)
  # GY x WAAS.pdf
}

\author{Tiago Olivoto \email{tiagoolivoto@gmail.com}}

\seealso{\code{\link{plot.WAASBY}}}


