\name{plot.scores}
\alias{plot.scores}
\title{Plot scores in different graphical interpretations}

\description{

  Plot scores of genotypes and environments in different graphics. \code{1 = PC1 x PC2}, \code{2 =  GY x PC1},  \code{3 = GY x WAASB}, and \code{4 = Nominal yield x EPCA1}.

}

\usage{
\method{plot}{scores}(x, type = 1, polygon = FALSE, file.type = "pdf",
   export = FALSE, file.name = NULL, theme = theme_waasb(),
   axis.expand = 1.1,  width = 8, height = 7, x.lim = NULL,
   x.breaks = waiver(), x.lab = NULL, y.lab = NULL,
   y.lim = NULL, y.breaks = waiver(), shape.gen = 21,
   shape.env = 23, size.shape = 2.5, size.tex.lab = 12,
   size.tex.pa = 3.5, size.line = 0.5, size.segm.line = 0.5,
   leg.lab = c("Gen", "Env"), line.type = "solid", line.alpha = 0.9,
   col.line = "gray", col.gen = "blue", col.env = "darkgreen",
   col.alpha = 0.9, col.segm.gen = "transparent",
   col.segm.env = "darkgreen", resolution = 300, ...)

}

\arguments{

  \item{x}{The object \code{WAASB} or \code{WAAS.AMMI}}

  \item{type}{Three types of graphics can be generated: \code{1 = PC1 x PC2}, default, to make inferences related to the interaction effects; \code{2 = GY x PC1} to make inferences related to stability and productivity; \code{3 = GY x WAASB}, and \code{4 = Nominal yield x Environment PC1}.}

  \item{polygon}{Logical argument. If \code{TRUE}, a polygon is drawn when \code{type 1}.}

  \item{file.type}{The type of file to be exported. Valid parameter if \code{export = T|TRUE}.  Default is \code{"pdf"}. The graphic can also be exported in \code{*.tiff} format by declaring  \code{file.type = "tiff"}.}

  \item{export}{Export (or not) the plot. Default is \code{FALSE}.}

  \item{file.name}{The name of the file for exportation, default is \code{NULL}, i.e. the files are automatically named.}

  \item{theme}{The graphical theme of the plot. Default is `theme = theme_waasb()`. Please, see `?WAASB::theme_waasb`. An own theme can be applied using the arguments: `theme = theme_waasb() + theme(some stuff here)`. For more details, please, see `?ggplot2::theme`}

  \item{axis.expand}{Multiplication factor to expand the axis limits by to enable fitting of labels. Default is \code{1.1}.}

  \item{width}{The width "inch" of the plot. Default is \code{8}.}

  \item{height}{The height "inch" of the plot. Default is \code{7}.}

  \item{x.lim}{The range of x-axis. Default is \code{NULL} (maximum and minimum values of the data set). New arguments can be inserted as \code{x.lim = c(x.min, x.max)}.}

  \item{x.breaks}{The breaks to be plotted in the x-axis. Default is \code{authomatic breaks}. New arguments can be inserted as \code{x.breaks = c(breaks)}}

  \item{x.lab}{The label of x-axis. Each plot has a default value. New arguments can be inserted as \code{x.lab = "my label"}.}

  \item{y.lab}{The label of y-axis. Each plot has a default value. New arguments can be inserted as \code{y.lab = "my label"}.}

  \item{y.lim}{The range of x-axis. Default is \code{NULL}. The same arguments than \code{x.lim} can be used.}

  \item{y.breaks}{The breaks to be plotted in the x-axis. Default is \code{authomatic breaks}. The same arguments than \code{x.breaks} can be used.}

  \item{shape.gen}{The shape for genotype indication in the biplot. Default is \code{21} (circle). Values must be between \code{21-25}: \code{21} (circle), \code{22} (square), \code{23} (diamond), \code{24} (up triangle), and \code{25} (low triangle).}

  \item{shape.env}{The shape for environment indication in the biplot. Default is \code{23} (diamond). The same arguments than \code{"shape.gen"}.}

  \item{size.shape}{The size of the shape (both for genotypes and environments). Default is \code{2.5}.}

  \item{size.tex.lab}{The size of the text in the axes text and labels. Default is \code{12}.}

  \item{size.tex.pa}{The size of the text of the plot area. Default is \code{3.5}.}

  \item{size.line}{The size of the line that indicate the means in the biplot. Default is \code{0.5}.}

  \item{size.segm.line}{The size of the segment that start in the origin of the biplot and end in the scores values. Default is \code{0.5}.}


  \item{leg.lab}{The labs of legend. Default is \code{Gen} and \code{Env}. }

  \item{line.type}{The type of the line that indicate the means in the biplot. Default is \code{"solid"}. Other values that can be attributed are: \code{"blank"}, no lines in the biplot, \code{"dashed", "dotted", "dotdash", "longdash", and "twodash"}.}

  \item{line.alpha}{The alpha value that combine the line with the background to create the appearance of partial or full transparency. Default is \code{0.4}. Values must be between "0" (full transparency) to "1" (full color).}

  \item{col.line}{The color of the line that indicate the means in the biplot. Default is \code{"gray"}}

  \item{col.gen}{The shape color for genotypes. Must be one value or a vector of colours with the same length of the number of genotypes. Default is \code{"blue"}. Other values can be attributed. For example, \code{"transparent"}, will make a plot with only an outline around the shape area.}

  \item{col.env}{The shape color for environments. Default is \code{"red"}. The same usability than \code{"col.gen"}.}

  \item{col.alpha}{The alpha value for the color. Default is \code{0.9}. Values must be between \code{0} (full transparency) to \code{1} (full color).}

  \item{col.segm.gen}{The color of segment for genotypes.Default is \code{"transparent"}. Parameter valid for \code{type = 1} and \code{type = 2} graphics. This segment start in the origin of the biplot and end in the scores values.}

  \item{col.segm.env}{The color of segment for environments. Default is \code{"gray50"}. The same usability than \code{"col.segm.gen"} }

  \item{resolution}{The resolution of the plot. Parameter valid if \code{file.type = "tiff"} is used. Default is \code{300} (300 dpi)}

  \item{...}{Other arguments of the function}

}

\details{
The plots type 1 and 2 have the same interpretation than those used in traditional-usage AMMI analysis (well know as AMMI2 and AMMI1, respectively). In the plot type 3, the scores of both genotypes and environments are plotted considering the response variable and the WAASB (stability index that considers all significant principal component axis of traditional AMMI models or all principal component axis estimated with BLUP-interaction effects. Plot type 4 may be used to better understand the well known "which-won-where" pattern, facilitating the recommendation of appropriate genotypes targeted for specific environments, thus allowing the exploitation of narrow adaptations.
}

\examples{
library(METAAB)
library(ggplot2)
scores = WAASB(data_ge,
               resp = GY,
               gen = GEN,
               env = ENV,
               rep = REP)
# PC1 x PC2
plot.scores(model$GY,
            type = 1,
            polygon = TRUE)

# GY x PC1
plot.scores(model$GY,
            type = 2,
            col.env = "olivedrab",
            col.gen = "orange2",
            x.lab = "My own x label")

# GY x WAASB
plot.scores(model$GY,
            type = 3,
            size.tex.pa = 2,
            size.tex.lab = 16)

}

\author{Tiago Olivoto \email{tiagoolivoto@gmail.com}}

\seealso{\code{\link{plot.eigen}}, \code{\link{plot.validation.AMMIF}}}


