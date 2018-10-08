\name{plot.validation.AMMIF}
\alias{plot.validation.AMMIF}
\title{Plot the RMSPD of all AMMI-family tested models}
\usage{

\method{plot}{validation.AMMIF}(x, violin = TRUE, export = FALSE,  x.lab = NULL,
     y.lab = NULL, file.type = "pdf", file.name = NULL,
     theme = theme_waasb(),width = 6, height = 6, resolution = 300,
     col.violin = "gray90", col.boxplot = "gray70", width.boxplot = 0.2,
     x.lim = NULL, x.breaks = waiver(), ...)
}

\description{
Boxplot showing the Root Means Square Prediction Difference of all AMMI family models obtained in a cross-validation procedure.

}

\arguments{

  \item{x}{An object of class \code{validation.AMMIF}.}

  \item{violin}{Define if a violine plot is used with boxplot. Default is "TRUE"}

    \item{export}{Export (or not) the plot. Default is \code{T}.}

   \item{x.lab}{The label of x-axis. New arguments can be inserted as \code{x.lab = "my x label"}.}

  \item{y.lab}{The label of y-axis. New arguments can be inserted as \code{y.lab = "my y label"}.}

  \item{file.type}{The type of file to be exported. Default is \code{pdf}, Graphic can also be exported in \code{*.tiff} format by declaring  \code{file.type = "tiff"}.}

  \item{file.name}{The name of the file for exportation, default is \code{NULL}, i.e. the files are automatically named.}

  \item{theme}{The graphical theme of the plot. Default is `theme = theme_waasb()`. Please, see `?WAASB::theme_waasb`. An own theme can be applied using the arguments: `theme = theme_waasb() + theme(some stuff here)`. For more details, please, see `?ggplot2::theme`}


  \item{width}{The width "inch" of the plot. Default is \code{6}.}

  \item{height}{The height "inch" of the plot. Default is \code{6}. }

    \item{resolution}{The resolution of the plot. Parameter valid if \code{file.type = "tiff"} is used. Default is \code{300} (300 dpi)}

    \item{col.violin}{Parameter valid if  \code{violin = T}. Define the color of the violin plot. Default is "gray90.}


   \item{col.boxplot}{Define the color for boxplot. Default is "gray70".}


   \item{width.boxplot}{The width of boxplots. Default is \code{0.2}.}


   \item{x.lim}{The range of x-axis. Default is \code{NULL} (maximum and minimum values of the data set). New arguments can be inserted as \code{x.lim = c(x.min, x.max)}.}

  \item{x.breaks}{The breaks to be plotted in the x-axis. Default is \code{authomatic breaks}. New arguments can be inserted as \code{x.breaks = c(breaks)}}

  \item{...}{Other arguments of the function}

}

\details{
Five statistics are shown in this type of plot. The lower and upper hinges correspond to the first and third quartiles (the 25th and 75th percentiles). The upper whisker extends from the hinge to the largest value no further than 1.5 * IQR from the hinge (where IQR is the inter-quartile range). The lower whisker extends from the hinge to the smallest value at most 1.5 * IQR of the hinge. Data beyond the end of the whiskers are considered outlying points.

}

\author{Tiago Olivoto \email{tiagoolivoto@gmail.com}}

\examples{
library(WAASB)
validation = validation.AMMIF(WAASBdata, GY, GEN, ENV, REP,
                              nboot = 100, nrepval = 2)

#A plot without a violin plot exported to a *.pdf file
plot(validation, violin = F, export = T)
}

\seealso{\code{\link{plot.scores}}, \code{\link{plot.WAASBY}}}
