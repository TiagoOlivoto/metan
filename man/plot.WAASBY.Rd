\name{plot.WAASBY}
\alias{plot.WAASBY}
\title{Plot WAASBY values for genotype ranking}

\description{
Plot heat maps with genotype ranking in two ways.

}

\usage{
\method{plot}{WAASBY}(x, export = F , file.type = "pdf", file.name = NULL,
    theme = theme_waasb(), width = 6, height = 6, size.shape = 3.5,
    col.shape = c("blue", "red"), y.lab = "Genotypes", x.breaks = waiver(),
    resolution = 300, ...)}


\arguments{

  \item{x}{The \code{WAASBY object} }

  \item{export}{Export (or not) the plot. Default is \code{T}.}

  \item{file.type}{The type of file to be exported. Default is \code{pdf}, Graphic can also be exported in \code{*.tiff} format by declaring  \code{file.type = "tiff"}.}

  \item{file.name}{The name of the file for exportation, default is \code{NULL}, i.e. the files are automatically named.}

  \item{theme}{The graphical theme of the plot. Default is `theme = theme_waasb()`. Please, see `?WAASB::theme_waasb`. An own theme can be applied using the arguments: `theme = theme_waasb() + theme(some stuff here)`. For more details, please, see ` ?ggplot2::theme`}

  \item{width}{The width "inch" of the plot. Default is \code{8}.}

  \item{height}{The height "inch" of the plot. Default is \code{7}. }

  \item{size.shape}{The size of the shape in the plot. Default is \code{3.5}.}

  \item{col.shape}{A vector of length 2 that contains the color of shapes for genotypes above and below of the mean, respectively. Default is \code{c("blue", "red")}. }

  \item{y.lab}{The label of the y-xix in the plot. Default is \code{"Genotypes"}.}

  \item{x.breaks}{The breaks to be plotted in the x-axis. Default is \code{authomatic breaks}. New arguments can be inserted as \code{x.breaks = c(breaks)}}

  \item{resolution}{The resolution of the plot. Parameter valid if \code{file.type = "tiff"} is used. Default is \code{300} (300 dpi)}

  \item{...}{Other arguments of the function}

  }



\examples{
library(WAASB)
# plotting the WAASBY values based on blup procedure
#model = WAASBYratio(WAASBdata, resp = "GY", increment = 5, saveWAASY = 60)
#plot.WAASBY(model)

# plotting the WAASY value based on traditional AMMI model
#model2 = WAASratio.AMMI(WAASBdata, resp = "GY", increment = 5, saveWAASY = 60)
#plot.WAASBY(model2)

}

\author{Tiago Olivoto \email{tiagoolivoto@gmail.com}}


\seealso{\code{\link{plot.scores}}}
