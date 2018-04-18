\name{plot.WAASBY}
\alias{plot.WAASBY}
\title{Plot WAASBY values for genotype ranking}
\usage{
\method{plot}{WAASBY}(data, export , file.type,  width , height , size.lab,
            size.shape , leg.pos, col.shape, y.lab, x.breaks, resolution)
}

\description{
Plot heat maps with genotype ranking in two ways.

}
\arguments{

  \item{data}{The \code{WAASBYratio object} }

  \item{export}{Export (or not) the plot. Default is \code{T}.}

  \item{file.type}{The type of file to be exported. Default is \code{pdf}, Graphic can also be exported in \code{*.tiff} format by declaring  \code{file.type = "tiff"}. Editable vector-based graphics can be saved into a *.pptx file using \code{file.type = "pptx"}}

  \item{width}{The width "inch" of the plot. Default is \code{8}.}

  \item{height}{The height "inch" of the plot. Default is \code{7}. }

  \item{size.lab}{The size of the labels in the plot. Default is \code{12}.}

  \item{size.shape}{The size of the shape in the plot. Default is \code{3.5}.}

  \item{leg.pos}{A vector of length 2 to configure, respectively, the horizontal and vertical,  position of the legend in the plot. Default is \code{c(0.9, 0.1)} (bottom-right. }

  \item{col.shape}{A vector of length 2 that contains the color of shapes for genotypes above and below of the mean, respectively. Default is \code{c("blue", "red")}. }

  \item{y.lab}{The label of the y-xix in the plot. Default is \code{"Genotypes"}.}

  \item{x.breaks}{The breaks to be plotted in the x-axis. Default is \code{authomatic breaks}. New arguments can be inserted as \code{x.breaks = c(breaks)}}

  \item{resolution}{The resolution of the plot. Parameter valid if \code{file.type = "tiff"} is used. Default is \code{300} (300 dpi)}
  }

\details{
The first type of heatmap shows the genotype ranking depending on the number of principal component axis used for estimating the WAASB index. An euclidian distance-based dendrogram is used for grouping the genotype ranking for both genotypes and principal component axis. The second type of heatmap shows the genotype ranking depending on the WAASB/GY ratio. The ranks obtained with a ratio of 100/0 considers exclusively the stability for the genotype ranking. On the other hand, a ratio of 0/100 considers exclusively the productivity for the genotype ranking.  Four clusters are estimated (1) unproductive and unstable genotypes; (2) productive, but unstable genotypes; (3) stable, but unproductive genotypes; and (4), productive and stable genotypes.

}


\examples{
library(WAASB)
# plotting the WAASBY values based on blup procedure
model = WAASBYratio(WAASBdata, resp = "GY", increment = 5, saveWAASY = 60)
plot.WAASBY(model)

# plotting the WAASY value based on traditional AMMI model
model2 = WAASratio.AMMI(WAASBdata, resp = "GY", increment = 5, saveWAASY = 60)
plot.WAASBY(model2)

}

\author{Tiago Olivoto \email{tiagoolivoto@gmail.com}}


\seealso{\code{\link{plot.scores}}}
