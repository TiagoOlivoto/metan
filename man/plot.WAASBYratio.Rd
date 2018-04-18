\name{plot.WAASBYratio}
\alias{plot.WAASBYratio}
\title{Plot heat maps with genotype ranking}
\usage{
\method{plot}{WAASBYratio}(data, type, export = FALSE, file.type, width, height, labSize,
                           margins, y.lab, key.lab, resolution)
}

\description{
Plot heat maps with genotype ranking in two ways.

}
\arguments{

  \item{data}{The \code{WAASYratio object} }

  \item{type}{\code{1 = Heat map Ranks}: this graphic shows the genotype ranking considering the WAAS estimated with different numbers of Principal Components; \code{2 = Heat map WAASY-GY ratio}: this graphic shows the genotype ranking considering the different combinations in the WAAS/GY ratio. }

  \item{export}{Export (or not) the plot. Default is \code{TRUE}.}

  \item{file.type}{If \code{export = TRUE} define the type of file to be exported. Default is \code{pdf}, Graphic can also be exported in \code{*.tiff} format by declaring  \code{file.type = "tiff"}.}

  \item{width}{The width "inch" of the plot. Default is \code{8}.}

  \item{height}{The height "inch" of the plot. Default is \code{7}. }

  \item{labSize}{The label size of the plot. It is suggested attribute 1}

  \item{margins}{Numeric vector of length 2 containing the margins for column and row names, respectively. Default is \code{c(5, 4)}.}

  \item{y.lab}{The label of y axis. Default is "Genotypes".}

  \item{key.lab}{The label of color key. Default is "Genotype ranking".}

  \item{resolution}{Valid parameter if \code{file.type = "tiff"}. Define the resolution of the plot. Default is "300".}

}

\examples{
library(WAASB)
model = WAASBYratio(WAASBdata, increment = 5, saveWAASY = 60)

# plot the heatmap considering the ranks in different number of PC axis.
plot.WAASBYratio(model, type = 1)

# plot the heatmap considering the ranks in WAASY/GY ratio and export to a *.tiff file
plot.WAASBYratio(model, type = 2, export = TRUE, file.type = "tiff")

}

\author{Tiago Olivoto \email{tiagoolivoto@gmail.com}}


\seealso{\code{\link{plot.scores}}}
