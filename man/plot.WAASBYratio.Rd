\name{plot.WAASBYratio}
\alias{plot.WAASBYratio}
\title{Plot heat maps with genotype ranking}
\usage{
\method{plot}{WAASBYratio}(x, type, export = FALSE, file.type = "pdf",
     file.name = NULL, width = 6, height = 5, size.lab = 1,
     margins = c(5, 4),  y.lab = NULL, x.lab = NULL,
     key.lab = "Genotype ranking", resolution = 300, ...)
}

\description{
Plot heat maps with genotype ranking in two ways.

}
\arguments{

  \item{x}{The \code{WAASYratio object} }

  \item{type}{\code{1 = Heat map Ranks}: this graphic shows the genotype ranking considering the WAAS estimated with different numbers of Principal Components; \code{2 = Heat map WAASY-GY ratio}: this graphic shows the genotype ranking considering the different combinations in the WAAS/GY ratio. }

  \item{export}{Export (or not) the plot. Default is \code{TRUE}.}

  \item{file.type}{If \code{export = TRUE} define the type of file to be exported. Default is \code{pdf}, Graphic can also be exported in \code{*.tiff} format by declaring  \code{file.type = "tiff"}.}

  \item{file.name}{The name of the file for exportation, default is \code{NULL}, i.e. the files are automatically named.}

  \item{width}{The width "inch" of the plot. Default is \code{8}.}

  \item{height}{The height "inch" of the plot. Default is \code{7}. }

  \item{size.lab}{The label size of the plot. It is suggested attribute 1}

  \item{margins}{Numeric vector of length 2 containing the margins for column and row names, respectively. Default is \code{c(5, 4)}.}

  \item{y.lab}{The label of y axis. Default is "Genotypes".}

  \item{x.lab}{The label of x axis. Default is "Number of axes".}

  \item{key.lab}{The label of color key. Default is "Genotype ranking".}

  \item{resolution}{Valid parameter if \code{file.type = "tiff"}. Define the resolution of the plot. Default is "300".}

  \item{...}{Other arguments of the function}

}

\details{
The first type of heatmap shows the genotype ranking depending on the number of principal component axis used for estimating the WAASB index. An euclidian distance-based dendrogram is used for grouping the genotype ranking for both genotypes and principal component axis. The second type of heatmap shows the genotype ranking depending on the WAASB/GY ratio. The ranks obtained with a ratio of 100/0 considers exclusively the stability for the genotype ranking. On the other hand, a ratio of 0/100 considers exclusively the productivity for the genotype ranking.  Four clusters are estimated (1) unproductive and unstable genotypes; (2) productive, but unstable genotypes; (3) stable, but unproductive genotypes; and (4), productive and stable genotypes.

}

\examples{
library(METAAB)
model = WAASBYratio(WAASBdata, resp = "GY")

# plot the heatmap considering the ranks in different number of PC axis.
plot.WAASBYratio(model, type = 1)

# plot the heatmap considering the ranks in WAASY/GY ratio and export to a *.tiff file
plot.WAASBYratio(model, type = 2, export = TRUE, file.type = "tiff")

}

\author{Tiago Olivoto \email{tiagoolivoto@gmail.com}}


\seealso{\code{\link{plot.scores}}}
