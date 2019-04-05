\name{plot.MTSI}
\alias{plot.MTSI}
\title{Plot the multi-trait stability index}
\usage{
\method{plot}{MTSI}(x,  SI = 15, radar = TRUE, size.point = 2,
     col.sel = "red", col.nonsel = "black", ...)
}

\description{
Plot the multitrait stability index using ggplot-generated graphics. .


}
\arguments{

  \item{x}{An object of class \code{MTSI} }

  \item{SI}{An integer [0-100]. The selection intensity in percentage of the total number of genotypes.}

  \item{radar}{Logical argument. If true (default) a radar plot is generated after using \code{coord_polar()}.}

  \item{size.point}{The size of the point in graphic.}

  \item{col.sel}{The colour for selected genotypes.}

  \item{col.nonsel}{The colour for nonselected genotypes.}

  \item{...}{Other arguments of the function.}

}


\author{Tiago Olivoto \email{tiagoolivoto@gmail.com}}


\examples{
library(METAAB)
MTSI_MODEL = WAASB(data_ge,
                   env = ENV,
                   gen = GEN,
                   rep = REP,
                   resp = c(GY, HM))

MTSI_index = MTSI(MTSI_MODEL)
plot(MTSI_index)

# Alternatively using the pipe operator \%>\%
library(dplyr)
MTSI_index2 = data_ge \%>\%
              WAASB(ENV, GEN, REP, c(GY, HM)) \%>\%
              MTSI()
plot(MTSI_index2)


}
