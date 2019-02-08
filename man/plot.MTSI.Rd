\name{plot.MTSI}
\alias{plot.MTSI}
\title{Plot the multitrait stability index}
\usage{
\method{plot}{MTSI}(x, SI = 15, radar = TRUE, ...)
}

\description{
Plot the multitrait stability index using ggplot-generated graphics. .


}
\arguments{

  \item{x}{An object of class \code{MTSI} }

  \item{SI}{An integer [0-100]. The selection intensity in percentage of the total number of genotypes.}

  \item{radar}{Logical argument. If true (default) a radar plot is generated after using \code{coord_polar()}.}


  \item{...}{Other arguments of the function}

}


\author{Tiago Olivoto \email{tiagoolivoto@gmail.com}}

