\name{rbind_fill}
\alias{rbind_fill}
\title{Helper function for combining data.frames by row, filling missing values}
\usage{
rbind_fill(..., fill = ".")
}
\description{
Helper function similar to \code{plyr::rbind.fill} that combines data.frames by row and fills with \code{.} missing values
}

\arguments{

\item{...}{Input dataframes.}

\item{fill}{What use to fill? Default is \code{"."}}

}

\author{Tiago Olivoto \email{tiagoolivoto@gmail.com}}

\examples{
df1 = data.frame(v1 = c(1, 2), v2 = c(2, 3))
df2 = data.frame(v3 = c(4, 5))
rbind_fill(df1, df2)
rbind_fill(df1, df2, fill = "NA")
}


