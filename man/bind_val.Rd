\name{bind_val}
\alias{bind_val}
\title{Bind cross-validation objects.}
\usage{
bind_val(...)
}
\description{
Helper function that combines objects of class \code{validation.AMMIF} or \code{validation.blup}. It is useful when looking for a boxplot containing the RMSPD values of both AMMI and BLUP models.
}

\arguments{
  \item{...}{Input objects of class \code{validation.AMMIF} or \code{validation.blup}.}
}

\author{Tiago Olivoto \email{tiagoolivoto@gmail.com}}

\examples{
\dontrun{
library(METAAB)
AMMI = validation.AMMIF(data_ge,
                        resp = GY,
                        gen = GEN,
                        env = ENV,
                        rep = REP,
                        nboot = 100,
                        nrepval = 2)
BLUP = validation.blup(data_ge,
                        resp = GY,
                        gen = GEN,
                        env = ENV,
                        rep = REP,
                        nboot = 100,
                        nrepval = 2)
bind_data = bind_val(AMMI, BLUP)
plot(bind_data)
}
}


