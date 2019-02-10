\name{autoplot}
\alias{autoplot.WAAS.AMMI}
\alias{autoplot.WAASB}
\title{Set of plots residuals for models of class WAAS.AMMI and WAASB}
\usage{
\method{autoplot}{WAAS.AMMI}(x, conf = 0.95, labels = FALSE, theme = theme_waasb(),
         alpha = 0.2, fill.hist = "gray", col.hist = "black",
         col.point = "black", col.line = "red", col.lab.out = "red",
         size.lab.out = 2.5, size.tex.lab = 10, size.shape = 1.5,
         bins = 30, which = c(1:4), mfrow = c(2, 2), ...)


\method{autoplot}{WAASB}(x, type = "res", conf = 0.95, out = "print",
         labels = FALSE, theme = theme_waasb(),
         alpha = 0.2, fill.hist = "gray", col.hist = "black",
         col.point = "black", col.line = "red", col.lab.out = "red",
         size.lab.out = 2.5, size.tex.lab = 10, size.shape = 1.5,
         bins = 30, which = c(1:4), mfrow = c(2, 2), ...)
}

\description{
Residual plots for a output model of class "WAAS.AMMI" and "WAASB". Six sets of plots are produced: (1) Residuals vs fitted, (2) normal Q-Q plot for the residuals, (3) scale-location plot (standardized residuals vs Fitted Values), (4) standardized residuals vs Factor-levels, (5) Histogram of raw residuals and (6) standardized residuals vs observation order. For a "WAASB" object, normal Q-Q plot for random effects may also be obtained declaring \code{type = "re"}

}

\arguments{

  \item{x}{An object of class \code{WAAS.AMMI} or \code{WAASB}.}

  \item{type}{If \code{type = "re"}, normal Q-Q plots for the random effects are obtained.}

  \item{conf}{Level of confidence interval to use in the Q-Q plot (0.95 by default). }

    \item{out}{How the output is returned. Must be one of the "print" (default) or "return".}

  \item{labels}{Logical argument. If \code{TRUE} labels the points outside confidence interval limits.}

  \item{theme}{The theme to use in the graphics. Default is \code{theme = theme_waasb()}.}


  \item{alpha}{The transparency of confidence band in the Q-Q plot. Must be a number between 0 (opaque) and 1 (full transparency).}

  \item{fill.hist}{The color to fill the histogram. Default is "gray".}

  \item{col.hist}{The color of the border of the the histogram. Default is "black".}

  \item{col.point}{The color of the points in the graphic. Default is "black".}

  \item{col.line}{The color of the lines in the graphic. Default is "red".}

  \item{col.lab.out}{The color of the labels for the "outlying" points.}

  \item{size.tex.lab}{The size of the text in axis text and labels.}

  \item{size.shape}{The size of the shape in the plots.}

  \item{bins}{The number of bins to use in the histogram. Default is 30.}

  \item{which}{Which graphics should be plotted. Default is \code{which = c(1:4)} that means that the first four graphics will be plotted.}

  \item{mfrow}{Allow creating a simple multi-paneled plot. Must be a vector of length 2. The first argument is the number of rows and the second the number of columns. Default is \code{mfrow = c(2, 2)}.}

  \item{...}{Additional parameter for the function}




  }


\examples{
library(METAAB)
model = WAAS.AMMI(WAASBdata, GY, GEN, ENV, REP)
autoplot(model)

model2 = WAASB(WAASBdata, c(GY, HM), GEN, ENV, REP)
lapply(model_oat, function(x){
  autoplot.WAASB(x, which = 2, out = "return")
})

}

\author{Tiago Olivoto \email{tiagoolivoto@gmail.com}}

