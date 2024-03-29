% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/corr_coef.R
\name{network_plot}
\alias{network_plot}
\title{Network plot of a correlation matrix}
\usage{
network_plot(
  model,
  min_cor = NULL,
  show = c("signif", "all"),
  p_val = 0.05,
  legend = c("full", "range"),
  colours = c("red", "white", "blue"),
  legend_width = 1,
  legend_height = 15,
  legend_position = c("right", "left", "top", "bottom"),
  curved = TRUE,
  angle = 90,
  curvature = 0.5,
  expand_x = 0.25,
  expand_y = 0.25
)
}
\arguments{
\item{model}{A model computed with \code{\link[=corr_coef]{corr_coef()}} or a symmetric matrix, often
produced with \code{\link[stats:cor]{stats::cor()}}.}

\item{min_cor}{Number to indicate the minimum value of correlations to plot
(0-1 in absolute terms). By default, all the correlations are plotted when
\code{model} is a matrix, and significant correlations (p-value < 0.05) when
\code{model} is an object computed with \code{\link[=corr_coef]{corr_coef()}}.}

\item{show}{The correlations to be shown when \code{model} is an object computed
with \code{\link[=corr_coef]{corr_coef()}}. Either \code{"signif"} (default) to show only significant
correlations or \code{"all"} to show all the correlations.}

\item{p_val}{The p-value to indicate significant correlations. Defaults to
\code{0.05}.}

\item{legend}{The type of legend. Either \code{"full"} (ranges from -1 to +1) or
\code{"range"} (ranges according to the data range). Defaults to \code{"full"}.}

\item{colours}{A vector of colors to use for n-color gradient.}

\item{legend_width}{The width of the legend (considering \code{position = "right"})}

\item{legend_height}{The height of the legend (considering \code{position = "right"})}

\item{legend_position}{The legend position. Defaults to \code{"right"}.}

\item{curved}{Shows curved paths? Defaults to \code{TRUE}.}

\item{angle}{A numeric value between 0 and 180, giving an amount to skew the
control points of the curve. Values less than 90 skew the curve towards the
start point and values greater than 90 skew the curve towards the end
point.}

\item{curvature}{A numeric value giving the amount of curvature. Negative
values produce left-hand curves, positive values produce right-hand curves,
and zero produces a straight line.}

\item{expand_x, expand_y}{Vector of multiplicative range expansion factors. If
length 1, both the lower and upper limits of the scale are expanded
outwards by mult. If length 2, the lower limit is expanded by \code{mult[1]} and
the upper limit by \code{mult[2]}.}
}
\value{
A \code{ggplot} object
}
\description{
Produces a network plot of a correlation matrix or an object computed with
\code{\link[=corr_coef]{corr_coef()}}. Variables that are more highly correlated appear closer
together and are joined by stronger (more opaque) and wider paths.  The proximity of the
points is determined using multidimensional clustering, also known as
principal coordinates analysis (Gower, 1966). The color of the paths also
indicates the sign of the correlation (blue for positive and red for
negative).
}
\examples{
cor <- corr_coef(iris)
network_plot(cor)
network_plot(cor,
             show = "all",
             curved = FALSE,
             legend_position = "bottom",
             legend = "range")

}
\references{
Gower, J.C. 1966. Some Distance Properties of Latent Root and Vector Methods
Used in Multivariate Analysis. Biometrika 53(3/4): 325–338.
\doi{10.2307/2333639}
}
