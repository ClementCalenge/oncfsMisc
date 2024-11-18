\name{poissonnessPlot}
\alias{poissonnessPlot}
\title{
  Draw a Poissonness Plot
}
\description{
  This function draws a Poissonness plot as described in Hoaglin and Tukey (2006).
}
\usage{
poissonnessPlot(x, truncate=-1, outer=1.5,...)
}
\arguments{
  \item{x}{a numeric vector for which the Poissonness plot is desired.}
  \item{truncate}{A value for truncation (for example, \code{truncate=0}
    allows to compare the observed distribution with a Poisson
    distribution truncated at 0).}
  \item{outer}{numeric value used to control "outer intervals" on the
    plot (see details). }
  \item{\dots}{additional arguments to be passed to the function \code{plot}.}
}
\details{
  The Poissonness plot can be used to determine whether a variable
  \code{x} is drawn from a Poisson distribution.  Let nk be the number of
  elements of \code{x} taking the value k.  A function phi(nk) is
  represented as a function of k, and is expected to be linear if
  \code{x} is Poisson distributed.  The points (k, phi(nk)) are
  represented as red points on the plot.  To account for the sampling
  variation of nk, Hoaglin and Tukey recommend to calculate a modified
  value \code{nk*}, and to calculate a symetric confidence interval on
  phi(nk*).  The 95\% confidence intervals are represented as thick line
  segments (and the centers nk* are represented as open circle).
  Hoaglin and Tukey recommend to "expand" these confidence intervals to
  account for multiple testing, by multiplying these intervals by 1.5
  (represented by thin segments). We allow for the possibility to change
  this value with the argument  \code{outer}.\cr

  The red line correspond to the linear regression of phi(nk*) as a
  function of k. Then, the function tries to find a line passing through
  all thin segments (blue dashed line). If it does not succeed, no blue
  line is drawn on the plot.

}
\references{
  Hoaglin, D.C and Tukey (2006) Checking the shape of discrete
  distribution. pp. 345--416 in Hoaglin, D.C, Mosteller, F. and Tukey,
  J.W. \emph{Exploring data tables, trends, and shapes}. John Wiley and
  sons. 
}
\author{
  Clement Calenge \email{clement.calenge@ofb.gouv.fr}
}
\examples{
## For reproducibility
set.seed(980)

## Poisson distribution
po <- rpois(500, 2)
poissonnessPlot(po)

## Poisson distribution with parameter lower than 1
poissonnessPlot(rpois(5000, 0.3))

## Poisson distribution with parameter from a log-normal distribution
x <- sapply(1:1000, function(x) rpois(1, exp(1+rnorm(1, 0,0.3))))
poissonnessPlot(x)


}
\keyword{distribution}
\keyword{aplot}
