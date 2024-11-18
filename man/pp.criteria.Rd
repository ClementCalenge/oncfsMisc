\name{pp.criteria}
\alias{pp.criteria}
\alias{pp.holes}
\alias{pp.centralmass}
\alias{pp.FriedmanTukey}
\alias{pp.Friedman1987}
\title{
  Several interesting criteria for one dimensional projection pursuit
}
\description{
  These functions implement several interesting criteria for the one
  dimensional projection pursuit (see Details).
}
\usage{

pp.holes(x, p)

pp.centralmass(x, p)

pp.FriedmanTukey(x, R, prop=0.1, f=expression(R-r))

pp.Friedman1987(x, J)

}
\arguments{
  \item{x}{The vector containing the scores of individuals on a
    given variable.}
  \item{p}{The dimension of the space in which the criterion is to be
    maximized (i.e. the number of columns of the data.frame passed to
    the function \code{projpurs}).}
  \item{R}{maximum distance under which the point-to-point distance is
    calculated (see details)}
  \item{prop}{the proportion of points excluded of the calculation of
    the trimmed standard deviation (see details)}
  \item{f}{R expression describing the monotonic functions of the
    differences between the point to point distance \code{r} and the
    maximum distance \code{R}}
  \item{J}{the degree of the Legendre polynomial expansion used to
    capture the distribution of the data in this criterion.}
}
\details{
  \code{pp.holes} has been used by Cook and Swayne (2007, p. 30) to
  identify the projections with few points in the centre of the plot
  (this criterion is quickly optimized, but will miss a structure if
  groups are not separated at the origin -- this is the case of the
  dataset \code{deug} actually, see examples).\cr
  
  \code{pp.centralmass} has been used by Cook and Swayne (2007, p. 30)
  to identify projections with a lot of point at the origin of the
  plot (this criterion is quickly optimized, and is useful to identify
  outliers).\cr

  \code{pp.FriedmanTukey} is the original criterion of Friedman and
  Tukey (1973). This criterion is made of two components: (i) a global
  variance measured with a trimmed standard deviation (i.e. a standard
  deviation calculated on the data after removal of a given proportion
  of points, controlled by \code{prop}), and (ii) a measure of local
  density of points (the point to point distance \code{r} is calculated,
  considering only the distances lower than \code{R}. Then, the density
  of points is calculated by the sum of monotonic functions of
  \code{R-r}.  By default, the function $f$ is simply the difference
  \code{R-r}, but other functions can be defined).  This criterion is
  generally slow to optimize.\cr
  
  \code{pp.Friedman1987} is the criterion of Friedman (1987), and can be
  used to identify the directions where the deviation from a normal
  distribution occurs at the core of the distribution and not in the
  tails.  Very useful, but: (i) very slow and (ii) affected by
  outliers. This criterion is meaningful only when the argument
  \code{sphere} of the function \code{projpurs} is set to \code{TRUE}.
  This criterion relies on a Legendre polynomials expansion of the
  distribution of \code{x}, and \code{J} controls the maximum degree of
  these polynomials.

}
\value{
  A numeric value containing the value of the criterion.
}
\references{
  Friedman, J. and Tukey, J. 1973. A projection pursuit for exploratory
  data analysis. \emph{IEEE Transactions on Computers}, c23: 881--890.

  Friedman, J. 1987. Exploratory projection pursuit. \emph{Journal of
  the American Statistical Association}, 82: 249--266.  
}
\author{
  Clement Calenge \email{clement.calenge@ofb.gouv.fr}
}
\seealso{
  \code{\link{projpurs}} for additional details on the one-dimensional
  projection pursuit algorithm.
}
\examples{
\dontrun{
## Example with the deug dataset
data(deug)
tab <- deug$tab


## Example with pp.centralmass (p = 9) (par is unused here)
foo <- function(x, par)
     pp.centralmass(x, 9)
set.seed(12) ## to make the calculations reproducible
prp <- projpurs(tab, foo)
plot(prp)


## Example with pp.holes
foo <- function(x, par)
     pp.holes(x, 9)
set.seed(12) ## to make the calculations reproducible
prp <- projpurs(tab, foo)
plot(prp)


## Example with pp.FriedmanTukey
set.seed(12) ## to make the calculations reproducible
foo <- function(x, par)
     pp.FriedmanTukey(x, 0.5)
prp <- projpurs(tab, foo, sphere=TRUE)


## pp.Friedman1987 is demonstrated in the help page of projpurs

}
}
\keyword{multivariate}
