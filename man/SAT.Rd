\name{SAT}
\alias{SAT}
\alias{invSAT}
\title{
  Solid Angle Transform
}
\description{
  These functions implement the solid angle transform described by
  Friedman and Tukey (1973).  This transform has been recommended to
  allow the optimization in the projection pursuit approach (see
  details).
}
\usage{
SAT(eta)

invSAT(X)
}
\arguments{
  \item{eta}{The point in the unit hypercube in the (N-1)-dimensional
    space.} 
  \item{X}{The point on the surface of the sphere in the N-dimensional
    space.} 
}
\details{
  
  A given direction in the N-dimensional space is defined by a normed
  vector, i.e. a point on a sphere.  To find the direction that
  maximizes a given criterion is a constrained problem (the solution
  should be of length one).  Friedman and Tukey (1973) indicate that the
  SAT approach reversibly maps the surface of a unit sphere in a
  C-dimensional space to a (C-1)-dimensional unit hypercube. Thus this
  transform can be used to transform a constrained optimization problem
  into an unconstrained optimization problem (see \code{?projpurs}).

}
\note{
  The formulas described by Friedman and Tukey (1973) are incorrect for
  the case "odd N".  We have corrected these formulas in this
  implementation. 
}
\references{
  Friedman, J. and Tukey, J. 1973. A projection pursuit for exploratory
  data analysis. \emph{IEEE Transactions on Computers}, c23: 881--890.
}
\author{
  Clement Calenge \email{clement.calenge@oncfs.gouv.fr}
}
\seealso{
  \code{\link{projpurs}} for further details on the projection pursuit
  approach.
}
\examples{

## generates a point in a N-1 dimensional unit hypercube
## with N=9
eta <- runif(8)

## Map it on a sphere in a N-dimensional space
X <- SAT(eta)
sum(X^2) ## the length is 1
length(X) ## 9

## and reversibly
invSAT(X)
eta


}
\keyword{multivariate}
