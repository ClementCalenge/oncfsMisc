\name{mvrspherunif}
\alias{mvrspherunif}
\title{
  Random sampling of directions in a N-dimensional space
}
\description{
  This function samples random directions in a N-dimensional space.
}
\usage{
mvrspherunif(n=1, dim)
}
\arguments{
  \item{n}{The number of vectors that should be simulated.}
  \item{dim}{The dimension of the space where the directions should be sampled.} 
}
\details{
  A direction in a N-dimensional space is defined by a unit vector,
  i.e. by the coordinates of a point on the surface of a unit sphere.
  The function \code{mvspherunif} samples \code{n}  points uniformly on
  this surface.  This is done by sampling point in a N-dimensional
  hypercube, and by removing the vectors with a length greater than 1. 
}
\value{
  a matrix with \code{n} columns and \code{dim} rows, containing the
  random vectors.
}
\author{
  Clement Calenge \email{clement.calenge@oncfs.gouv.fr}
}
\examples{

## Sample 5 random vectors in a 10-dimensional space
a <- mvrspherunif(5, 10)

## the vectors are 10-dimensional, and of length 1
crossprod(a)

}
\keyword{multivariate}
