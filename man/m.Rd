\name{m}
\alias{m}
\alias{year}
\alias{month}
\alias{mat2spdf}
\title{
  Various functions used daily.
}
\description{
  \code{m} is just a quick alias for the function \code{as.matrix}.
  \code{year} and \code{month} extract the year and the month
  respectively from a \code{POSIXt} object.  \code{mat2spdf} converts a
  two column matrix to a SpatialPolygons object with one polygon.
  \code{docbind(x)} is a shortcut for \code{do.call(cbind, x)}, and
  \code{dorbind(x)} is a shortcut for \code{do.call(rbind, x)}.
}
\usage{
m(x, ...)
year(x)
month(x)
mat2spdf(x, id="A")
docbind(x)
dorbind(x)
}
\arguments{
  \item{x}{an R object}
  \item{id}{the ID of the polygon}
  \item{\dots}{additionnal arguments to \code{as.matrix}}
}
\value{
  \code{m} returns a matrix. \code{year} and \code{month} return a
  vector. \code{mat2spdf} return a \code{SpatialPolygons} object.
}
\seealso{
  \code{\link{as.matrix}}, \code{\link{POSIXt}}.
}
\examples{
## No examples
}
\keyword{math}
