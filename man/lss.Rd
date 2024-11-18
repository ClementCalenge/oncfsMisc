\name{lss}
\alias{lss}
\title{
  List objects belonging to a given class
}
\description{
  This function lists the objects belonging to a given class present in
  the environment (as the function \code{ls}, but restricted to a given
  class). 
}
\usage{
lss(cls = "data.frame", ...)
}
\arguments{
  \item{cls}{ The class required}
  \item{\dots}{ additionnal elements passed to the function \code{ls}}
}
\author{
  Clement Calenge \email{clement.calenge@ofb.gouv.fr}
}
\seealso{
  \code{\link{ls}}
}
\examples{
foo <- function(x) 2*x

lss("function")

}
\keyword{ data }
