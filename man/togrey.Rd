\name{togrey}
\alias{togrey}
\title{
  quickly transforms a vector into a vector of grey levels for colors mapping
}
\description{
  quickly transforms a vector into a vector of grey levels for colors mapping
}
\usage{
togrey(x)
}
\arguments{
  \item{x}{a numeric vector}
}
\value{
  a vector of grey levels
}
\author{
  Clement Calenge \email{clement.calenge@oncfs.gouv.fr}
}
\examples{
plot(1:10, pch=16, col=togrey(1:10))
}
\keyword{ aplot }
