\name{wildboar}
\alias{wildboar}
\title{
  Draws wild boar at specified coordinates
}
\description{
  This function draws wild boar at specified coordinates on an existing device.
}
\usage{
wildboar(x, y = x, cex = 1, colo = "brown")
}
\arguments{
  \item{x}{a numeric vector giving the x coordinates where the wild boar should be drawn}
  \item{y}{a numeric vector giving the y coordinates where the wild boar should be drawn}
  \item{cex}{numeric value.  The size of the wildboar}
  \item{colo}{character value.  The color of the wildboar}
}
\details{
  The function should be used on an existing graphical device.
}
\value{
  The wildboar is displayed on the graph.
}
\author{
  Clement Calenge \email{clement.calenge@oncfs.gouv.fr}
}
\examples{

plot(1:10, 1:10, ty="n")
wildboar(1:10, cex=3)
}
\keyword{ aplot }
