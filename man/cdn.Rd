\name{cdn}
\alias{cdn}
\title{
  Spatial smoothing with a loess regression.
}
\description{
  This function smooths the value of a variable measured at points with
  a loess regression over space.
}
\usage{
cdn(xy, z, limits = NULL, span = 0.25, le = 100, degree = 2, \dots)
}
\arguments{
  \item{xy}{An object inheriting the class \code{SpatialPoints} from the package sp.}
  \item{z}{A numerical vector containing the values measured at points
    stored in \code{xy}}
  \item{limits}{Optionally, an object inheriting the class \code{SpatialPolygons}
    from the package sp, definining the limits for the interpolation.}
  \item{span}{the span parameter of the loess regression (see \code{?loess})}
  \item{le}{the resolution of the smoothed map (a (le x le) grid is
    built and the variable is predicted for each pixel).}
  \item{degree}{the degree of the loess regression (see \code{?loess})}
  \item{\dots}{additional arguments passed to the function \code{loess}}
}
\details{
  The loess prediction will be carried out only within the range of
  the x and y coordinates of the points in \code{xy} (and within a
  shorter range if \code{limits} is fixed).  In other words, no
  extrapolation will be carried out.
}
\value{
  An object of class \code{SpatialPixelsDataFrame}
}
\author{
  Clement Calenge \email{clement.calenge@ofb.gouv.fr}
}
\seealso{
  \code{\link{loess}}
}
\examples{
## Unsound example just to demonstrate how to use the function
## in practice
## (unsound because the variable has already measured along the pixels
## of a grid, so that loess is not really required here).

data(meuse)
coordinates(meuse) <- ~x+y

r1 = cbind(c(180114, 180553, 181127, 181477, 181294, 181007, 180409, 
180162, 180114), c(332349, 332057, 332342, 333250, 333558, 333676, 
332618, 332413, 332349))
sr1=Polygons(list(Polygon(r1)),"r1")
sr=SpatialPolygons(list(sr1))

aa <- cdn(meuse, meuse[[1]], sr, le=500)
image(aa)

}
\keyword{ aplot }
\keyword{ spatial }
