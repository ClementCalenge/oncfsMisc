\name{miscgrid}
\alias{miscgrid}
\alias{s.arrowg}
\alias{s.corcircleg}
\alias{s.labelg}
\alias{s.classg}
\alias{s.valueg}
\alias{rastergrid}
\alias{polygrid}
\title{
  Grid based plotting functions (including spatial object mapping)
}
\description{
  \code{rastergrid} can be used to plot an object of class
  \code{SpatialPixelsDataFrame}.  \code{polygrid} can be used to plot an
  object of class \code{SpatialPolygonsDataFrame}.  The other functions
  are grid-based functions corresponding to plotting functions of the
  package ade4 (i.e. \code{s.arrowg} is the grid-based version of
  \code{s.arrow}, etc.).
}
\usage{
s.arrowg(dfxy, xax = 1, yax = 2, new = TRUE, summarize = 3, axes = TRUE,
         grid = TRUE, xlim = NULL, ylim = NULL, name = "vpc1")

s.corcircleg(dfxy, xax = 1, yax = 2, new = TRUE, summarize = 3, name = "vprt")

s.classg(dfxy, fac, xax = 1, yax = 2, new = TRUE, axes = TRUE, grid = TRUE,
         xlim = NULL, ylim = NULL,
         name = "glou", biplot = NULL, summarize = 3)

s.labelg(dfxy, xax = 1, yax = 2, new = TRUE, biplot = NULL,
         summarize = 3, label = TRUE,
         axes = TRUE, grid = TRUE, xlim = NULL, ylim = NULL,
         name = "datavp")

s.valueg(dfxy, value, xax = 1, yax = 2, new = TRUE, cexmax = 2,
         biplot = NULL, summarize = 3, axes = TRUE, grid = TRUE,
         xlim = NULL, ylim = NULL, name = "datavp")

rastergrid(z, interpolate = FALSE, ...)

polygrid(x, vpb, fill = "white", border = "black", pbg = "white", bb = bbox(x))

}
\arguments{
  \item{dfxy}{a data.frame}
  \item{xax,yax}{the column number of \code{dfxy} that will be used as x
    and y axes.}
  \item{new}{logical value: whether the function will draw a new plot or
    will add the result to an already existing plot.}
  \item{summarize}{By default, the variable names are truncated to
    \code{summarize} characters before plotting.  Alternatively, if
    \code{summarize} is negative, the original variable names will be
    plotted.}
  \item{axes}{logical value: whether to add axes to the plot.}
  \item{grid}{logical value: whether to add a grid to the plot.}
  \item{xlim,ylim}{vectors of length two containing the x and y limits
    of the plots.}
  \item{cexmax}{The maximum size for the values to be plotted.}
  \item{name}{the name of the viewport created by the function}
  \item{label}{logical value: whether the labels should be added.}
  \item{fac}{a factor defining the classes.}
  \item{value}{a vector defining the values to be plotted.}
  \item{biplot}{optionally, a data.frame containing the coordinates of
    arrows that should be added to the plot.}
  \item{z}{an object of class SpatialPixelsDataFrame (package sp).}
  \item{interpolate}{logical value. Corresponds to the parameter
    \code{interpolate} of the function \code{grid.raster} of the package
    grid.}
  \item{x}{an object of class SpatialPolygons (package sp).}
  \item{vpb}{the viewport on which the SpatialPolygons should be
    plotted.}
  \item{fill}{The colors of the polygon.}
  \item{border}{the colors of the border.}
  \item{pbg}{The colour of the background.}
  \item{bb}{The bounding box of the plot.}
  \item{\dots}{additionnal arguments to be passed
    to the function \code{raster.grid} used internally by
    \code{rastergrid}.} 
}
\value{
  These functions return, invisibly, the viewport containing the plot
}
\author{
  Clement Calenge \email{clement.calenge@oncfs.gouv.fr}
}
\seealso{
  \code{\link[ade4]{s.value}}, \code{\link[ade4]{s.label}},
  \code{\link[ade4]{s.class}}, \code{\link[ade4]{s.corcircle}},
  \code{\link[ade4]{s.arrow}}
}
\examples{
library(ade4)
data(deug)
pc <- dudi.pca(deug$tab, nf=2, scannf = FALSE)
s.corcircleg(pc$co)
s.labelg(pc$li)

}
\keyword{spatial}
\keyword{aplot}
\keyword{hplot}
