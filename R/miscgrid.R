
.polygrid0 <- function(x,y, vpb, gpb, bb)
{
    aa1 <- c(diff(bb[1,]),diff(bb[2,]))
    if (aa1[1]>aa1[2]) {
        denom <- aa1[1]
        numx <- bb[1,1]
        numy <- bb[2,1]-denom/2+aa1[2]/2
    } else {
        denom <- aa1[2]
        numx <- bb[1,1]-denom/2+aa1[1]/2
        numy <- bb[2,1]
    }
    grid.polygon((x-numx)/denom,
                 (y-numy)/denom, vp=vpb, gp=gpb)
}


.polyringhole <- function (Sr, fill, border = NULL,
                          pbg, vpb, bb)
{
    if (!is(Sr, "Polygons"))
        stop("Not an Polygons object")
    pO <- slot(Sr, "plotOrder")
    polys <- slot(Sr, "Polygons")
    for (i in pO) {
        xy <- slot(polys[[i]], "coords")
        if (!slot(polys[[i]], "hole"))
            .polygrid0(xy[,1], xy[,2], vpb, gpb = gpar(col=border, fill=fill), bb)
        else .polygrid0(xy[,1], xy[,2], vpb, gpb=gpar(col=border, fill=pbg), bb)
    }
}


### Cette fonction permet d'utiliser les fonctions de grid
### pour representer les spatialpolygonsdataframe
polygrid <- function (x, vpb, fill="white", border="black", pbg="white", bb=bbox(x))
{
    if (!is(x, "SpatialPolygons"))
        stop("Not a SpatialPolygons object")
    n <- length(slot(x, "polygons"))
    if (length(border) != n)
        border <- rep(border, n, n)
    polys <- slot(x, "polygons")
    pO <- slot(x, "plotOrder")
    if (length(fill) != n)
        fill <- rep(fill, n, n)
    for (j in pO) .polyringhole(polys[[j]], fill=fill[j], border = border[j],
                                pbg = pbg, vpb, bb)
}

### Cette fonction permet d'utiliser les fonctions de grid
### pour representer les images raster pixels ou grid
rastergrid <- function(z, interpolate=FALSE, ...)
{
    if (inherits(z, "SpatialPixelsDataFrame"))
        fullgrid(z) <- TRUE
    x <- as.image.SpatialGridDataFrame(z)$z
    x <- (x-min(na.omit(c(x))))/diff(range(na.omit(c(x))))
    x[is.na(x)] <- 1
    x <- t(x)
    x <- x[nrow(x):1,]
    grid.raster(image=x, interpolate = interpolate,...)
}
