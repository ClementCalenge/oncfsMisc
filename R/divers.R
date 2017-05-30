year <- function(x) as.POSIXlt(x)$year+1900

month <- function(x) as.POSIXlt(x)$mon+1

mat2spdf <- function(x, id="A")
{
    SpatialPolygons(list(Polygons(list(Polygon(x)), ID=id)))
}
