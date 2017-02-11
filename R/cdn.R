## fonction reste a tester
cdn <- function(xy,z, limits=NULL, span=0.25, le=100, degree=2, ...){

    if (!inherits(xy, "SpatialPoints"))
        stop("xy should inherit the class SpatialPoints")
    if (!is.null(limits)) {
        if (!inherits(limits, "SpatialPolygons"))
            stop("limits should inherits the class SpatialPolygon")
        limits <- as(limits, "SpatialPolygons")
        if (length(limits)>1)
            stop("when specified, limits should contain only one polygon")
    }

    ## preparation des donnees
    xy<-as.data.frame(coordinates(xy))
    names(xy)<-c("x","y")

    if (length(z)!=nrow(xy))
        stop("z should be of the same length as xy")
    w = cbind.data.frame(xy,z)

    ## On effectue la regression loess
    lo=loess(z~x+y,data=w,span=span, degree=degree, ...)

    ## On construit la grille de pixels sur laquelle on va
    ## predire les resultats
    rgx<-range(xy[,1])
    rgy<-range(xy[,2])
    ma<-max(c(rgx[2]-rgx[1], rgy[2]-rgy[1]))
    xg = seq(rgx[1],rgx[1]+ma,le=le)
    yg = seq(rgy[1],rgy[1]+ma,le=le)
    gr=expand.grid(xg, yg)
    names(gr)=names(xy)

    ## On effectue la prediction:
    mod = predict(lo,newdata=gr)

    ## On transforme le resultat en spatialPixelsDataFrame
    df <- data.frame(prediction=as.vector(mod))
    coordinates(df) <- cbind(gr[,1],gr[,2])
    gridded(df) <- TRUE

    ## On ne conserve que les zones a l'interieur de la zone d'etude
    ## (definie par poly4)
    if (!is.null(limits)) {
        ov <- over(df, limits)
        df <- df[!is.na(ov),]
    }

    ## le resultat
    return(df)
}
