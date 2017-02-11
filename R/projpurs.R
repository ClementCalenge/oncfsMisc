projpurs <- function (df, foo, par = NULL, centsca = TRUE, sphere = FALSE,
                      nrep = 100, maxit = 10000,
                      violent = FALSE, triesIfViolent=10)
{
    ## Check that df is a data.frame
    if (!inherits(df, "data.frame"))
        stop("df should be a data.frame")

    ## preliminary PCA
    pc <- dudi.pca(df, scannf=FALSE, nf=ncol(df))

    ## ta2b stores the centred and scaled table
    ta2b <- pc$tab

    ## The number of axes that can be found using this approach is
    ## the dimension of the space - 1 (the last one is determined by the
    ## previous ones)
    nf <- sum(pc$eig>0.0000001)-1

    ## If the data should be sphered
    if (sphere) {
        ta <- pc$l1
    } else {
        if (centsca) {
            ta <- pc$tab
        } else {
            ta <- df
        }
    }

    ## Checks that the function has only two parameters named x and par
    art <- names(as.list(args(foo)))
    if ((!all(art[1:2] == c("x", "par"))) | (length(art) != 3))
        stop("The function foo should have only two arguments named x and par")

    ## This is not used actually (the underlying c function will sample a random
    ## direction as starting point for the algorithm -- if init=0, it uses the best
    ## PCA axis, which may actually pose problem since in this case, the axis will contain
    ## zero coordinates, which are not correctly managed by the sat transform (try
    ## invSAT(c(1,0,0,0)))
    init <- 1

    ## Assigns the function to the environment that will be used for the evaluation of the
    ## function
    e1 <- new.env()
    assign("foo", foo, envir = e1)


    if (violent) {
        fin <- .Call("ProjectionPursuitViolent", ta, quote(foo(x, par)),
                     par, e1, nrep, maxit,
                     nf, triesIfViolent, PACKAGE="oncfsMisc")
    } else {
        fin <- .Call("ProjectionPursuit", ta, quote(foo(x, par)),
                     par, e1, nrep, maxit,
                     nf, init, PACKAGE="oncfsMisc")
    }

    ## Format the output
    fin <- lapply(fin, function(x) {
        re <- as.data.frame(x[[2]][[2]][[2]])
        names(re) <- letters[1:ncol(re)]
        rt = list(sim=x[[3]], obs=x[[4]])
        return(list(axe=x[[1]],
                    valpro=x[[2]][[2]][[1]],
                    vecpro=re,
                    rt=rt))
    })

    ## Reconstructs original axes
    c1 <- list()
    mat <- diag(rep(1,ncol(fin[[1]]$vecpro)))
    for (i in 1:nf) {
        U <- mat%*%as.matrix(fin[[i]]$vecpro)
        a <- fin[[i]]$axe
        v <- U%*%a
        c1[[i]] <- v
        aat <- a%*%t(a)
        mat <- U%*%(diag(rep(1, ncol(U))) - aat)
    }
    c1 <- as.data.frame(c1)
    ## if sphere, back to the original data
    if (sphere) {
        lambmud <- diag(1/sqrt(pc$eig))
        c1 <- as.matrix(pc$c1)%*%lambmud%*%as.matrix(c1)
        c1 <- as.data.frame(c1)
    }
    names(c1) <- paste("Axis", 1:ncol(c1), sep="")
    row.names(c1) <- names(pc$tab)

    ## Results
    res <- list()
    res$tab <- pc$tab
    res$c1 <- c1
    res$li <- as.data.frame(as.matrix(res$tab) %*% as.matrix(res$c1))
    res$cor <- as.data.frame(cor(ta2b, res$li))
    res$criterion <- foo
    res$par <- par
    res$call <- match.call()
    res$sphere <- sphere
    res$maxit <- maxit
    res$randtest <- lapply(fin, function(x) x$rt)
    res$val.crit <- sapply(res$randtest, function(x) x$obs)
    class(res) <- c("oprpu","ppdim1")

    return(res)
}



.showtext <- function(x,y,lab, colo="black") {
    re <- grid.rect(unit(x, "native"),
                    unit(y, "native"),
                    width=stringWidth(lab),
                    height=stringHeight(lab),
                    gp=gpar(fill="white", col=colo),
                    name=paste("rectlab.", lab,sep=""))
    te <- grid.text(lab, x=unit(x,"native"),
                    y=unit(y,"native"),
                    name=paste("lablab.", lab,sep=""),gp=gpar(fontsize=8, col=colo))
}



plot.oprpu <- function(x, xax = 1, yax = 2, nam="sprpu", summarize = 3,...)
{
    if (!inherits(x, "oprpu"))
        stop("x should inherit the class oprpu")


    grid.newpage()
    vp <- viewport(name="gen1")
    pushViewport(vp)
    princ <- viewport(x=0.34, y=0.01, width=0.65, height=0.65,
                      just=c("left","bottom"), name="princ")
    basgauche <- viewport(x=0.01, y=0.01, width=0.32, height=0.32,
                          just = c("left","bottom"), name="basgauche")
    milieugauche <- viewport(x=0.01, y=0.34, width=0.32, height=0.32,
                             just = c("left","bottom"), name="milieugauche")
    hautgauche <- viewport(x=0.01, y=0.67, width=0.32, height=0.32,
                             just = c("left","bottom"), name="hautgauche")
    hautmilieu <- viewport(x=0.34, y=0.67, width=0.32, height=0.32,
                           just = c("left","bottom"), name="hautmilieu")
    hautdroite <- viewport(x=0.67, y=0.67, width=0.32, height=0.32,
                           just = c("left","bottom"), name="hautdroite")

    pushViewport(princ)
    s.labelg(x$li, xax, yax, ...)

    ## quantile plot
    popViewport(2)
    pushViewport(basgauche)
    y <- sort(x$li[,yax])
    qu <- c(1:length(y))/length(y)
    quplot2 <- dataViewport(qu, y, extension=0.05,
                            name="quplot2")
    pushViewport(quplot2)


    ## ajout de la grille
    xgri <- c(seq(0, 1, by=0.2))
    ygri <- seq(min(y), max(y), length=6)
    tmp <- lapply(1:length(xgri), function(i) {
        grid.lines(unit(c(xgri[i],xgri[i]),"native"),
                   unit(c(min(y), max(y)), "native"),
                   gp=gpar(col="grey"),
                   name=paste("qp2grid.x.",i,sep=""))
    })
    tmp <- lapply(1:length(ygri), function(i) {
        grid.lines(unit(c(0, 1),"native"),
                   unit(c(ygri[i],ygri[i]), "native"),
                   gp=gpar(col="grey"),
                   name=paste("qp2grid.y.",i,sep=""))
    })
    grid.rect()

    ## Et des points:
    grid.points(qu, y, size=unit(0.5, "char"),
                name="qp2.points")
    grid.text(paste("Axis", yax), unit(min(qu), "native"), unit(max(y),"native"),
              gp=gpar(fontsize=15, col="red"), just = c("left","top"))

    ## Quantile plot du premier axe
    popViewport(2)
    pushViewport(milieugauche)
    y <- sort(x$li[,xax])
    qu <- c(1:length(y))/length(y)
    quplot1 <- dataViewport(qu, y, extension=0.05,
                            name="quplot1")
    pushViewport(quplot1)


    ## ajout de la grille
    xgri <- c(seq(0, 1, by=0.2))
    ygri <- seq(min(y), max(y), length=6)
    tmp <- lapply(1:length(xgri), function(i) {
        grid.lines(unit(c(xgri[i],xgri[i]),"native"),
                   unit(c(min(y), max(y)), "native"),
                   gp=gpar(col="grey"),
                   name=paste("qp1grid.x.",i,sep=""))
    })
    tmp <- lapply(1:length(ygri), function(i) {
        grid.lines(unit(c(0, 1),"native"),
                   unit(c(ygri[i],ygri[i]), "native"),
                   gp=gpar(col="grey"),
                   name=paste("qp1grid.y.",i,sep=""))
    })
    grid.rect()

    ## Et des points:
    grid.points(qu, y, size=unit(0.5, "char"),
                name="qp1.points")
    grid.text(paste("Axis", xax), unit(min(qu), "native"), unit(max(y),"native"),
              gp=gpar(fontsize=15, col="red"), just = c("left","top"))

    ## Resultat des simulations
    popViewport(2)
    pushViewport(hautgauche)
    if (inherits(x, "ppdim1")) {
        tmp <- lapply(1:length(x$randtest), function(i) {
            den <- density(x$randtest[[i]]$sim)
            xp <- c(den$x, rev(den$x))
            yp <- i+c(den$y/(3*max(den$y)), -rev(den$y/(3*max(den$y))))
            return(list(x=xp, y=yp))
        })
        yli <- (range(c(unlist(lapply(tmp, function(y) y$x)),
                        unlist(lapply(x$randtest, function(y) y$obs)))))
                                        #    yli[1] <- yli[1]-0.01*diff(range(yli))
                                        #    yli[2] <- yli[2]+0.01*diff(range(yli))
        vprt <- dataViewport(c(0,length(x$randtest)+1), yli, name="vprt")
        pushViewport(vprt)
        tmp2 <- lapply(tmp, function(y) {
            grid.polygon(unit(y$y, "native"), unit(y$x, "native"), gp=gpar(fill="grey"))
        })
        grid.points(unit(1:length(x$randtest), "native"),
                    unit(sapply(x$randtest, function(z) z$obs), "native"),
                    pch=21, gp=gpar(fill="red"))
        grid.rect()
        grid.text("Crit.", unit(length(x$randtest)+1, "native"), unit(yli[2],"native"),
                  gp=gpar(fontsize=15, col="red"), just = c("right","top"))
    } else {
        den <- density(c(x$randtest$sim, x$randtest$obs))
        xp <- c(den$x, den$x[1])
        yp <- c(den$y, den$y[1])
        vprt <- dataViewport(xp, yp, name="vprt")
        pushViewport(vprt)
        grid.polygon(unit(xp, "native"), unit(yp, "native"), gp=gpar(fill="grey"))
        grid.lines(unit(c(x$randtest$obs,x$randtest$obs), "native"),
                   unit(c(max(yp)/4,0), "native"),
                   arrow = arrow())
        grid.rect()
    }


    ## Cercle des correlations
    popViewport(2)
    pushViewport(hautmilieu)
    s.arrowg(x$cor, xax, yax, summarize)
    grid.text("$cor", 0, 1,
              gp=gpar(fontsize=15, col="red"), just = c("left","top"))


    ## Et les scores des colonnes
    popViewport(2)
    pushViewport(hautdroite)
    s.arrowg(x$c1, xax, yax, summarize)
    grid.text("$c1", 0, 1,
              gp=gpar(fontsize=15, col="red"), just = c("left","top"))

    ## Et on revient au graphe principal
    popViewport(2)
    pushViewport(princ)
    return(invisible(princ))
}


print.oprpu <- function (x, ...)
{
    cat("Object of class prpu\n")
    cat("Contains:\n")
    print(names(x))
}







s.arrowg <- function(dfxy, xax=1, yax=2, new=TRUE, summarize=3, axes=TRUE, grid=TRUE,
                     xlim = NULL, ylim = NULL, name="vpc1")
{
    if (!inherits(dfxy, "data.frame"))
        stop("dfxy should be a data.frame")

    dfxy <- dfxy[,c(xax,yax)]
    namm <- row.names(dfxy)
    if (summarize>0)
        namm <- substr(namm,1, summarize)

    ra <- apply(dfxy,2,range)
    ra[1,1] <- min(c(ra[1,1], -0.1))
    ra[1,2] <- min(c(ra[1,2], -0.1))
    di2 <- apply(ra,2,diff)
    di <- max(di2)
    whi <- which.min(di2)
    aju <- (di-di2[whi])/2
    ra[1,whi] <- ra[1,whi]-aju
    ra[2,whi] <- ra[2,whi]+aju
    ra[1,] <- ra[1,]-di*0.1
    ra[2,] <- ra[2,]+di*0.1


    if (!is.null(xlim)) {
        ra[,1] <- xlim
    }
    if (!is.null(ylim)) {
        ra[,2] <- ylim
    }

    vpc1 <- dataViewport(ra[,1], ra[,2], name=name, extension = 0)
    pushViewport(vpc1)
    if (new)
        grid.rect(gp=gpar(fill="white"))

    if (grid) {
        inte <- 0.25
        xgri <- c(rev(seq(0, ra[1,1], by=-inte)),
                  seq(0, ra[2,1], by=inte))
        ygri <- c(rev(seq(0, ra[1,2], by=-inte)),
                  seq(0, ra[2,2], by=inte))
        tmp <- lapply(1:length(xgri), function(i) {
            grid.lines(unit(c(xgri[i],xgri[i]),"native"),
                       unit(c(ra[1,2], ra[2,2]), "native"),
                       gp=gpar(col="grey"),
                       name=paste("c1grid.x.",i,sep=""))
        })
        tmp <- lapply(1:length(ygri), function(i) {
            grid.lines(unit(c(ra[1,1], ra[2,1]),"native"),
                       unit(c(ygri[i],ygri[i]), "native"),
                       gp=gpar(col="grey"),
                       name=paste("c1grid.y.",i,sep=""))
        })
        grid.rect()

    }

    if (axes) {
        if (ra[1,2]<0&ra[2,2]>0) {
            grid.lines(unit(c(ra[1,1], ra[2,1]),"native"),
                       unit(c(0,0), "native"), name="xaxiscco")
        }
        if (ra[1,1]<0&ra[2,1]>0) {
            grid.lines(unit(c(0,0), "native"),
                       unit(c(ra[1,2], ra[2,2]),"native"),
                       name="yaxiscco")
        }
    }


    tmp <- sapply(1:nrow(dfxy), function(i) grid.lines(unit(c(0,dfxy[i,1]),"native"),
                                                       unit(c(0,dfxy[i,2]),"native"),
                                                       arrow=arrow(length=unit(0.1,"native"))))
    tmp <- sapply(1:nrow(dfxy), function(i) {
        .showtext(dfxy[i,1], dfxy[i,2], namm[i])
    })

    return(invisible(vpc1))
}


s.corcircleg <- function(dfxy, xax=1, yax=2, new = TRUE, summarize=3, name="vprt")
{
    if (!inherits(dfxy, "data.frame"))
        stop("dfxy should be a data.frame")
    dfxy <- dfxy[,c(xax,yax)]
    namm <- row.names(dfxy)
    if (summarize>0)
        namm <- substr(namm,1, summarize)
    vpcco <- dataViewport(c(-1,1), c(-1,1), name=name)
    pushViewport(vpcco)
    if (new)
        grid.rect(gp=gpar(fill="white"))
    xgri <- c(seq(-1, 1, by=0.2))
    ygri <- seq(-1, 1, by=0.2)
    tmp <- lapply(1:length(xgri), function(i) {
        grid.lines(unit(c(xgri[i],xgri[i]),"native"),
                   unit(c(sin(-acos(xgri[i])), sin(acos(xgri[i]))), "native"),
                   gp=gpar(col="grey"),
                   name=paste("ccogrid.x.",i,sep=""))
    })
    tmp <- lapply(1:length(ygri), function(i) {
        grid.lines(unit(c(-cos(asin(ygri[i])), cos(asin(ygri[i]))),"native"),
                   unit(c(ygri[i],ygri[i]), "native"),
                   gp=gpar(col="grey"),
                   name=paste("ccogrid.y.",i,sep=""))
    })
    grid.rect()


    grid.lines(unit(c(-1, 1),"native"),
               unit(c(0,0), "native"), name="xaxiscco")
    grid.lines(unit(c(0,0), "native"),
               unit(c(-1, 1),"native"),
               name="yaxiscco")
    uuy <- sin(c(seq(pi, -pi, length=100)))
    uux <- cos(c(seq(pi, -pi, length=100)))
    grid.polygon(x=unit(uux,"native"),y=unit(uuy,"native"))

    tmp <- sapply(1:nrow(dfxy), function(i) grid.lines(unit(c(0,dfxy[i,1]),"native"),
                                                       unit(c(0,dfxy[i,2]),"native"),
                                                       arrow=arrow(length=unit(0.1,"native"))))

    tmp <- sapply(1:nrow(dfxy), function(i) {
        .showtext(dfxy[i,1], dfxy[i,2], namm[i])
    })

    return(invisible(vpcco))
}




# s.arrowg(biplot, xax, yax, new=FALSE, axes=FALSE, grid=FALSE)

s.classg <- function(dfxy, fac, xax=1, yax=2, new=TRUE, axes=TRUE, grid=TRUE,
                     xlim=NULL, ylim=NULL,
                     name="glou", biplot=NULL, summarize=3)
{
    if (!inherits(dfxy, "data.frame"))
        stop("dfxy should be a data.frame")
    dfxy <- dfxy[,c(xax,yax)]
    ra <- apply(dfxy,2,range)
    di2 <- apply(ra,2,diff)
    di <- max(di2)
    whi <- which.min(di2)
    aju <- (di-di2[whi])/2
    ra[1,whi] <- ra[1,whi]-aju
    ra[2,whi] <- ra[2,whi]+aju
    ra[1,] <- ra[1,]-di*0.1
    ra[2,] <- ra[2,]+di*0.1
    if (!is.null(xlim)) {
        ra[,1] <- xlim
    }
    if (!is.null(ylim)) {
        ra[,2] <- ylim
    }

    ## Le viewport
    vp2 <- dataViewport(xscale=ra[,1], yscale=ra[,2], extension=0,
                       name=name)
    pushViewport(vp2)
    if (new)
        grid.rect(gp=gpar(fill="white"))

    ## la grille
    if (grid) {
        inte <- round(di/5)
        xgri <- c(rev(seq(0, ra[1,1], by=-inte)),
                  seq(0, ra[2,1], by=inte))
        ygri <- c(rev(seq(0, ra[1,2], by=-inte)),
                  seq(0, ra[2,2], by=inte))
        tmp <- lapply(1:length(xgri), function(i) {
            grid.lines(unit(c(xgri[i],xgri[i]),"native"),
                       unit(c(ra[1,2], ra[2,2]), "native"),
                       gp=gpar(col="grey"),
                       name=paste("grid.x.",i,sep=""))
        })
        tmp <- lapply(1:length(ygri), function(i) {
            grid.lines(unit(c(ra[1,1], ra[2,1]),"native"),
                       unit(c(ygri[i],ygri[i]), "native"),
                       gp=gpar(col="grey"),
                       name=paste("grid.y.",i,sep=""))
        })
        grid.rect()
    }

    ## les axes:
    if (axes) {
        if (ra[1,2]<0&ra[2,2]>0) {
            grid.lines(unit(c(ra[1,1], ra[2,1]),"native"),
                       unit(c(0,0), "native"), name="xaxis")
        }
        if (ra[1,1]<0&ra[2,1]>0) {
            grid.lines(unit(c(0,0), "native"),
                       unit(c(ra[1,2], ra[2,2]),"native"),
                       name="yaxis")
        }
    }

    ## Les points sur le graphique
    mox <- sapply(1:nrow(dfxy), function(i) { mean(dfxy[fac==fac[i],1])})
    moy <- sapply(1:nrow(dfxy), function(i) { mean(dfxy[fac==fac[i],2])})
    grid.segments(unit(mox,"native"), unit(moy,"native"), unit(dfxy[,1],"native"), unit(dfxy[,2],"native"))
    co <- apply(dfxy,2,function(x)tapply(x,fac,mean))
    tmp <- sapply(1:nrow(co), function(i) {
        .showtext(co[i,1],co[i,2], row.names(co)[i])
    })

    ## Un eventuel biplot:
    if (!is.null(biplot)) {
        s.arrowg(biplot, xax, yax, new=FALSE, summarize, axes=FALSE, grid=FALSE, xlim=ra[,1],
                 ylim=ra[,2])
    }

    return(invisible(vp2))
}


s.labelg <- function(dfxy, xax=1, yax=2, new=TRUE, biplot=NULL, summarize = 3, label=TRUE,
                     axes = TRUE, grid = TRUE, xlim = NULL, ylim=NULL,
                     name="datavp")
{
    if (!inherits(dfxy, "data.frame"))
        stop("dfxy should be a data.frame")
    dfxy <- dfxy[,c(xax,yax)]
    ra <- apply(dfxy,2,range)
    di2 <- apply(ra,2,diff)
    di <- max(di2)
    whi <- which.min(di2)
    aju <- (di-di2[whi])/2
    ra[1,whi] <- ra[1,whi]-aju
    ra[2,whi] <- ra[2,whi]+aju
    ra[1,] <- ra[1,]-di*0.1
    ra[2,] <- ra[2,]+di*0.1
    if (!is.null(xlim)) {
        ra[,1] <- xlim
    }
    if (!is.null(ylim)) {
        ra[,2] <- ylim
    }

    ## Le viewport
    vp2 <- dataViewport(xscale=ra[,1], yscale=ra[,2], extension=0,
                       name=name)
    pushViewport(vp2)
    if (new)
        grid.rect(gp=gpar(fill="white"))

    ## la grille
    if (grid) {
        inte <- round(di/5)
        xgri <- c(rev(seq(0, ra[1,1], by=-inte)),
                  seq(0, ra[2,1], by=inte))
        ygri <- c(rev(seq(0, ra[1,2], by=-inte)),
                  seq(0, ra[2,2], by=inte))
        tmp <- lapply(1:length(xgri), function(i) {
            grid.lines(unit(c(xgri[i],xgri[i]),"native"),
                       unit(c(ra[1,2], ra[2,2]), "native"),
                       gp=gpar(col="grey"),
                       name=paste("grid.x.",i,sep=""))
        })
        tmp <- lapply(1:length(ygri), function(i) {
            grid.lines(unit(c(ra[1,1], ra[2,1]),"native"),
                       unit(c(ygri[i],ygri[i]), "native"),
                       gp=gpar(col="grey"),
                       name=paste("grid.y.",i,sep=""))
        })
        grid.rect()
    }

    ## les axes:
    if (axes) {
        if (ra[1,2]<0&ra[2,2]>0) {
            grid.lines(unit(c(ra[1,1], ra[2,1]),"native"),
                       unit(c(0,0), "native"), name="xaxis")
        }
        if (ra[1,1]<0&ra[2,1]>0) {
            grid.lines(unit(c(0,0), "native"),
                       unit(c(ra[1,2], ra[2,2]),"native"),
                       name="yaxis")
        }
    }
    ## Les points sur le graphique
    grid.points(dfxy[,1], dfxy[,2], pch=21,
                name="pointslab", size = unit(0.6, "char"), gp=gpar(fill="white"))

    if (!is.null(biplot)) {
        s.arrowg(biplot, xax, yax, new=FALSE, summarize, axes=FALSE, grid=FALSE, xlim=ra[,1],
                 ylim=ra[,2])
    }

    ##
    if (label) {
        tmp <- sapply(1:nrow(dfxy), function(i) {
            .showtext(dfxy[i,1],dfxy[i,2],row.names(dfxy)[i])
        })
    }


    return(invisible(vp2))
}







s.valueg <- function(dfxy, value, xax=1, yax=2, new=TRUE, cexmax=2, biplot=NULL,
                     summarize = 3, axes=TRUE, grid=TRUE, xlim = NULL, ylim = NULL,
                     name="datavp")
{

    dfxy <- dfxy[,c(xax,yax)]
    ra <- apply(dfxy,2,range)
    di2 <- apply(ra,2,diff)
    di <- max(di2)
    whi <- which.min(di2)
    aju <- (di-di2[whi])/2
    ra[1,whi] <- ra[1,whi]-aju
    ra[2,whi] <- ra[2,whi]+aju
    ra[1,] <- ra[1,]-di*0.1
    ra[2,] <- ra[2,]+di*0.1
    if (!is.null(xlim)) {
        ra[,1] <- xlim
    }
    if (!is.null(ylim)) {
        ra[,2] <- ylim
    }

    ## Le viewport
    vp2 <- dataViewport(xscale=ra[,1], yscale=ra[,2], extension=0,
                       name=name)
    pushViewport(vp2)
    if (new)
        grid.rect(gp=gpar(fill="white"))

    ## la grille
    if (grid) {
        inte <- round(di/5)
        xgri <- c(rev(seq(0, ra[1,1], by=-inte)),
                  seq(0, ra[2,1], by=inte))
        ygri <- c(rev(seq(0, ra[1,2], by=-inte)),
                  seq(0, ra[2,2], by=inte))
        tmp <- lapply(1:length(xgri), function(i) {
            grid.lines(unit(c(xgri[i],xgri[i]),"native"),
                       unit(c(ra[1,2], ra[2,2]), "native"),
                       gp=gpar(col="grey"),
                       name=paste("grid.x.",i,sep=""))
        })
        tmp <- lapply(1:length(ygri), function(i) {
            grid.lines(unit(c(ra[1,1], ra[2,1]),"native"),
                       unit(c(ygri[i],ygri[i]), "native"),
                       gp=gpar(col="grey"),
                       name=paste("grid.y.",i,sep=""))
        })
        grid.rect()
    }

    ## les axes:
    if (axes) {
        if (ra[1,2]<0&ra[2,2]>0) {
            grid.lines(unit(c(ra[1,1], ra[2,1]),"native"),
                       unit(c(0,0), "native"), name="xaxis")
        }
        if (ra[1,1]<0&ra[2,1]>0) {
            grid.lines(unit(c(0,0), "native"),
                       unit(c(ra[1,2], ra[2,2]),"native"),
                       name="yaxis")
        }
    }

    ## Les points sur le graphique
    ma <- max(abs(value))
    if (any(value>=0)) {
        vap <- value[value>=0]
        grid.points(dfxy[value>=0,1], dfxy[value>=0,2], pch=16,
                    name="pointslabp",
                    size = unit(cexmax*vap/ma, "char"), gp=gpar(fill="black"))
    }
    if (any(value<0)) {
        van <- value[value<0]
        grid.points(dfxy[value<0,1], dfxy[value<0,2], pch=22,
                    name="pointslabn",
                    size = unit(cexmax*abs(van/ma), "char"), gp=gpar(fill="white"))
    }

    if (!is.null(biplot)) {
        s.arrowg(biplot, xax, yax, new=FALSE, summarize, axes=FALSE, grid=FALSE, xlim=ra[,1],
                 ylim=ra[,2])
    }

    return(invisible(vp2))
}






pp.holes <- function(x, p)
{
    ## critere utilise par Cook et Swayne 2007 (p. 30) pour identifier
    ## des projections avec peu de points au centre du graphe.
    ## Critere assez rapide, mais peut louper une structure quand
    ## la separation entre groupes ne se produit pas a l'origine
    ## le cas du jeu de donnees deug$tab est un bon exemple
    ##
    ## p est la dimension de l'espace dans lequel on maximise
    ## le critere (nombre de colonnes de df, passe a prpu)
    ##
    (1-sum(exp(-(x^2)/2))/length(x))/(1-exp(-p/2))
}

pp.centralmass <- function(x, p)
{
    ## critere utilise par Cook et Swayne 2007 (p. 30) pour identifier
    ## des projections avec plein de points au centre du graphe.
    ## Critere assez rapide, et utile pour identifier des outliers.
    ##
    ## p est la dimension de l'espace dans lequel on maximise
    ## le critere (nombre de colonnes de df, passe a prpu)
    ##
    (sum(exp(-(x^2)/2)) - exp(-p/2))/(1-exp(-p/2))
}


pp.FriedmanTukey <- function(x, R, prop=0.1, f=expression(R-r))
{
    ## Critere originel de Friedman et Tukey (1974).
    ## Constitue de deux parties: (i) une variance globale mesuree
    ## par une "trimmed" standard deviation (un ecart-type calcule sur
    ## les donnees auxquelles on a vire une certaine proportion des points,
    ## proportion controlee par prop)
    ## une mesure de densite de points "locale": on calcule
    ## la distance des points deux a deux, et on ne considere que les
    ## distances inferieures a R. On utilise alors, pour ces distances,
    ## la somme des fonctions monotone des differences entre R et r.
    ## on peut aussi definir d'autres fonctions.
    ##
    ## Assez lent
    ##
    z <- sort(x)
    nok <- round(prop*length(z)/2)
    z <- z[(1+nok):(length(z)-nok)]
    sk <- sd(z)
    r <- c(outer(x,x,function(x,y) abs(x-y)))
    Iii <- (R-r)>0
    dk <- sum(eval(f)*Iii)
    return(sk*dk)
}


.legpol <- function(x, J)
{
    ## polynomes de Legendre, utilises pour pp.Friedman1987
    pol <- c(1,x)
    for (j in (2:J)) {
        pol[j+1] <- ((2*j-1)*x*pol[j] - (j-1)*pol[j-1])/j
    }
    return(pol)
}

pp.Friedman1987 <- function(x, J)
{
    ## critere de Friedman (1987). Met en exergue les ecarts a la normalite
    ## se produisant dans le coeur de la distribution et non dans les queues.
    ## En theorie parce qu'en pratique, est assez sensible aux
    ## outliers, comme le dit Ripley.
    ## Redoutablement efficace, mais d'une lenteur affligeante...
    ## n'a de sens que quand l'argument "sphere" est TRUE dans
    ## prpu. J est le degre max des polynomes utilises pour
    ## approcher la distribution
    z <- 2*pnorm(x) - 1
    lp <- do.call("rbind", lapply(1:length(z), function(i) {
        y <- z[i]
        return(.legpol(y, J))
    }))
    mo <- apply(lp,2,mean)^2
    return(sum(sapply(1:J, function(j) (2*j+1)*mo[j+1]))/2)
}



posse1995 <- function(df, foo, par=NULL, centsca =TRUE, sphere=FALSE, c=1,
                      tol=1e-8, maxhalf=50, maxit=100000, nrep=500,
                      violent = FALSE, triesIfViolent=10)
{
    if (!inherits(df, "data.frame"))
        stop("df should be a data.frame")
    art <- names(as.list(args(foo)))
    if ((!all(art[1:2] == c("x", "par"))) | (length(art) != 3))
        stop("The function foo should have only two arguments named x and par")

    pc <- dudi.pca(df, scannf = FALSE, nf = ncol(df))
    if (sphere) {
        tab <- pc$l1
    } else {
        if (centsca) {
            tab <- pc$tab
        } else {
            tab <- df
        }
    }
    e1 <- new.env()

    if (!violent) {
        resu <- .Call("trouveplan", tab, quote(foo(x,par)), par, e1,
                      c, maxit, tol, nrep, maxhalf, PACKAGE="oncfsMisc")
    } else {
        resu <- .Call("trouveplan", tab, quote(foo(x,par)), par, e1,
                      c, maxit, tol, nrep, maxhalf, PACKAGE="oncfsMisc")
        for (i in 1:triesIfViolent) {
            cat("try number: ",i,"\r")
            resub <- .Call("trouveplan", tab, quote(foo(x,par)), par, e1,
                           c, maxit, tol, nrep, maxhalf)
            if (resub[[2]][[1]]>resu[[2]][[1]]) {
                resu <- resub
            }
        }
    }
    l1 <- as.data.frame(as.matrix(tab)%*%do.call("cbind",resu[[1]]))
    res <- list()
    res$tab <- tab
    c1 <- do.call("cbind",resu[[1]])
    if (sphere) {
        lambmud <- diag(1/sqrt(pc$eig))
        c1 <- as.matrix(pc$c1) %*% lambmud %*%c1
    }
    res$c1 <- as.data.frame(c1)
    names(res$c1) <- c("Axis1","Axis2")
    row.names(res$c1) <- names(df)
    res$li <- l1
    names(res$li) <- c("Axis1","Axis2")
    res$cor <- as.data.frame(cor(df, res$li))
    names(res$cor) <- c("Axis1","Axis2")
    row.names(res$cor) <- names(df)
    res$criterion <- foo
    res$par <- par
    res$call <- match.call()
    res$sphere <- sphere
    res$maxit <- maxit
    res$randtest <- as.randtest(resu[[2]][[2]], resu[[2]][[1]])
    res$val.crit <- resu[[2]][[1]]
    class(res) <- c("oprpu", "ppdim2")
    return(res)
}






pp2.holes <- function(x, p)
{
    ## critere utilise par Cook et Swayne 2007 (p. 30) pour identifier
    ## des projections avec peu de points au centre du graphe.
    ## Critere assez rapide, mais peut louper une structure quand
    ## la separation entre groupes ne se produit pas a l'origine
    ## le cas du jeu de donnees deug$tab est un bon exemple
    ##
    ## p est la dimension de l'espace dans lequel on maximise
    ## le critere (nombre de colonnes de df, passe a prpu)
    ##
    (1-sum(exp(-(sum(apply(do.call(cbind,x),1,
                           function(y) sum(y^2))))/2))/length(x[[1]]))/(1-exp(-p/2))
}


pp2.centralmass <- function(x, p)
{
    ## critere utilise par Cook et Swayne 2007 (p. 30) pour identifier
    ## des projections avec plein de points au centre du graphe.
    ## Critere assez rapide, et utile pour identifier des outliers.
    ##
    ## p est la dimension de l'espace dans lequel on maximise
    ## le critere (nombre de colonnes de df, passe a prpu)
    ##
    x <- do.call(cbind,x)
    (sum(exp(-(sum(apply(x,1,function(y) sum(y^2))))/2)) -
     exp(-p/2))/(1-exp(-p/2))
}

pp2.FriedmanTukey <- function(x)
{
    x <- do.call(cbind,x)
    r <- unlist(lapply(1:nrow(x), function(i) {
        sapply(1:nrow(x), function(j) {
            sum((x[i,]-x[j,])^2)
        })}))
    R <- (2.29*(nrow(x)^(-1/5)))^2
    Iii <- (R-r)>0
    dk <- sum(((R-r)^3)*Iii)
    return(dk)
}
