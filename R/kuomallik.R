kuomallick <-
function(y, X, nsim=10000,
                      theta = "auto",
                      sd = 100, type = c("binomial","poisson","gaussian"),
                      gaussSigma2 = var(y), listevar = NULL)
{
    ## y est la variable reponse
    ## X est la matrice de variables explicatives
    ## nsim est le nombre de simulations
    ## theta est le vecteur theta original comprenant
    ## les coefficients gamma, l'intercept et les coefficients beta.
    ## Notons que quand type=3 (gaussian), les coefficients initiaux
    ## sont ceux d'une regression lineaire
    ## sd est l'ecart-type de la prior placee sur les coefficients
    ## type indique le type de modele
    ## gaussSigma2: dans le cas ou type = 3, c'est la variance residuelle
    ## listevar: une liste comprenant autant d'elements qu'il y a de variables
    ## chaque liste donnant le numero (indexe a partir de 0) des colonnes
    ## de X.
    if (!is.null(listevar)) {
        for (i in 1:length(listevar)) {
            listevar[[i]] <- as.integer(listevar[[i]]-1)
        }

        if (!all(unlist(lapply(X, is.numeric))))
            stop("When listevar is not null, X should contain only numeric variables")
        lv <- unlist(listevar)
        if (length(lv)!=ncol(X))
            stop("Either missing variables or virtual variables in listevar")
        if (!all(unlist(lv)==(0:(ncol(X)-1))))
            stop("Either missing variables or virtual variables in listevar\n(and do not forget: numbering should start at 1)")
        if (is.null(names(listevar)))
            stop("listevar should have named elements")
    } else {
        pk <- .prepkuomallick(y, X)
        X <- pk$X
        listevar <- pk$listevar
    }

    type <- match.arg(type)
    type <- which(c("binomial","poisson", "gaussian")==type)

    if ((length(theta)==1)&&(as.character(theta)=="auto")) {
        theta <- c(rep(1,length(listevar)), rnorm(ncol(X)+1, 0, sd))
    }

    if (length(theta)!=(ncol(X)+length(listevar)+1))
        stop("incorrect length for theta")
    if (type==3) {
        theta <- c(rep(1,length(listevar)),
                   coefficients(lm(y~., data=X)), gaussSigma2)
    }
    km <- .Call("kuomallik", y, X, nsim, theta, sd, type, listevar, PACKAGE="oncfsMisc")

    km <- as.data.frame(matrix(km, nrow=nsim))
    if (type!=3) {
        names(km) <- c(paste("Effect", names(listevar), sep="."), "Intercept", names(X))
    } else {
        names(km) <- c(paste("Effect", names(listevar), sep="."), "Intercept", names(X), "MMsigma2")
    }

    for (i in 1:length(listevar)) {
        listevar[[i]] <- listevar[[i]]+1
    }

    km <- list(resu=km, y=y, Xdesign=X, thetaori = theta, sd=sd, burnin=0,
               listevar=listevar)
    class(km) <- "kuomal"
    return(km)
}




summary.kuomal <-
function(object, ...)
{
    prob <- apply(object$resu[,1:length(object$listevar)], 2, mean)
    cat("### Credible intervals and median\n")
    med <- apply(object$resu[,(length(object$listevar)+1):(ncol(object$resu))],2,median)
    q2.5 <- apply(object$resu[,(length(object$listevar)+1):(ncol(object$resu))],2,quantile, 0.025)
    q97.5 <- apply(object$resu[,(length(object$listevar)+1):(ncol(object$resu))],2,quantile, 0.975)
    df <- data.frame(median=med, q2.5=q2.5, q97.5=q97.5)
    print(df, digits=4)
    cat("\n\n### Probability\n")
    print(prob)
    cat("\n\n### Best models\n")
    uu <- .associations.kuo(object)
    print(head(uu))
    return(invisible(list(summary=df, best.mod=uu)))
}


print.kuomal <-
function(x, ...)
{
    if (!inherits(x, "kuomal"))
        stop("x should be of class kuomal")
    summary(x)
}


removeBurnin <- function(x, n=500)
{
    x$resu <- x$resu[-c(1:n),]
    x$burnin <- x$burnin+n
    return(x)
}


.associations.kuo <- function(x)
{
    ii <- table(apply(x$resu, 1, function(y) paste(y[1:length(x$listevar)],
                                                   collapse="")))
    na <- names(ii)
    na <- unlist(lapply(strsplit(na, ""), function(y) paste(names(x$listevar)[y=="1"], collapse="-")))
    na <- na[order(ii, decreasing=TRUE)]
    ii <- ii[order(ii, decreasing=TRUE)]
    names(ii) <- na
    return(ii/nrow(x$resu))
}

plot.kuomal <- function(x, ...)
{
    par(mfrow=n2mfrow(ncol(x$Xdesign)+2), mar=c(0,0,2,0))
    tmp <- lapply((length(x$listevar)+1):(ncol(x$resu)), function(i) {
        plot(1:nrow(x$resu), x$resu[,i], ty="l", axes=FALSE, main=names(x$resu)[i])
        box()
    })
}


evalconv.effects <- function(x, ngroups = 5)
{
    rr <- nrow(x$resu)
    re <- round(rr/ngroups)
    uu <- rep(1:ngroups, each=re)
    uu <- c(uu, rep(ngroups, rr - length(uu)))
    sp <- split(x$resu[,1:length(x$listevar)], factor(uu))
    do.call("rbind",lapply(sp, function(x) apply(x,2,mean)))
}

.prepkuomallick <- function(y, X)
{
    Xmod <- as.data.frame(model.matrix(lm(y~., data=X))[,-1])
    listevar <- lapply(names(X), function(x) grep(x, names(Xmod))-1)

    lv <- unlist(listevar)
    if (length(lv)!=ncol(Xmod))
        stop("Non convenient column names (grep fails)")
    if (!all(unlist(lv)==(0:(ncol(Xmod)-1))))
        stop("Non convenient column names (grep fails)")
    listevar <- lapply(listevar, function(x) as.integer(x))
    names(listevar) <- names(X)
    return(list(X=Xmod, listevar=listevar))
}
