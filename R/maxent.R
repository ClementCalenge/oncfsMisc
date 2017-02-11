## Calculations of the characteristics fj from the environmental variables
prepmaxent <- function(dfo, thr=NULL, lithr=NULL, nthr=10, rare=0.1,
                       PreviousObject=NULL)
{
    if (!is.data.frame(dfo))
        stop("dfo should be a data.frame")

    if (any(sapply(dfo, is.character)))
        stop("Character string variables are not allowed in this function.\n please transform them into factor before")

    if (!is.null(PreviousObject)&((!is.null(thr))|(!is.null(lithr))))
        warning("When PreviousObject is set, other arguments will not be taken into account")

    if (!is.null(PreviousObject)) {
        liattr <- attr(PreviousObject,"ForReuse")
        if (is.null(liattr))
            stop("non convenient object for Previous object")
        thr <- liattr$thr
        lithr <- liattr$lithr
        cons <- liattr$cons
        nthr <- liattr$nthr
        rare <- liattr$rare
    } else {
        cons <- NULL
    }

    ## Transformation of the factors into disjonctive table
    fac <- sapply(dfo, is.factor)
    dfob <- dfo
    if (sum(fac)>0) {
        xf <- dfo[, fac, drop=FALSE]
        xf <- as.data.frame(do.call(cbind, lapply(1:ncol(xf), function(i) {
            x <- xf[,i]
            ddd <- data.frame(x)
            names(ddd) <- names(xf)[i]
            ca <- acm.disjonctif(ddd)
            names(ca) <- gsub("\\.","_",names(ca))
            return(ca)
        })))
        dfob <- cbind(dfob[,!fac],xf)
    }

    ## Definition of the interaction and squared variables (including the factors)
    cb <- as.data.frame(cbind(dfob, do.call("cbind", lapply(1:ncol(dfob), function(i) {
        do.call("cbind", lapply(1:ncol(dfob), function(j) {
            dfob[,i]*dfob[,j]
        }))
    }))))
    na <- c(names(dfob), unlist(lapply(1:ncol(dfob), function(i) {
        unlist(lapply(1:ncol(dfob), function(j) {
            if (i!=j) {
                return(paste(names(dfob)[i],".",names(dfob)[j], sep=""))
            }
            return(paste(names(dfob)[i],"2", sep=""))
        }))
    })))
    names(cb) <- na

    ## If there are any threshold features to define
    if (!is.null(thr)) {

        ## List of variables for which thresholds should be defined
        xb <- dfo[,thr, drop=FALSE]

        if (any(!sapply(xb, is.numeric)))
            stop("Threshold variables can only be defined from numeric variables")

        ## if no list of threshold is passed as arguments, generates
        ## one (based on the number of thresholds passed as arguments)
        if (is.null(lithr)) {
            lithr <- lapply(1:ncol(xb), function(i) {
                a <- xb[,i]
                sq <- quantile(a, seq(0,1,length=nthr))
            })
        }

        ## then, for each variable to be transformed into threshold, builds
        ## a matrix of dummy variables storing these threshold variables
        xthr <- do.call(cbind, lapply(1:ncol(xb), function(i) {
            a <- xb[,i]
            sq <- lithr[[i]]
            sq[1] <- -1e12
            sq[length(sq)] <- 1e12
            ddd <- data.frame(cut(a, sq))
            names(ddd) <- paste(names(xb)[i], "_thr",sep="")
            ca <- acm.disjonctif(ddd)
            names(ca) <- gsub("\\.","_",names(ca))
            return(ca)
        }))

        ## And pastes it to the matrix of other characteristics
        cb <- cbind(cb, xthr)
    }

    ## Which variables in cb should be kept for model fitting
    if (is.null(cons)) {
        ## By default, all
        cons <- rep(TRUE, ncol(cb))

        ## But if there are replicated columns, delete them
        for (i in 1:(ncol(cb)-1)) {
            for (j in (i+1):ncol(cb)) {
                if (all(cb[,i]==cb[,j]))
                    cons[j] <- FALSE
            }
        }

        ## If there are columns chacterized by a null variance, delete them
        cons[apply(cb,2,var)<1e-8] <- FALSE

        ## For binary characteristics, checks that these characteristics are not
        ## too rare as defined
        binaires <- apply(cb,2,function(x) length(unique(x))==2)
        rares <- apply(cb,2,function(x) sum(x)<(rare*length(x)))
        cons[rares&binaires] <- FALSE
    }

    ## Delete all the variables not kept
    cb <- cb[,cons]

    ## rescale the variables
    if (is.null(PreviousObject)) {
        rescale <- lapply(cb, function(x) c(min(x),max(x)))
    } else {
        rescale <- liattr$rescale
    }
    na <- names(cb)
    cb <- as.data.frame(lapply(1:length(cb), function(i) {
        as.double((cb[,i]-rescale[[i]][1])/(rescale[[i]][2]-rescale[[i]][1]))
    }))

    names(cb) <- na
    attr(cb, "ForReuse") <- list(thr=thr, lithr=lithr, cons=cons, nthr=nthr,
                                 rare=rare, rescale=rescale)

    return(cb)
}



## Function interfacing the C code allowing the fit of the maxent model per se
maxentcore <- function(df, y, beta=0, betaj = NULL, verbose=FALSE,
                       critconv=0.0000000001, maxIter=10000,
                       typealgo=c("BFGS","Sequential"))
{
    if (nrow(df)!=length(y))
        stop("Dimensions of df and y do not match")
    if (!all(sort(unique(y))==c(0,1)))
        stop("y should contain only 0 and 1")
    if (!inherits(df, "data.frame"))
        stop("df should be a data.frame")
    typealgo <- match.arg(typealgo)
    if (!all(sapply(df, is.numeric)))
        stop("df should only contain numeric variables.\n Please consider the use of prepmaxent bfore.")

    ## preliminary transformation of the explanatory variables and of the coefficients
    if (is.null(betaj)) {
        betaj <- unlist(apply(df,2,function(x) {
            beta*sqrt(var(x[y>0.5])/sum(y))
        }))
    }

    ## fits the model
    if (typealgo=="Sequential") {
        modele <- .Call("maxentSequential", df, df[y==1,], as.double(betaj),
                        as.double(y), as.double(critconv), as.integer(verbose),
                        as.integer(maxIter))
    } else {
        lambda <- rep(0, ncol(df))
        modele <- .Call("maxentBFGS", df, as.double(betaj), as.integer(maxIter),
                        as.double(critconv), as.integer(verbose), as.numeric(lambda),
                        as.double(y))
    }
    names(modele) <- c("coefficients", "number.of.iterations", "convergence")
    names(modele$coefficients) <- names(df)

    ## penalized likelihood at the solution
    vpen <- .Call("vraisemblance_penaliseer", as.double(betaj), as.double(modele$coefficients),
                  df,  as.double(y))

    modele$penalized.likelihood <- vpen[1]
    modele$likelihood <- vpen[2]
    modele$penalty <- vpen[2]-vpen[1]
    modele$data <- list(X=df, y=y)
    modele$call <- match.call()
    modele$fitted.values <- apply(as.matrix(df),1,function(x) exp(sum(x*modele$lambda)))
    modele$betaj <- betaj
    class(modele) <- "maxent"
    return(modele)
}

print.maxent <- function(x, ...)
{
    if (!inherits(x, "maxent"))
        stop("x should be an object of class \"maxent\"")
    cat("Object of class \"maxent\".\nThis object is a list with the following elements:\n")
    cat("   $coefficients: the coefficients associated to the features\n")
    cat("   $number.of.iterations: the number of iterations of the algorithm\n")
    cat("   $convergence: whether the algorithm has converged\n")
    cat("   $penalized.likelihood: the value of the penalized likelihood at the solution\n")
    cat("   $likelihood: the value of the unpenalized likelihood at the solution \n")
    cat("   $penalty: the value of the penalty at the solution\n")
    cat("   $data: the data (X is the matrix of features, and y defines the response\n")
    cat("   $call: the call to the function\n")
    cat("   $fitted.values: the fitted values\n")
    cat("   $betaj: the value of the tuning parameters for each feature\n")
}


predict.maxent <- function(object, newdata = NULL, type=c("link","response"), ...)
{
    if (!inherits(object, "maxent"))
        stop("object should be an object of class \"maxent\"")
    type <- match.arg(type)
    ## Checks that newdata
    if (is.null(newdata)) {
        newdata <- object$data$X
    } else {
        if (!inherits(newdata, "data.frame"))
            stop("newdata should be a data.frame")
        if (ncol(newdata) != length(object$coefficients))
            stop("The number of columns in newdata does not correspond\n to the number of coefficients in the model.\nAre you sure that newdata contains the *features* and not the *variables*.\nSee prepmaxent otherwise.")
        if (!all(sapply(newdata, is.numeric)))
            stop("All variables in newdata should be numeric")
    }
    aa <- as.vector(as.matrix(newdata)%*%coef(object))
    if (type=="response")
        aa <- exp(aa)
    return(aa)
}


coef.maxent <- function(object, ...)
{
    if (!inherits(object, "maxent"))
        stop("object should be of class maxent")
    return(object$coefficients)
}

