calculeR2 <- function(P, hc)
{
    RSQ <- rep(0, nrow(P))
    SQTot <- sum(scale(P, scale = FALSE)^2)
    for (i in 1:nrow(P)) {
        Cla <- as.factor(cutree(hc, i))
        RSQ[i] <- sum(t((t(sapply(1:ncol(P),
                                  function(i)
                            tapply(P[,i], Cla, mean)))-
                         apply(P, 2, mean))^2) * as.vector(table(Cla)))/SQTot
    }
    return(RSQ)
}

exploreClust <- function(ds, P=NULL, meth=c("clust","R2","break","cut"), n=NULL,  hang=-1, ...)
{
    meth <- match.arg(meth)
    if (is.null(P)) {
        if (any(c("cut","R2")==meth))
            stop("meth cannot be cut or R2 if P is null")
        if ((ncol(P)>2)&(meth=="cut")) {
            warning("only the first two columns of P used for the plot of cut")
            P <- P[,1:2]
        }
    }

    opar <- par(mfrow=c(3,2))
    on.exit(par(opar))


    if (meth=="clust") {
        plot(hclust(ds, method="single"), hang=hang, main="single", ...)
        plot(hclust(ds, method="complete"),hang=hang, main="complete", ...)
        plot(hclust(ds, method="average"), hang=hang, main="average", ...)
        plot(hclust(ds, method="mcquitty"), hang=hang, main="wpgma", ...)
        plot(hclust(ds, method="centroid"), hang=hang, main="centroid", ...)
        plot(hclust(ds, method="ward.D2"), hang=hang, main="ward", ...)
    }


    if (meth=="cut") {
        if (is.null(n))
            stop("n is required when method == \"cut\"")
        plot(P, col=(cutree(hclust(ds, method="single"), n)), main="single", ...)
        plot(P, col=(cutree(hclust(ds, method="complete"), n)), main="complete", ...)
        plot(P, col=(cutree(hclust(ds, method="average"), n)), main="average", ...)
        plot(P, col=(cutree(hclust(ds, method="mcquitty"), n)), main="wpgma", ...)
        plot(P, col=(cutree(hclust(ds, method="centroid"), n)), main="centroid", ...)
        plot(P, col=(cutree(hclust(ds, method="ward.D2"), n)), main="ward", ...)
    }

    if (meth=="R2") {
        plot(1:nrow(P), calculeR2(P, hclust(ds, method="single")),
             ty="b", main="single", xlab="Nb of clusters", ylab="R2", ...)
        plot(1:nrow(P), calculeR2(P, hclust(ds, method="complete")),
             ty="b", main="complete", xlab="Nb of clusters", ylab="R2",...)
        plot(1:nrow(P), calculeR2(P, hclust(ds, method="average")),
             ty="b", main="average", xlab="Nb of clusters", ylab="R2",...)
        plot(1:nrow(P), calculeR2(P, hclust(ds, method="mcquitty")),
             ty="b", main="mcquitty", xlab="Nb of clusters",ylab="R2", ...)
        plot(1:nrow(P), calculeR2(P, hclust(ds, method="centroid")),
             ty="b", main="centroid", xlab="Nb of clusters", ylab="R2",...)
        plot(1:nrow(P), calculeR2(P, hclust(ds, method="ward.D2")),
             ty="b", main="ward.D2", xlab="Nb of clusters", ylab="R2",...)
    }

    if (meth=="break") {
        plot(rev(hclust(ds, method="single")$height),
             ty="b", main="single", ylab="height",...)
        plot(rev(hclust(ds, method="complete")$height),
             ty="b", main="complete", ylab="height",...)
        plot(rev(hclust(ds, method="average")$height),
             ty="b", main="average", ylab="height",...)
        plot(rev(hclust(ds, method="mcquitty")$height),
             ty="b", main="mcquitty", ylab="height",...)
        plot(rev(hclust(ds, method="centroid")$height),
             ty="b", main="centroid", ylab="height",...)
        plot(rev(hclust(ds, method="ward.D2")$height),
             ty="b", main="ward", ylab="height",...)
    }
    return(invisible(NULL))
}
