boxplot.pca <-
function(x, ...)
{
    if (!inherits(x,"pca"))
        stop("x should be of class pca")
    y <- c(as.matrix(x$tab))
    x2 <- factor(rep(names(x$tab),each=nrow(x$tab)))
    opar <- par(mar=c(5,5,4,2))
    boxplot(y~x2, horizontal=TRUE, las=1, col="grey", ...)
    par(opar)
}
