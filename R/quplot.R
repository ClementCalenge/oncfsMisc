quplot <-
function(x, ylab="", main="", ...) {
    x <- na.omit(x)
    plot(c(1:length(x))/length(x), sort(x), xlim=c(0,1),
         xlab="Quantiles", ylab=ylab, main=main, ty="n")
    axis(3, labels=FALSE)
    axis(4, labels=FALSE)
    r <- range(x)
    abline(v=seq(0,1,by=0.2), h=seq(r[1],r[2], by=diff(r)/5),
           col=grey(0.9), lwd=2)
    points(c(1:length(x))/length(x), sort(x), ...)
}
