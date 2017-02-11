.ha12 <-
function()
{
        par(mar=c(0,0,0,0))
        symbols(0,0,circles=1.1, inches=FALSE, lwd=3)
        symbols(0,0,circles=1, inches=FALSE, lwd=3, add=TRUE)
        heures <- c(3:1,12:4)
        text(0.9*cos(seq(0,2*pi, by=2*pi/12)), 0.9*sin(seq(0,2*pi, by=2*pi/12)),
             heures, cex=2)
        hourc<-as.numeric(substr(date(),12,13))%%12
        minc<-as.numeric(substr(date(),15,16))
        secc<-as.numeric(substr(date(),18,19))
        angh <- (pi/2) - ((hourc+minc/60)/12)*2*pi
        segments(0,0, 0.4*cos(angh), 0.5*sin(angh), lwd=10)
        angm <- (pi/2) - ((minc/60))*2*pi
        segments(0,0, 0.8*cos(angm), 0.8*sin(angm), lwd=10)
        angs <- (pi/2) - ((secc/60))*2*pi
        segments(0,0, 0.8*cos(angs), 0.8*sin(angs), lwd=4, col="red")
        points(0,0, cex=5, pch=16, col="blue")
}

startclock <- function()
{
    oncoptions(runha=TRUE)
    return(invisible(NULL))
}

clock <-
function() {
    .ha12()
    if (oncoptions()$runha)
        tcltk::tcl("after", 1000, clock)
}

stopclock <- function()
{
    oncoptions(runha=FALSE)
    graphics.off()
}
