spm <- function(x, density=FALSE, ajoutQ=FALSE, add0=TRUE, avaxes=TRUE)
{
    if (!inherits(x,"data.frame"))
        stop("x should be a data.frame")

    na <- names(x)

    par(mfrow=c(ncol(x), ncol(x)))
    for (i in 1:ncol(x)) {
        for (j in 1:ncol(x)) {
            if (i==j) {
                var <- x[,i]
                if (!avaxes) {
                    par(mar=c(0.1,0.1,0.1,0.1))
                } else {
                    par(mar=c(2, 2, 0.2, 0.2))
                }
                ra <- range(var)
                plot(c(1:length(var))/length(var),
                     sort(var), axes=avaxes)
                abline(v=c(0,0.25,0.5,0.75,1), col="grey")
                if (add0)
                    abline(h=0)
                box()
                text(0.5, max(var), na[i], pos=1)
            } else {
                if (!avaxes) {
                    par(mar=c(0.1,0.1,0.1,0.1))
                } else {
                    par(mar=c(2, 2, 0.2, 0.2))
                }
                plot(x[,j], x[,i],
                     axes=avaxes)
                qu1 <- quantile(x[,j])
                if (ajoutQ)
                    abline(v=qu1, col="grey")
                qu2 <- quantile(x[,i])
                if (ajoutQ)
                    abline(h=qu2, col="grey")

                if (density) {
                    k <- MASS::kde2d(x[,j],x[,i])
                    cnt <- grDevices::contourLines(k$x, k$y, k$z)
                    n <- length(cnt)
                    cols <- rev(colorspace::sequential_hcl(n))
                    for( k in seq_len(n) ) lines(cnt[[k]], col=cols[k], lwd=2)
                }
                if (add0)
                    abline(h=0, v=0)

                box()

            }
        }
    }

}
