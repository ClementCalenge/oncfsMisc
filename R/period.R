period <- function(yt, ncons=1L)
{
    n <- length(yt)
    yt <- (yt - mean(yt))/sd(yt)
    ct <- do.call("data.frame",
                  lapply(c(1:floor(n/2)),
                         function(o) cos(2*pi * o/n * (1:n))))
    st <- do.call("data.frame",
                  lapply(c(1:floor(n/2)),
                         function(o) sin(2 * pi * o/n * (1:n))))
    cct <- apply(ct, 2, function(x) sum(x * yt)^2)
    cst <- apply(st, 2, function(x) sum(x * yt)^2)
    plot(1:floor(n/2), (cct + cst)/n, ty = "l",
         xlab="Number of cycles during the duration of the study",
         ylab="Periodogram")
    spe <- (cct + cst)/n
    fre <- c(1:floor(n/2))
    names(spe) <- NULL
    cyc <- fre[order(spe, decreasing=TRUE)][1:ncons]
    spb <- spe[order(spe, decreasing=TRUE)][1:ncons]
    cat("The number of cycles during the duration of the study is:\n")
    print(data.frame(frequency=cyc, Peri=spb))
    res <- data.frame(fre,spe)
    res <- res[order(res[,2], decreasing=TRUE),]
    attr(res,"yt") <- yt
    attr(res,"ncons") <- ncons
    class(res) <- "period"
    return(invisible(res))
}


fitmod <- function(per, increaseres=1L)
{
    if (!inherits(per, "period"))
        stop("per should be of class period")
    increaseres <- round(increaseres)
    if (increaseres < 0.5)
        stop("increaseres should be strictly greater than 0")

    yt <- attr(per,"yt")
    nc <- attr(per,"ncons")
    n <- length(yt)
    df <- do.call("data.frame", lapply(1:nc, function(i) {
        cyc <- per$fre[i]
        st <- sin((2*pi*cyc/n)*c(1:n))
        ct <- cos((2*pi*cyc/n)*c(1:n))
        re <- data.frame(ct,st)
        names(re) <- paste(c("ct","st"),i,sep="")
        return(re)
    }))
    df$yt <- yt
    mod <- lm(yt~.-1, data=df)
    dfpr <- do.call("data.frame", lapply(1:nc, function(i) {
        cyc <- per$fre[i]
        st <- sin((2*pi*cyc/n)*c(1:(n*increaseres))/increaseres)
        ct <- cos((2*pi*cyc/n)*c(1:(n*increaseres))/increaseres)
        re <- data.frame(ct,st)
        names(re) <- paste(c("ct","st"),i,sep="")
        return(re)
    }))
    return(data.frame(t=c(1:(n*increaseres))/increaseres, prediction=predict(mod, newdata=dfpr)))
}
