randomRotationMatrix <- function(dim = 3, method=c("A","B"))
{
    if (dim<2)
        stop("dim should be equal to 2 or more")
    method <- match.arg(method)

    if (method=="A") {
        ## Generates the matrix
        m <- matrix(rnorm(dim^2), nrow=dim)
        ## QR decomposition:
        deqr <- qr(m)
        ## The matrices Q and R
        Q <- qr.Q(deqr)
        R <- qr.R(deqr)
        ## We require the matrix R to have positive elements on the diagonal
        ## see Stewart 1980 and Ozols 2009
        lambda <- diag(diag(R)/abs(diag(R)))
        mat <- Q%*%lambda
        ## Note that if we calculate R'= lambda%*%R, Then mat%*%R' = m, as required
        ## This is the QR decomposition relying on R matrix with positive elemnts
        ## on the diagonal
    }
    if (method=="B") {
        theta <- runif(1,0,2*pi)
        b <- sample(c(-1,1), 1)
        mat <- matrix(c(cos(theta), sin(theta), -b*sin(theta), b*cos(theta)),
                      ncol=2, byrow = TRUE)

        if (dim>2) {
            for (i in (1:(dim-2))) {
                rb <- rbind(c(1, rep(0,ncol(mat))),cbind(rep(0, ncol(mat)), mat))
                v <- mvrspherunif(1, ncol(rb))
                e1 <- c(1, rep(0,length(v)-1))
                g <- e1-v
                aa <- g/sqrt(sum(g^2))
                mat <- (diag(rep(1, length(v)))-2*aa%*%t(aa))%*%rb
            }
        }
    }
    return(mat)
}
