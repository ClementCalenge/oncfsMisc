invSAT <- function(X)
{
    n <- length(X)
    eta <- 1:(n-1)
    ev <- (n%%2==0)
    if (ev) {
        for (i in 1:(n/2)) {
            pr <- atan2(X[2*i], X[2*i-1])
            eta[2*i-1] <- ifelse(pr<0, pr+2*pi, pr)/(2*pi)
        }
        z <- 1
        for (i in 1:((n/2)-1)) {
            if (i>1)
                z <- z*(eta[2*(i-1)]^(1/(n-2*(i-1))))
            eta[2*i] <- (sin(acos(X[2*i]/(z*(sin(2*pi*eta[2*i-1]))))))^(n-2*i)
        }
        return(eta)
    } else {
        ## for all until n-3
        for (i in 1:((n-3)/2)) {
            pr <- atan2(X[2*i], X[2*i-1])
            eta[2*i-1] <- ifelse(pr<0, pr+2*pi, pr)/(2*pi)
        }
        z <- 1
        for (i in 1:(((n-3)/2))) {
            if (i>1)
                z <- z*(eta[2*(i-1)]^(1/(n-2*(i-1))))
            eta[2*i] <- (sin(acos(X[2*i]/(z*(sin(2*pi*eta[2*i-1]))))))^(n-2*i)
        }
        ##
        pr <- 1
        for (j in 1:((n-3)/2)) {
            pr <- pr*eta[2*j]^(1/(n-2*j))
        }
        eta[n-2] <- (1 - X[n]/(pr))/2
        uu <- pr*sqrt((4*eta[n-2])-(4*eta[n-2]^2))
        uu <- atan2(X[n-1]/uu, X[n-2]/uu)
        uu <- ifelse(uu<0, uu+2*pi, uu)
        eta[n-1] <- uu/(2.0*pi)
        return(eta)

    }
}



SAT <- function(eta)
{
    n <- length(eta)+1
    X <- 1:n

    ev <- (n%%2==0)
    if (ev) {

        for (i in 1:(n/2-1)) {
            pr <- 1
            if (i>1) {
                for (j in 1:(i-1))
                    pr <- pr*eta[2*j]^(1/(n-2*j))
            }
            pr <- pr*cos(asin(eta[2*i]^(1/(n-2*i))))*sin(2*pi*eta[2*i-1])
            X[2*i] <- pr
        }
        pr <- 1
        for (j in 1:(n/2-1))
            pr <- pr*eta[2*j]^(1/(n-2*j))
        X[n] <- pr*sin(2*pi*eta[n-1])
        for (i in 1:(n/2)) {
            X[2*i-1] <- X[2*i]/tan(2*pi*eta[2*i-1])
        }
        return(X)

    } else {
        pr <- 1
        for (j in 1:((n-3)/2)) {
            pr <- pr*eta[2*j]^(1/(n-2*j))
        }
        X[n-2] <- pr*sqrt(4*eta[n-2]-(4*eta[n-2]^2))*cos(2*pi*eta[n-1])
        X[n-1] <- pr*sqrt(4*eta[n-2]-(4*eta[n-2]^2))*sin(2*pi*eta[n-1])
        X[n] <- pr*(1-2*eta[n-2])
        for (i in 1:((n-3)/2)) {
            pr <- 1
            if (i>1) {
                for (j in 1:(i-1)) {
                    pr <- pr*eta[2*j]^(1/(n-2*j))
                }
            }
            X[2*i] <- pr*cos(asin(eta[2*i]^(1/(n-2*i))))*sin(2*pi*eta[2*i-1])
        }
        for (i in 1:((n-3)/2)) {
            X[2*i-1] <- X[2*i]/tan(2*pi*eta[2*i-1])
        }
        return(X)
    }
}
