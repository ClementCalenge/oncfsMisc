poissonnessPlot <- function(x, truncate=-1, outer=1.5,...)
{
    ## Troncation
    z <- round(na.omit(x))
    z <- z[z>truncate]
    nk <- unclass(table(z))
    k <- sort(unique(z))
    nko <- nk

    ## improvement
    ph <- nk/sum(nk)
    nk <- nk-0.67-0.8*ph
    if (any(nko==1))
        nk[nko==1] <- 1/exp(1)

    hk <- rep(0,length(nk))
    hk[nko>1] <- 1.96*sqrt(1-ph[nko>1])/sqrt(nk[nko>1] -
                                                 (0.47 + 0.25*ph[nko>1])*
                                                     sqrt(nko[nko>1]))
    ## intervalles de confiances (deux niveaux)
    lsup <- hk
    linf <- hk
    lsup2 <- hk
    linf2 <- hk
    lsup[nko>1] <- log(nk[nko>1])+hk[nko>1]+log(factorial(k[nko>1]))-log(sum(nko))
    linf[nko>1] <- log(nk[nko>1])-hk[nko>1]+log(factorial(k[nko>1]))-log(sum(nko))
    lsup2[nko>1] <- log(nk[nko>1])+outer*hk[nko>1]+log(factorial(k[nko>1]))-log(sum(nko))
    linf2[nko>1] <- log(nk[nko>1])-outer*hk[nko>1]+log(factorial(k[nko>1]))-log(sum(nko))
    lsup[nko==1] <- -1 + 2.717 - 2.3/sum(nko) +log(factorial(k[nko==1]))- log(sum(nko))
    linf[nko==1] <- -1-2.677 +log(factorial(k[nko==1]))-log(sum(nko))
    lsup2[nko==1] <- -1 + outer*(2.717 - 2.3/sum(nko)) +log(factorial(k[nko==1]))- log(sum(nko))
    linf2[nko==1] <- -1- outer*(2.677) +log(factorial(k[nko==1]))-log(sum(nko))

    ## On autorise données manquantes et pas entières
    y <- log(unclass(nk))+log(factorial(k))-log(sum(nk))
    yp <- log(unclass(nko))+log(factorial(k))-log(sum(nk))
    ylim <- range(c(y,lsup2,linf2))

    if (any(nko==1)) {
        plot(k, y, ty="n", xlab="k", ylab="log(k!n_k)/N", ylim=ylim, ...)
        segments(k, linf, k, lsup, lwd=4)
        segments(k, linf2, k, lsup2, lwd=1)
        points(k[nko!=1], y[nko!=1], col="black", lwd = 2)
        points(k[nko==1], y[nko==1], pch=8, cex=1.5, col="black", lwd = 2)
    } else {
        plot(k, y, xlab="k", ylab="log(k!n_k/N)", ylim=ylim)
        segments(k, linf, k, lsup, lwd=4)
        segments(k, linf2, k, lsup2, lwd=1)
    }
    points(k,yp, pch=16, col="red")

    ## Calcul de la droite. Régression linéaire
    mod <- lm(y~k)
    co <- coefficients(mod)
    pred <- function(co) {
        co[1]+co[2]*k
    }

    ## On recherche les deux intervalles les plus petits
    v <- function(co) {
        predit <- pred(co)
        all(predit>=linf2&predit<=lsup2)
    }
    vi <- function(co) {
        predit <- pred(co)
        predit>=(linf2-(lsup2-linf2)/1000)
    }
    vs <- function(co) {
        predit <- pred(co)
        predit<=lsup2+(lsup2-linf2)/1000
    }

    if (v(co)) {
        coval <- co
    } else {
        di <- (lsup2-linf2)
        dis <- sort(di)
        quels <- (di < (dis[3]+dis[2])/2)
        kq <- k[quels]
        iq <- linf2[quels]
        sq <- lsup2[quels]

        pc <- function(x1,y1,x2,y2) {
            pente <- (y2-y1)/(x2-x1)
            ordori <- y1-pente*x1
            return(c(ordori, pente))
        }
        ok <- FALSE
        iqq <- iq
        sqq <- sq
        while (!ok) {
            droite1 <- pc(kq[1],iqq[1], kq[2], sqq[2])
            droite2 <- pc(kq[1],sqq[1], kq[2], iqq[2])
            droitemil <- pc(kq[1],(iqq[1]+sqq[1])/2, kq[2], (iqq[2]+sqq[2])/2)
            i1 <- vi(droite1)
            s1 <- vs(droite1)
            i2 <- vi(droite2)
            s2 <- vs(droite2)
            im <- vi(droitemil)
            sm <- vs(droitemil)
            ## cas où ça ne marchera pas:
            c1 <- any((!i1)&(!i2)) ## valeur inférieure à borne inf pour les deux extremes
            c2 <- any((!s1)&(!s2)) ## valeur supérieure à borne sup pour les deux extremes
            if ((c1+c2)>0) {
                coval <- NA  ## si un cas merdique, on sort
                ok <- TRUE
            } else {
                ## sinon, trois cas de figure
                ## la droite du milieu passe partout, on s'arrête
                if (v(droitemil)) {
                    coval <- droitemil
                    ok <- TRUE
                }
                ## La droite du milieu est en dessous des erronés à droite de kq
                ## et au dessus à gauche de kq
                if (any(k<kq[1])) {
                    deb <- sm[k<kq[1]]
                } else {
                    deb <- logical(0)
                }
                if (any(k>kq[2])) {
                    fin <- im[k>kq[2]]
                } else {
                    fin <- logical(0)
                }
                gg <- c(deb,fin)
                ## Cette droite devient la nouvelle borne inf
                if (any(gg)) {
                    iqq[1] <- (iqq[1]+sqq[1])/2
                    sqq[2] <- (iqq[2]+sqq[2])/2
                } else {
                    sqq[1] <- (iqq[1]+sqq[1])/2
                    iqq[2] <- (iqq[2]+sqq[2])/2
                }
            }
        }
    }

    ## la droite estimée passe-t-elle par tous les intervalles?
    ## Itération
    abline(co[1],co[2], col="red", lwd=2)
    if (all(!is.na(coval)))
        abline(coval[1],coval[2], col="blue", lty=2, lwd=3)
}
