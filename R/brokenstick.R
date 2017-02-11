brokenstick <-
function(pc, niter=1000)
{
    ## Approche pourrie: le principe. On a p valeurs propres
    ## On considere le segment 0 1 que l'on decoupe aleatoirement en p sous-segments
    ## Pour ce faire, on tire au sort p-1 valeurs dans une loi uniforme (0,1).
    ## Ces p-1 valeurs decoupent le segment 0,1 en p segments.
    ##
    ## On ordonne alors les segments dans l'ordre decroissant. La taille de ces segments suit
    ## une "broken stick distribution". Cette distribution peut etre utilisee
    ## comme "benchmark" pour comparer la distribution
    ## L'esperance d'une telle loi pour le jeme segment (dans l'ordre decroissant)
    ## est (eq. 6.49 dans Legendre, voir page 410)
    ## E(l_j) = 1/p sum_(x=j)^p 1/x
    ## avec l_j la longueur du jeme segment
    ##
    ## Ici, on simule cette distribution niter fois, et on represente les distributions
    ## pour chaque axe.
    ## Legendre recommande de ne s'interesser qu'a l'esperance, sur des bases
    ## empiriques, comme toujours. Frontier (1976) fait de meme.
    ## Pas vraiment de justifications theoriques a ca.
    ## En outre, ne prend pas du tout en compte le nombre d'individus dans l'ACP
    ## Par exemple, si p=2 (deux variables), sous le broken stick, on a une chance
    ## sur 4 d'avoir au moins 75% de la variance exprimee par le premier axe. Pire,
    ## quand p = 2, E(l_1) = (1/2) * ( 1 + (1/2)) = 0.75.
    ## C'est l'esperance qui est de 75%. En d'autres termes, meme si l'on a 10000 individus
    ## dans l'ACP, avoir 70% de variations sur le premier axe ne conduira
    ## pas a affirmer que l'on conserve le premier axe...
    ncolmax <- length(pc$eig)
    res <- t(sapply(1:niter, function(i) {
        sort(diff(sort(c(0,runif(length(pc$eig)-1, 0, sum(pc$eig)), 1))),decreasing=TRUE)
    }))
    pollig <- function(r, mat, percent)
    {
        q1 <- (100-percent)/200
        q2 <- 1-q1
        pola <- t(sapply(1:ncol(mat), function(i) {
            quantile(unlist(mat[,i]), prob=c(q1,q2))
        }))
        rbind(cbind(r, pola[,1]), cbind(rev(r), rev(pola[,2])))
    }
    plot(1:ncolmax, pc$eig[1:ncolmax], xlab="Axis number", ylab="Eigenvalue", ty="n",
         ylim=c(0,max(c(pc$eig, c(res)))))
    ## On affiche les randomisations
    polygon(pollig(1:ncolmax, res, 95), col="lightgrey", border=NA)
    polygon(pollig(1:ncolmax, res, 90), col="grey", border=NA)
    polygon(pollig(1:ncolmax, res, 80), col="darkgrey", border=NA)
    ## Et on affiche les valeurs propres observees
    points(1:ncolmax, pc$eig[1:ncolmax], pch=16, ty="b")
    points(1:ncolmax, sum(pc$eig)*(1/ncolmax)*rev(cumsum(rev((1/(1:ncolmax))))), col="red",
           ty="b")
    points(1:ncolmax, apply(res,2,mean), col="green", ty="b")
    ## et on renvoie le resultat
    invisible(list(obs=pc$eig[1:ncolmax], sim=res))
}
