examACP <- function(pc, ncolmax=3, niter=100)
{
    if (!inherits(pc, "pca"))
        stop("pc should inherit the class \"pca\"")

    pc <- ade4::redo.dudi(pc, newnf=length(pc$eig))
    l1 <- as.matrix(pc$l1)
    c1 <- as.matrix(pc$c1)
    eig <- diag(sqrt(pc$eig))
    res <- t(sapply(1:niter, function(r) {
        veceig <- sapply(1:ncolmax, function(j) {
            ## Reconstruction du tableau (une partie structurelle+une partie residuelle)
            ## stru <- l1[,1:j]%*%eig[1:j,1:j]%*%t(c1[,1:j])
            ## Mais la partie structurelle, on s'en fiche, ce qui nous interesse
            ## est la partie residuelle:
            res <- l1[,j:ncol(l1)]%*%eig[j:ncol(l1),j:ncol(l1)]%*%t(c1[,j:ncol(l1)])
            ## On randomise alors les observations des residus, pour "casser" les correlations
            res <- as.data.frame(do.call("cbind", lapply(1:ncol(res), function(k) sample(res[,k]))))
            ## et on fait l'ACP sur le tableau des residus (non centree -- le tableau l'est
            ## deja -- et surtout non reduite)
            eig2 <- ade4::dudi.pca(res, scannf=FALSE, scale=FALSE, center=FALSE)$eig[1]
            ## La premiere valeur propre correspond a la jeme de l'acp
            return(eig2)
        })
        return(veceig)
    }))

    ## On met le truc en forme
    plot(1:ncolmax, pc$eig[1:ncolmax], xlab="Axis number", ylab="Eigenvalue", ty="n",
         ylim=c(0,max(pc$eig)))

    pollig <- function(r, mat, percent)
    {
        q1 <- (100-percent)/200
        q2 <- 1-q1
        pola <- t(sapply(1:ncol(mat), function(i) {
            quantile(unlist(mat[,i]), prob=c(q1,q2))
        }))
        rbind(cbind(r, pola[,1]), cbind(rev(r), rev(pola[,2])))
    }
    ## On affiche les randomisations
    polygon(pollig(1:ncolmax, res, 100), col="lightgrey", border=NA)
    polygon(pollig(1:ncolmax, res, 95), col="grey", border=NA)
    polygon(pollig(1:ncolmax, res, 90), col="darkgrey", border=NA)

    ## Et on affiche les valeurs propres observees
    points(1:ncolmax, pc$eig[1:ncolmax], pch=16, ty="b")
    ## et on renvoie le resultat
    invisible(list(obs=pc$eig[1:ncolmax], sim=res))
}
