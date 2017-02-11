mvrspherunif <- function(n=1, dim)
{
    sapply(1:n, function(j) {
        ta <- lapply(1:dim, function(i) i)
        aa <- .Call("initrand", ta, PACKAGE="oncfsMisc")
        return(aa)
    })
}
