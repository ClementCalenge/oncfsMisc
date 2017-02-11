## An environment useful to store the options
.oncfsEnv <- new.env()

.onLoad <- function(lib, pkg)
{
    environment(.oncfsEnv) <- asNamespace("oncfsMisc")
    assign(".oncoptions", list(runha=TRUE),
           envir=.oncfsEnv)

}




oncoptions <- function(...)
{
    olde <- get(".oncoptions", envir=.oncfsEnv)
    class(olde) <- "optonc"
    oo <- list(...)
    if (length(oo)>0) {
        if (is.list(oo[[1]])) {
            if (inherits(oo[[1]], "optonc"))
                oo <- oo[[1]]
        }
        newe <- olde
        for (i in names(oo)) {
            newe[[i]] <- oo[[i]]
        }
        assign(".oncoptions", newe, envir=.oncfsEnv)
        invisible(olde)
    } else {
        return(olde)
    }
}
