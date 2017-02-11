lss <-
function(cls="data.frame", ...) {
    ## Permet de trouver des objets heritant d'une certaine classe
    ls(envir=.GlobalEnv, ...)[
       sapply(cls, function(y)
              sapply(ls(envir=.GlobalEnv), function(x) inherits(get(x), y)))]
}
