togrey <-
function(x)
{
    x <- -x
    grey(0.8*(x-min(x))/diff(range(x))+0.1)
}
