##' AlignmentPairsList
##'
##' @export
##' @rdname AlignmentPairsList-class
##'
##' @importFrom methods new
##'
setMethod("AlignmentPairsList", "list",
          function(obj) new("AlignmentPairsList", listData = obj))


##' as.data.frame
##'
##' @description Convert AlignmentPairsList to data.frame.
##'
##' @param x AlignmentPairsList object
##' @param ... additional parameters to lapply
##'
##' @return data.frame
##'
##' @export
##'
setMethod("as.data.frame", signature = "AlignmentPairsList",
          function(x, ...) {
    do.call("rbind", lapply(x, as.data.frame, ...))
})



##' @importFrom ggplot2 autoplot ggplot geom_point geom_boxplot geom_violin
autoplot.AlignmentPairsList <- function(object, aes, which="point", ...) {
    data <- do.call("rbind", lapply(object, as.data.frame))
    p <- ggplot(data, {{ aes }})
    which <- match.arg(which, c("point", "boxplot", "violin"))
    if (which == "point")
        p <- p + geom_point(...)
    else if (which == "boxplot")
        p <- p + geom_boxplot(...)
    else if (which == "violin")
        p <- p + geom_violin(...)
    p
}


##' plot
##'
##' @description plot an AlignmentPairsList
##'
##' @export
##' @importFrom graphics plot
plot.AlignmentPairsList <- function(x, ...) {
    print(autoplot(x, ...))
}
