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
##' @param .id name of id column added by dplyr::bind_rows
##'
##' @return data.frame
##'
##' @importFrom dplyr bind_rows
##'
##' @export
##'
setMethod("as.data.frame", signature = "AlignmentPairsList",
          function(x, ..., .id="id") {
    dplyr::bind_rows(lapply(x, as.data.frame, ...), .id=.id)
})



##' autoplot.AlignmentPairsList
##'
##' @importFrom ggplot2 autoplot ggplot geom_point geom_boxplot
##'     geom_violin geom_density facet_wrap
##'
##' @param object AlignmentPairsList
##' @param aes aes mapping
##' @param vars variable mapping to facet plots
##' @param ... additional parameters to ggplot function
##' @param which which plot to make. 'grid' option makes a scatter
##'     plot with marginal densities
##'
##' @export
##'
autoplot.AlignmentPairsList <- function(object, aes, vars, ..., which="point") {
    data <- as.data.frame(object)
    p <- ggplot(data, {{ aes }})
    which <- match.arg(which, c("point", "boxplot", "violin", "density"))
    if (which == "point")
        p <- p + geom_point(...)
    else if (which == "boxplot")
        p <- p + geom_boxplot(...)
    else if (which == "violin")
        p <- p + geom_violin(...)
    else if (which == "density")
        p <- p + geom_density(...)
    if (!missing(vars))
        p <- p + facet_wrap( {{ vars }} )
    p
}


##' plot
##'
##' @description plot an AlignmentPairsList
##'
##' @param x object to plot
##' @param ... additional arguments for autoplot
##'
##' @export
##' @importFrom graphics plot
##'
plot.AlignmentPairsList <- function(x, ...) {
    print(autoplot(x, ...))
}
