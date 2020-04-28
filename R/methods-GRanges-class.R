##' @rdname subseqByRef
##' @export
##'
setMethod("subseqByRef", c("GRanges", "DNAStringSet"),
          function(x, ref, ...) {
    subseq(ref[seqnames(x)], start = start(x), end = end(x), ...)
})


##' @importFrom ggplot2 ggplot geom_point facet_wrap theme
##'     element_blank element_text unit
.windowPlot <- function(data, aes, vars,
                        ...) {
    stopifnot(("window" %in% colnames(data)))
    if (is.factor(data$window.size))
        data$window.size <- as.numeric(as.character(data$window.size))
    data <- droplevels(data)
    p <- ggplot(data, {{ aes }}) + geom_point(...)
    if (!missing(vars))
        p <- p + facet_wrap( {{ vars }}, ncol = 1, strip.position = "left")
    p <- p + theme(legend.position = "none", axis.text.x = element_blank(),
                   strip.text.y = element_text(size = 8), panel.spacing = unit(0.1, "lines"))
    p
}


##' autoplot.GRanges
##'
##' @param object GRanges object
##' @param aes aesthetic mapping
##' @param vars faceting variables
##' @param ... additional parameters to geom_point
##'
##' @importFrom ggplot2 autoplot
##' @importFrom GenomicRanges GRanges
##' @export
##'
autoplot.GRanges <- function(object, aes, vars,
                             ...) {
    data <- as.data.frame(object)
    p <- .windowPlot(data, aes, vars, ...)
    p
}

##' plot
##'
##' @description plot a GRanges object
##'
##' By assumption the GRanges is a representation of a windowed
##' analysis as obtained by \code{\link[ripr]{windowScore}}.
##'
##' @param x object to plot
##' @param ... parameters for autoplot function
##'
##' @export
##' @importFrom graphics plot
##' @importFrom GenomicRanges GRanges
##'
plot.GRanges <- function(x, ...) {
    print(autoplot(x, ...))
}


##' autoplot.GRangesList
##'
##' @param object GRanges object
##' @param aes aesthetic mapping
##' @param vars faceting variables
##' @param ... additional parameters to geom_point
##' @param .id column name for bind_rows. If the input is a named
##'     GRangesList, the names will be added as an additional column
##'
##' @importFrom ggplot2 autoplot
##' @importFrom dplyr bind_rows
##' @importFrom GenomicRanges GRangesList
##'
##' @export
##'
autoplot.GRangesList <- function(object, aes, vars,
                                 ..., .id="id") {
    data <- dplyr::bind_rows(lapply(object, as.data.frame), .id=.id)
    p <- .windowPlot(data, aes, vars, ...)
    p
}

##' plot
##'
##' @description plot a GRangesList object
##'
##' By assumption the GRangesList contains GRanges objects that are
##' representations of windowed analyses as obtained by
##' \code{\link[ripr]{windowScore}}.
##'
##' @param x object to plot
##' @param ... parameters for autoplot function
##'
##' @export
##' @importFrom graphics plot
##' @importFrom GenomicRanges GRangesList
##'
plot.GRangesList <- function(x, ...) {
    print(autoplot(x, ...))
}
