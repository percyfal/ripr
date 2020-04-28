##' @rdname subseqByRef
##' @export
##'
setMethod("subseqByRef", c("GRanges", "DNAStringSet"),
          function(x, ref, ...) {
    subseq(ref[seqnames(x)], start = start(x), end = end(x), ...)
})


##' @importFrom ggplot2 ggplot geom_point
.windowPlot <- function(data, aes, vars,
                        ...) {
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
##'
##' @importFrom ggplot2 autoplot
##'
##' @export
##'
autoplot.GRanges <- function(object, aes, vars,
                             ...) {
    stopifnot(("window" %in% names(mcols(object))))
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
##' @export
##' @importFrom graphics plot
##'
plot.GRanges <- function(x, ...) {
    print(autoplot(x, ...))
}


##' autoplot.GRangesList
##'
##' @importFrom ggplot autoplot
##' @importFrom dplyr bind_rows
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
##' @export
##' @importFrom graphics plot
##'
plot.GRangesList <- function(x, ...) {
    print(autoplot(x, ...))
}

