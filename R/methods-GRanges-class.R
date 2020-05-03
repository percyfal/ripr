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
    data[[.id]] <- factor(data[[.id]], levels=unique(data[[.id]]))
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


##'
##' @importFrom S4Vectors DataFrame
##' @importFrom Biostrings alphabetFrequency
##'
##' @export
##' @rdname windowGC
##'
setMethod("windowGC", c("GRanges", "DNAStringSet"),
          function(windows, ref) {
    DataFrame(gc = unlist(lapply(subseqByRef(windows, ref), function(x) {
        sum(alphabetFrequency(x)[c("G", "C")]) / length(x)
    })))
})


##'
##' @importFrom S4Vectors DataFrame
##'
##' @export
##' @rdname windowRIP
##'
setMethod("windowRIP", c("GRanges", "DNAStringSet"),
          function(windows, ref, ...) {
    dots <- list(...)
    arglist <- list(x = AlignmentItem(windows), ref = ref)
    arglist <- append(arglist, dots[which(names(dots) %in% names(formals(ripr::count)))])
    .rip <- do.call(calculateRIP, arglist)
    mcols(.rip)[, c("rip.product", "rip.substrate", "rip.composite")]
})


##'
##' @importFrom IRanges findOverlaps
##' @importFrom S4Vectors DataFrame
##'
##' @export
##' @rdname windowRepeatContent
##'
setMethod("windowRepeatContent", c("GRanges", "AlignmentItem", "DNAStringSet"),
          function(windows, obj, ref) {
    .res = DataFrame(repeat.obs = rep(0, length(windows)))
    hits <- findOverlaps(windows, obj, ignore.strand = TRUE)
    obs <- c(by(as.data.frame(hits), as.data.frame(hits)$queryHits,
                function(h){
        sum(width(intersect(ranges(obj[h$subjectHits]),
                            reduce(ranges(windows[h$queryHits])))))}))
    .res$repeat.obs[as.integer(names(obs))] <- obs
    .res
})

##'
##' @rdname windowScore
##' @export
##'
##' @param lambda list of functions that can be applied to the tuple
##'     x, ref, and the generated windows
##'
##'
setMethod("windowScore", c("GRanges", "AlignmentPairs", "DNAStringSet"),
          function(windows, x, ref, which = c("rip", "repeat.content", "gc"), ...,
                   lambda = list()) {
    dots <- list(...)
    which <- match.arg(which, c("rip", "repeat.content", "gc.content"), several.ok = TRUE)
    if ("rip" %in% which) {
        message("Calculating rip scores")
        mcols(windows) <- cbind(mcols(windows), windowRIP(windows, ref))
    }
    if ("repeat.content" %in% which) {
        message("Calculating repeat content")
        mcols(windows) <- cbind(mcols(windows), windowRepeatContent(windows, query(x), ref))
    }
    if ("gc.content" %in% which) {
        message("Calculating gc content")
        mcols(windows) <- cbind(mcols(windows), windowGC(windows, ref))
    }
    for (f in names(lambda)) {
        message("Applying ", f)
        windows <- lambda[[f]](x, ref, windows)
    }
    windows
})


##' expectedWindowScores
##'
##' @description Calculated expected window scores for a windowed value
##'
##' @param windows ranges on which to operate
##' @param col column name containing values in windows
##' @param ... additional parameters
##'
##' @export
##' @rdname expectedWindowScores
##'
##' @return vector of expected values
##'
expectedWindowScores <- function(windows, col, ...) UseMethod("expectedWindowScores", windows)



##'
##' @importFrom GenomeInfoDb seqinfo seqnames
##' @importFrom BiocGenerics width
##'
##' @export
##' @rdname expectedWindowScores
##'
##' @param seqinfo use seqinfo to define sequence lengths by which to
##'     normalize scores when calculating expected values
##' @param by.seqnames normalize scores by sequence name, e.g. within
##'     chromosome instead of across entire genome
##'
expectedWindowScores.GRanges <- function(windows, col, seqinfo=NULL, by.seqnames=FALSE, ...) {
    ## FIXME: not using seqinfo at the moment; could be useful if one
    ## wants to normalize by ot
    if (!is.null(seqinfo))
        seqinfo <- seqinfo(windows)
    obs <- mcols(windows)[[col]]
    if (by.seqnames) {
        frac <- tapply(windows, seqnames(windows),
                       function(x) {sum(obs, na.rm=TRUE) / sum(width(x))})
        exp <- as.numeric(frac[as.integer(match(seqnames(windows), names(frac)))] * width(windows))
    } else {
        frac <- sum(obs, na.rm = TRUE) / sum(width(windows))
        exp <- frac * width(windows)
    }
    exp
}
