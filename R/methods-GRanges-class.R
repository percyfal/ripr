##' @rdname subseqByRef
##' @export
##'
setMethod("subseqByRef", c("GRanges", "DNAStringSet"),
          function(x, ref, ...) {
    subseq(ref[seqnames(x)], start = start(x), end = end(x), ...)
})


##'
##' plot
##'
##' @description Plot a GRanges object on a DNAStringSet
##'
##' @importFrom GenomicRanges GRanges
##'
setMethod("plot", c("GRanges", "DNAStringSetOrMissing"),
          function(x, y, mapping=aes(), ..., size=1, null.which=NULL) {
    stopifnot(inherits(x, "GRanges"))
    message("running plot")
    ## if (!missing(y)) {
    ##     if (!is.null(null.which))
    ##         null.which <- "shuffle"
    ##     null.which <- match.arg(null.which, c("shuffle", "frequency"), several.ok=TRUE)
    ##     ## Need function to sample positions based on alignmentpair
    ##     nullseq <- lapply(null.which, function(z) {shuffleSeq(y, method=z)})
    ##     names(nullseq) <- null.which
    ## }
    x.df <- as.data.frame(x)
    message(head(x.df))
    p <- ggplot(x.df, mapping=mapping, ...) + geom_point(size=size)
    p
})
