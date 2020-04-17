##' @rdname subseqByRef
##' @export
##'
setMethod("subseqByRef", c("GRanges", "DNAStringSet"),
          function(x, ref, ...) {
    subseq(ref[seqnames(x)], start = start(x), end = end(x), ...)
})
