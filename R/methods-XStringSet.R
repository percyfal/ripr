##'
##' @export
##' @rdname count
##'
setMethod("count", "XStringSet",
          function(x, width = 2, step = 1, exclude = c("X", "-"), ...) {
    seq <- gsub(paste(exclude, collapse = "|"), "", x)
    oligonucleotideFrequency(DNAStringSet(seq), width = width, step = step, ...)
})

##' slidingWindows
##'
##' @description Generate a compressedGRangesList object corresponding
##'     to a sliding window
##'
##' @param x DNAStringSet object on which to generate sliding windows
##' @param width window size
##' @param step step size
##'
##' @importFrom IRanges IRanges
##' @importFrom GenomicRanges GRanges
##' @importFrom GenomeInfoDb seqinfo
##'
##' @return compressedGRangesList object
##'
##' @seealso \code{\link[IRanges]{slidingWindows}}
##'
setMethod("slidingWindows", "DNAStringSet", function(x, width, step = 1L) {
    ## Convert DNAStringSet to GRanges based on seqinfo
    y <- GRanges(seqnames = seqinfo(x)@seqnames, ranges = IRanges(start = 1, end = seqinfo(x)@seqlengths))
    IRanges::slidingWindows(y, width = width, step = step)
})
