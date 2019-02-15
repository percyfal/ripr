##' count
##'
##' @export
##' @rdname count
##'
##' @param width width of oligonucleotide
##' @param step window step
##' @param exclude Exclude characters from calculation
##'
setMethod("count", "XStringSet",
          function(x, width = 2, step = 1, exclude = c("X", "-"), ...) {
    seq <- gsub(paste(exclude, collapse = "|"), "", x)
    oligonucleotideFrequency(DNAStringSet(seq), width = width, step = step, ...)
})

##' slidingWindows
##'
##'
##'
##' @param x
##' @param width
##' @param step
##'
##' @importFrom IRanges slidingWindows
##' @importFrom GenomicRanges seqinfo
##'
##' @return compressedGRangesList object
##'
setMethod("slidingWindows", "DNAStringSet", function(x, width, step = 1L) {
    ## Convert DNAStringSet to IRanges based on seqinfo
    y <- GRanges(seqnames = seqinfo(x)@seqnames, ranges = IRanges(start = 1, end = seqinfo(x)@seqlengths))
    slidingWindows(y, width = width, step = step)
})
