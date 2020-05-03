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
    y <- GRanges(seqnames = seqinfo(x)@seqnames,
                 ranges = IRanges(start = 1, end = seqinfo(x)@seqlengths))
    IRanges::slidingWindows(y, width = width, step = step)
})



##' shuffleSeq
##'
##' @description shuffle a sequence
##'
##' @param x an XStringSet object
##' @param method method used for shuffling. Either shuffle or frequency.
##' @param collapse collapse shuffled sequence to one long pseudo-chromosome
##' @param ... additional parameters
##'
##' @return DNAStringSet with shuffled sequences
##'
##' @export
##' @rdname shuffleSeq
##'
setGeneric("shuffleSeq", signature = c("x"),
           function(x, method="shuffle", collapse=TRUE, ...)
    standardGeneric("shuffleSeq"))


##'
##' @rdname shuffleSeq
##' @export
##'
##' @importFrom Biostrings oligonucleotideFrequency
##' @importFrom S4Vectors endoapply
##'
setMethod("shuffleSeq", "DNAStringSet",
          function(x, method="shuffle", collapse=TRUE, ...) {
    method <- match.arg(method, c("shuffle", "frequency"))
    if (collapse)
        x <- DNAStringSet(unlist(x))
    if (method == "shuffle") {
        sequence <- endoapply(x, sample)
    } else if (method == "frequency") {
        nt.freq <- lapply(x,
                          function(y) {
            prop.table(as.array(oligonucleotideFrequency(y, width=1, step=1)))})
        sequence <- DNAStringSet(lapply(seq_along(x), function(i) {
            DNAString(
                paste(sample(c("A", "C", "G", "T"),
                             prob = nt.freq[[i]], length(x[[i]]),
                             replace = TRUE), collapse = ""))
        }))
    }
    sequence
})


##' Make windows from reference sequence
##'
##' @description Make windows from reference sequence
##'
##' @param obj input sequence
##' @param ... additional parameters
##'
##' @return GRanges
##'
##' @export
##' @rdname makeWindows
##'
makeWindows <- function(obj, ...) UseMethod("makeWindows", obj)

##' makeWindows.DNAStringSet
##'
##' @description Make windows from DNAStringSet
##'
##' @param window.size window size
##' @param window.step window step
##'
##' @export
##' @rdname makeWindows
##'
##' @importFrom IRanges slidingWindows
##'
##' @return GRanges
##'
makeWindows.DNAStringSet <- function(obj, window.size = 10000L, window.step = NULL, ...) {
    if (is.null(window.step)) window.step <- window.size
    windows <- unlist(slidingWindows(obj, width = window.size, step = window.step))
    windows$window <- unlist(seq_along(seqnames(windows)))
    windows$window.size <- window.size
    windows$window.step <- window.step
    seqinfo(windows) <- seqinfo(obj)
    windows
}
