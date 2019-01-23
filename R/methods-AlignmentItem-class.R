##' AlignmentItem
##'
##' AlignmentItem constructor
##'
##'
##' @export
##' @rdname AlignmentItem
##'
##'
##' @param seqnames sequence names
##' @param ranges IRanges
##' @param strand sequence strand
##' @param ... Additional parameters to pass to GRanges constructor
##' @param seqlengths sequence lengths
##' @param seqinfo Seqinfo object
##' @param bases bases from alignment end point to end of sequence
##' @param sequences XStringSet sequences
##'
AlignmentItem <- function(seqnames=NULL, ranges=NULL, strand=NULL,
                          ..., seqlengths=NULL, seqinfo=NULL, bases=NULL,
                          sequence = NULL) {
    gr <- GRanges(seqnames = seqnames, ranges = ranges, strand = strand,
                  ..., seqlengths = seqlengths, seqinfo = seqinfo)
    if (is.null(bases))
        bases <- rep(NA, length(gr))
    if (is.null(sequence))
        sequence <- BStringSet(rep("", length(gr)))
    ai <- new("AlignmentItem", gr, bases = bases, sequence = sequence)
    ai
}

##' count
##'
##' @export
##' @rdname count
##'
##' @param width width of oligonucleotide
##' @param step window step
##' @param exclude Exclude characters from calculation
##'
setMethod("count", "AlignmentItem",
          function(x, width = 2, step = 1, exclude = c("X", "-"), ...) {
    count(x@sequence, width, step, exclude, ...)
})

##'
##' @export
##' @rdname RIPProductIndex
##'
##'
setMethod("RIPProductIndex", "AlignmentItem",
          function(x, ...) {
    counts <- count(x)
    TpA <- counts[, "TA"]
    ApT <- counts[, "AT"]
    TpA / ApT
})

##'
##' @export
##' @rdname RIPSubstrateIndex
##'
##'
setMethod("RIPSubstrateIndex", "AlignmentItem",
          function(x, ...) {
    counts <- count(x)
    CpA <- counts[, "CA"]
    TpG <- counts[, "TG"]
    ApC <- counts[, "AC"]
    GpT <- counts[, "GT"]
    (CpA + TpG) / (ApC + GpT)
})

##'
##' @export
##' @rdname RIPCompositeIndex
##'
##'
setMethod("RIPCompositeIndex", "AlignmentItem",
          function(x, ...) {
    RIPProductIndex(x) - RIPSubstrateIndex(x)
})
