##' @section Constructor:
##'
##' \code{AlignmentItem(seqnames = NULL, ranges = NULL, strand = NULL,}
##' \code{..., seqlengths=NULL, seqinfo=NULL, bases=NULL, sequence=NULL)}:
##'   Creates an AligmentItem object.
##'
##' The constructor uses the same arguments as
##' \code{\link[GenomicRanges]{GRanges}} constructor and adds two new arguments:
##'
##' \code{bases} \code{NULL} or an integer vector with bases from
##'   alignment end to end of sequence
##'
##' \code{sequence} \code{NULL} or an XStringSet object containing the
##'   bases of the alignment item
##'
##'
##' @param seqnames sequence names
##' @param ranges IRanges
##' @param strand sequence strand
##' @param ... Additional parameters to pass to GRanges constructor
##' @param seqlengths sequence lengths
##' @param seqinfo Seqinfo object
##' @param bases bases from alignment end point to end of sequence
##' @param sequence XStringSet sequence representing alignment
##'
##' @examples
##'
##' AlignmentItem(
##'   S4Vectors::Rle("chr1"),
##'   ranges=IRanges::IRanges(names=c("chr1", "chr1"), start=c(10, 20), end=c(30, 22)),
##'   bases=as.integer(c(20,18)), strand=c("+", "+")
##' )
##'
##' @return AlignmentItem
##'
##' @export
##' @rdname AlignmentItem-class
##'
##' @importFrom GenomicRanges GRanges
##'
##' @seealso \code{\link[GenomicRanges]{GRanges}}
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

##'
##' @export
##' @rdname count
##'
setMethod("count", "AlignmentItem",
          function(x, width = 2, step = 1, exclude = c("X", "-"), ...) {
    count(x@sequence, width, step, exclude, ...)
})

##' subseqByRef
##'
##' @rdname subseqByRef
##' @export
##'
setMethod("subseqByRef", c("AlignmentItem", "DNAStringSet"),
          function(x, ref, ...) {
    subseq(ref[seqnames(x)], start = start(x), end = end(x), ...)
})

##' Convert AlignmentItem to data.frame.
##'
##' @param x AlignmentItem object
##' @param sequences include sequences column or not
##' @param metadata include metadata or not
##' @param ... additional arguments to as.data.frame
##'
##' @return data.frame
##'
##' @importFrom GenomicRanges GRanges
##'
##' @export
##'
setMethod("as.data.frame", "AlignmentItem",
          function(x, sequences = FALSE, metadata = FALSE, ...) {
    mcols_df <- as.data.frame(GRanges(x), ...)
    mcols_df[, "bases"] <- x@bases
    if (sequences)
        mcols_df[, "sequence"] <- x@sequence
    if (metadata) {
        md <- metadata(x)
        mcols_df[, names(md)] <- md
    }
    data.frame(mcols_df,
               stringsAsFactors = FALSE)
})

##' calculateRIP
##'
##' @export
##' @rdname calculateRIP
##'
setMethod("calculateRIP", c("AlignmentItem", "DNAStringSetOrMissing"),
          function(x, ref = NULL, ...) {
    if (missing(ref)) {
        rip.product = RIPProductIndex(x, ...)
        rip.substrate = RIPSubstrateIndex(x, ...)
        rip.composite = RIPCompositeIndex(x, ...)
    } else {
        rip.product = RIPProductIndex(subseqByRef(x, ref), ...)
        rip.substrate = RIPSubstrateIndex(subseqByRef(x, ref), ...)
        rip.composite = RIPCompositeIndex(subseqByRef(x, ref), ...)
    }
    mcols(x)$rip.product <- rip.product
    mcols(x)$rip.substrate <- rip.substrate
    mcols(x)$rip.composite <- rip.composite
    x
})



##' sample
##'
##' @description randomly sample positions on a sequence based on ranges
##'
##' @export
##' @rdname sample
##'
##' @importFrom IRanges ranges
##' @importFrom Biostrings DNAStringSet
##'
sample <- function(x, size, replace=FALSE, prob=NULL, sequence=NULL, ...) {
    if (is.null(sequence))
        return(base::sample(x, size, replace=replace, prob=prob))
    stopifnot(inherits(x, "AlignmentItem"))
    stopifnot(inherits(sequence, "DNAStringSet"))
    sequence <- unlist(sequence)
    n <- length(x)
    w <- sample(width(x), n, replace=TRUE)
    start <- sample(length(sequence), n)
    end <- start + w
    i <- which(end <= length(sequence))
    message("Dropping ", sum(end>length(sequence)), " ranges that extended beyond sequence end")
    AlignmentItem(ranges=IRanges(start=start[i], end=end[i]),
                  seqname=seqnames(x)[i],
                  sequence=DNAStringSet(lapply(i, function(j) {subseq(sequence, start=start[j], end=end[j])})))
}


##'
##' @importFrom GenomeInfoDb genome
##'
##' @export
##'
setMethod("genome", "AlignmentItem",
          function(x) {
    genome(seqinfo(x))
})
