##' @section Constructor:
##'
##' \code{AlignmentItem(seqnames=NULL, ranges=NULL, strand=NULL, ...,
##'   seqlengths=NULL, seqinfo=NULL, bases=NULL, sequence=NULL)}:
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
##' @param sequences XStringSet sequences
##'
##'
##' @examples
##'
##' AlignmentItem(
##'   Rle("chr1"),
##'   ranges=IRanges(names=c("chr1", "chr1"), start=c(10, 20), end=c(30, 22)),
##'   bases=as.integer(c(20,18)), strand=c("+", "+")
##' )
##'
##' @return AlignmentItem
##'
##' @export
##' @rdname AlignmentItem-class
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
##'
##' @param x AlignmentItem object
##' @param sequences include sequences column or not
##' @param metadata include metadata or not
##' @param ... additional arguments to as.data.frame
##'
##' @return data.frame
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
          function(x, ref = NULL, sequence = FALSE, metadata = FALSE, ...) {
    if (missing(ref))
        cbind(as.data.frame(x, sequence = sequence, metadata = metadata),
              rip.product = RIPProductIndex(x, ...),
              rip.substrate = RIPSubstrateIndex(x, ...),
              rip.composite = RIPCompositeIndex(x, ...))
    else
        cbind(as.data.frame(x, sequence = sequence, metadata = metadata),
              rip.product = RIPProductIndex(subseqByRef(x, ref), ...),
              rip.substrate = RIPSubstrateIndex(subseqByRef(x, ref), ...),
              rip.composite = RIPCompositeIndex(subseqByRef(x, ref), ...))
})
