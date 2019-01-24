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

##' subseqByRef
##'
##' Retrieve subsequences from reference
##'
##' @param x AlignmentItem
##' @param ref DNAStringSet
##' @param ...
##' @return
##' @author Per Unneberg
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
##' @param ...
##'
##' @return data.frame
##' @author Per Unneberg
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
##' @param x
##' @param ref
##' @param ...
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
