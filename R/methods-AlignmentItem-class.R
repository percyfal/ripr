##' AlignmentItem
##'
##' AlignmentItem constructor
##'
##' @title AlignmentItem
##' @export
##' @rdname AlignmentItem
##'
##'
AlignmentItem <- function(seqnames=NULL, ranges=NULL, strand=NULL,
                          ..., seqlengths=NULL, seqinfo=NULL, bases=NULL,
                          sequences = NULL) {
    gr <- GRanges(seqnames = seqnames, ranges = ranges, strand = strand,
                  ..., seqlengths = seqlengths, seqinfo = seqinfo)
    if (is.null(bases))
        bases <- rep(NA, length(gr))
    if (is.null(sequences))
        sequences <- BStringSet(rep("", length(gr)))
    ai <- new("AlignmentItem", gr, bases = bases, sequence = sequences)
    ai
}
