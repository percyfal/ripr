##' RepeatAlignmentItem
##'
##' RepeatAlignmentItem constructor
##'
##'
##' @export
##' @rdname RepeatAlignmentItem
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
RepeatAlignmentItem <- function(seqnames=NULL, ranges=NULL, strand=NULL,
                          ..., seqlengths=NULL, seqinfo=NULL, bases=NULL,
                          sequence = NULL, repeat_class = NULL) {
    ai <- AlignmentItem(seqnames = seqnames, ranges = ranges,
                        strand = strand, ..., seqlengths = seqlengths,
                        seqinfo = seqinfo, bases = bases, sequence = sequence)
    rai <- new("RepeatAlignmentItem", ai,
               repeat_class = as.factor(repeat_class))
    rai
}
