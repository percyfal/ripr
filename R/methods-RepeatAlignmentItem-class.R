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

##' Convert RepeatAlignmentItem to data.frame.
##'
##'
##' @param x RepeatAlignmentItem object
##' @param sequences include sequences column or not
##' @param metadata include metadata or not
##' @param ...
##'
##' @return data.frame
##'
##' @export
##'
setMethod("as.data.frame", "RepeatAlignmentItem",
          function(x, sequences = FALSE, metadata = FALSE, ...) {
    mcols_df <- as.data.frame(GRanges(x), ...)
    mcols_df[, "bases"] <- x@bases
    mcols_df[, "repeat_class"] <- x@repeat_class
    if (sequences)
        mcols_df[, "sequence"] <- x@sequence
    if (metadata) {
        md <- metadata(x)
        mcols_df[, names(md)] <- md
    }
    data.frame(mcols_df,
               stringsAsFactors = FALSE)
})
