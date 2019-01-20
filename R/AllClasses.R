setClassUnion("characterOrNA", c("character", "logical"))
setClassUnion("integerOrMissing", c("integer", "missing", "logical"))
setClassUnion("XStringSetOrMissing", c("XStringSet", "BStringSet", "missing", "logical"))

.valid.AlignmentItem <- function(object)
{
    if (length(object)) {
        if (!(length(object) == length(object@bases)))
            return("'bases' slot must be of same length as ranges")
        if (!(length(object) == length(object@sequence)))
            return("'sequence' slot must be of same length as ranges")
    }
}

##' Representation of an alignment item
##'
##' Representation of an alignment item. AlignmentItem subclasses and
##' extends GenomicRanges by adding two additional slots, bases and
##' sequence.
##'
##'
##' @export
##' @rdname AlignmentItem-class
##'
##' @import GenomicRanges
##' @import Biostrings
##'
setClass("AlignmentItem",
         representation = representation(
             bases = "integerOrMissing",
             sequence = "XStringSetOrMissing"
         ),
         contains = "GRanges",
         validity = .valid.AlignmentItem)


setMethod(GenomicRanges:::extraColumnSlotNames, "AlignmentItem",
          function(x) {
    c("bases", "sequence")
})


.valid.RepeatAlignmentItem <- function(object)
{
    if (length(object)) {
        if (!(length(object) == length(object@bases)))
            return("'bases' slot must be of same length as ranges")
        if (!(length(object) == length(object@sequence)))
            return("'sequence' slot must be of same length as ranges")
        if (!(length(object) == length(object@repeat_class)))
            return("'repeat_class' slot must be of same length as ranges")
    }
}

##' Representation of a repeat library alignment item
##'
##' RepeatAlignmentItem subclasses and extends AlignmentItem by adding
##' one additional slot for repeat_class. Note that the RepeatMasker
##' library fasta header uses the repeat name (here seqname) and the
##' repeat class to identify the sequence, concatenated by a #. Hence,
##' in order to match an entry in RepeatAlignmentItem to a library, a
##' helper function to convert between naming systems is needed.
##'
##' @export
##' @rdname RepeatAlignmentItem-class
##'
setClass("RepeatAlignmentItem",
         representation = representation(
             repeat_class = "factor"
         ),
         contains = c("AlignmentItem"),
         validity = .valid.RepeatAlignmentItem)


setMethod(GenomicRanges:::extraColumnSlotNames, "RepeatAlignmentItem",
          function(x) {
    c("bases", "sequence", "repeat_class")
})

setMethod("initialize", "RepeatAlignmentItem", function(.Object, ...) {
    .Object <- callNextMethod()
    validObject(.Object)
    .Object
})



.valid.AlignmentPairs <- function(object) {
    cols <- c("score", "divergence", "deletions", "insertions", "linkage_id", "query", "subject")
    if (length(object)) {
        if (!(all(cols %in% colnames(elementMetadata(object)))))
            return(paste0("AlignmentPairs object must have elementMetadata columns ", paste(cols, collapse=",")))
    }
}


##' Representation of an alignment pair
##'
##' Pairs subclass
##'
##'
##' @export
##' @rdname AlignmentPairs-class
##'
##' @import S4Vectors
##'
setClass("AlignmentPairs",
         contains = c("Pairs"),
         validity = .valid.AlignmentPairs)


setMethod("initialize", "AlignmentPairs", function(.Object, ...) {
    .Object <- callNextMethod()
    validObject(.Object)
    .Object
})
