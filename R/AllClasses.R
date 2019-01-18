.characterOrNA <- setClassUnion(".characterOrNA", c("character", "logical"))
.integerOrMissing <- setClassUnion(".integerOrMissing", c("integer", "missing", "logical"))
.XStringSetOrMissing <- setClassUnion(".XStringSetOrMissing", c("XStringSet", "missing", "logical"))

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
             bases = ".integerOrMissing",
             sequence = ".XStringSetOrMissing"
         ),
         contains = "GRanges",
         validity = .valid.AlignmentItem)


setMethod(GenomicRanges:::extraColumnSlotNames, "AlignmentItem",
          function(x) {
    c("bases", "sequence")
})




##' Representation of a repeat library alignment item
##'
##' RepeatAlignmentItem subclasses and extends AlignmentItem by adding
##' two additional slots, repeat_family and repeat_class, that hold
##' information about the repeat.
##'
##' @export
##' @rdname RepeatAlignmentItem-class
##'
setClass("RepeatAlignmentItem",
         contains = c("AlignmentItem"),
         representation = representation(
             repeat_family = "factor",
             repeat_class = "factor"
         ),
         prototype = prototype(
             repeat_family = factor(),
             repeat_class = factor()
         ))

setMethod(GenomicRanges:::extraColumnSlotNames, "RepeatAlignmentItem",
          function(x) {
            c("bases", "sequence", "repeat_family", "repeat_class")
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
