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
##' RangedSummarizedExperiment subclass
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
##' AlignmentItem subclass
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

##' Representation of an alignment pair
##'
##' Pairs subclass
##'
##' @export
##' @rdname AlignmentPairs-class
##'
##' @import S4Vectors
##'
setClass("AlignmentPairs", contains = c("Pairs"))

setMethod("initialize", "AlignmentPairs", function(.Object, ...) {
    .Object <- callNextMethod()
    validObject(.Object)
    .Object
})
