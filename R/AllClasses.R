setClassUnion("characterOrMissing", c("character", "missing", "logical"))
setClassUnion("factorOrcharacter", c("factor", "character"))
setClassUnion("integerOrMissing", c("integer", "missing", "logical"))
setClassUnion("XStringSetOrMissing", c("XStringSet", "BStringSet", "missing", "logical"))
setClassUnion("DNAStringSetOrMissing", c("DNAStringSet", "missing", "logical", "NULL"))

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
##' @description Representation of an alignment item. AlignmentItem
##'     subclasses and extends GenomicRanges by adding two additional
##'     slots, bases and sequence.
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
##' @description RepeatAlignmentItem subclasses and extends
##'     AlignmentItem by adding one additional slot for repeat_class.
##'     Note that the RepeatMasker library fasta header uses the
##'     repeat name (here seqname) and the repeat class to identify
##'     the sequence, concatenated by a #. Hence, in order to match an
##'     entry in RepeatAlignmentItem to a library, a helper function
##'     to convert between naming systems is needed.
##'
##' @export
##' @rdname RepeatAlignmentItem-class
##'
##' @import GenomicRanges
##' @import Biostrings
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

##' @importFrom methods validObject callNextMethod
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
##' @description Pairs subclass
##'
##'
##' @details The AlignmentPairs class extends the
##'     [S4Vectors:Pairs-class]
##'
##' @export
##' @rdname AlignmentPairs-class
##'
##' @import S4Vectors
##'
##' @seealso \link[ripr]{AlignmentPairs-class}
##'
setClass("AlignmentPairs",
         contains = c("Pairs"),
         validity = .valid.AlignmentPairs)


##' @importFrom methods validObject callNextMethod
setMethod("initialize", "AlignmentPairs", function(.Object, ...) {
    .Object <- callNextMethod()
    validObject(.Object)
    .Object
})


##' List of AlignmentPairs instances
##'
##' @description Subclass of S4Vectors SimpleList, where each entry is
##'     an AlignmentPairs object
##'
##' @export
##' @rdname AlignmentPairsList-class
##'
##' @import S4Vectors
##'
setClass("AlignmentPairsList",
         contains = "SimpleList",
         prototype = prototype(elementType = "AlignmentPairs")
         )

##' List of AlignmentItem instances
##'
##' @description Subclass of S4Vectors SimpleList, where each entry is
##'     an AlignmentItem object
##'
##' @export
##' @rdname AlignmentItem-class
##'
##' @import S4Vectors
##'
setClass("AlignmentItemList",
         contains = "SimpleList",
         prototype = prototype(elementType = "AlignmentItem")
         )
