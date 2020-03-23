##' AlignmentPairs constructor
##'
##' @param query Alignment query
##' @param subject Alignment subject
##' @param ... Arguments to pass to constructor
##' @return AlignmentPair object
##'
##' @export
##' @rdname AlignmentPairs
##'
setGeneric("AlignmentPairs", signature=c("query", "subject"),
           function(query, subject, ...)
    standardGeneric("AlignmentPairs"))


##' query generic
##'
##' Generic function to retrieve query
##'
##' @param x Object from which to retrieve query
##' @param ... Arguments to pass to function
##'
##' @export
##' @rdname query
##'
setGeneric("query", function(x, ...) standardGeneric("query"))

##' divergence
##'
##' Get divergence from an alignment
##'
##' @export
##' @rdname divergence
##'
##' @param x Object from which to retrieve divergence
##' @param ... Arguments to pass to function
##'
setGeneric("divergence", function(x, ...) standardGeneric("divergence"))

##' deletions
##'
##' @export
##' @rdname deletions
##'
##' @param x Object from which to retrieve deletions
##' @param ... Arguments to pass to function
##'
##'
setGeneric("deletions", function(x, ...) standardGeneric("deletions"))

##' insertions
##'
##' @export
##' @rdname insertions
##'
##' @param x Object from which to retrieve insertions
##' @param ... Arguments to pass to function
##'
##'
setGeneric("insertions", function(x, ...) standardGeneric("insertions"))

##' linkage_id
##'
##' Get repeatmasker linkage_id from alignment
##'
##' @export
##' @rdname linkage_id
##'
##' @param x Object from which to retrieve linkage_id
##' @param ... Arguments to pass to function
##'
##'
setGeneric("linkage_id", function(x, ...) standardGeneric("linkage_id"))

##' count
##'
##' Count nucleotides in sequence
##'
##' @export
##' @rdname count
##'
##' @param x Count nucleotides in object x
##' @param ... Arguments to pass to function
##'
setGeneric("count", signature = "x",
           function(x, width = 2, step = 1, exclude = c("X", "-"), ...)
    standardGeneric("count"))


##' count
##'
##' Count nucleotides in sequence
##'
##' @export
##' @rdname subseqByRef
##'
setGeneric("subseqByRef", signature = c("x", "ref"),
           function(x, ref, ...)
    standardGeneric("subseqByRef"))



##' AlignmentPairsList
##'
##' List of AlignmentPairs items
##'
##' @param obj object to convert
##' ##' @param ...
##' @return
##' @author Per Unneberg
##' @export
##' @rdname AlignmentPairsList
##'
setGeneric("AlignmentPairsList", signature = c("obj"),
           function(obj, ...)
    standardGeneric("AlignmentPairsList"))





##' calculateRIP
##'
##' calculate RIP scores
##'
##' @param x
##' @param sequences
##' @param metadata
##' @param ...
##'
##' @export
##' @rdname calculateRIP
##'
setGeneric("calculateRIP", signature = c("x", "ref"),
           function(x, ref, sequence = FALSE,
                    metadata = FALSE, ...)
    standardGeneric("calculateRIP"))


##' windowScore
##'
##' @param x
##' @param ref
##' @param width
##' @param step
##' @param which
##'
##' @export
##' @rdname windowScore
##'
setGeneric("windowScore", signature = c("x", "ref"),
           function(x, ref, window.size = 10000L, window.step = NULL, which = "rip", ...)
    standardGeneric("windowScore"))
