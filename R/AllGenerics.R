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


##' RIPProductIndex
##'
##' Calculate the product RIP index, defined as TpA/ApT
##'
##' @export
##' @rdname RIPProductIndex
##'
##' @param x Calculate RIPProductIndex on x
##' @param ... Arguments to pass to function
##'
setGeneric("RIPProductIndex", signature = "x",
           function(x, ...)
    standardGeneric("RIPProductIndex"))

##' RIPSubstrateIndex
##'
##' Calculate the substrate RIP index, defined as (CpA + TpG) / (ApC + GpT)
##'
##' @export
##' @rdname RIPSubstrateIndex
##'
##' @param x Calculate RIPCompositeIndex on x
##' @param ... Arguments to pass to function
##'
setGeneric("RIPSubstrateIndex", signature = "x",
           function(x, ...)
    standardGeneric("RIPSubstrateIndex"))


##' RIPCompositeIndex
##'
##' Calculate the composite RIP index, defined as the
##' product index minus the substrate index
##'
##' @export
##' @rdname RIPCompositeIndex
##'
##' @param x Calculate RIPCompositeIndex on x
##' @param ... Arguments to pass to function
##'
setGeneric("RIPCompositeIndex", signature = "x",
           function(x, ...)
    standardGeneric("RIPCompositeIndex"))
