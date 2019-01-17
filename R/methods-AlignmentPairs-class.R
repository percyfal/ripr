##' AlignmentPairs
##'
##' AlignmentPairs constructor
##'
##' @title AlignmentPairs
##' @export
##' @rdname AlignmentPairs
##'
setMethod("AlignmentPairs", signature=c("AlignmentItem", "AlignmentItem"),
          definition=function(query, subject, ...){
    if (!missing(...)) {
        elementMetadata <- DataFrame(query, subject)
    } else {
        elementMetadata <- new("DataFrame", nrows=length(query))
    }
    elementMetadata$query <- query
    elementMetadata$subject <- subject
    new("AlignmentPairs", first=1:length(query), second=1:length(subject),
        elementMetadata = elementMetadata)
})



setGeneric("query", function(x, ...) standardGeneric("query"))
##' Get query
##'
##' Get query of an AlignmentPairs object
##' @title query
##' @param x AlignmentPairs object to inspect
##' @param ... additional parameters
##' @return AlignmentItem
##'
##' @export
##' @rdname query
##'
setMethod("query", "AlignmentPairs", function(x) mcols(x)$query)

##' Get subject
##'
##' Get subject of an AlignmentPairs object
##' @title subject
##' @param x AlignmentPairs object to inspect
##' @param ... additional parameters
##' @return AlignmentItem
##'
##' @export
##' @rdname subject
##'
setMethod("subject", "AlignmentPairs", function(x) mcols(x)$subject)
