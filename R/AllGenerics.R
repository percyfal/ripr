##' Read satsuma output
##'
##'
##' @param con File connection
##' @param ... Arguments to pass to data access functions
##'
##' @export
##' @rdname readSatsuma
##'
setGeneric("readSatsuma", function(con, ...)
    standardGeneric("readSatsuma"))

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
