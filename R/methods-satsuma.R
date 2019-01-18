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

##' readSatsuma
##'
##' Read satsuma one-liner output
##'
##' @export
##' @rdname readSatsuma
##'
##'
setMethod("readSatsuma", signature(con = "list"),
          function(con) {
    NULL
})
