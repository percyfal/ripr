##' AlignmentPairsList
##'
##' @export
##' @rdname AlignmentPairsList-class
##'
##' @importFrom methods new
##'
setMethod("AlignmentPairsList", "list",
          function(obj) new("AlignmentPairsList", listData = obj))



##' as.data.frame
##'
##' @description Convert AlignmentPairsList to data.frame.
##'
##' @param x AlignmentPairsList object
##' @param sequences include sequences column or not
##' @param ... additional parameters to lapply
##'
##' @return data.frame
##'
##' @export
##'
setMethod("as.data.frame", signature = "AlignmentPairsList",
          function(x, ...) {
    do.call("rbind", lapply(x, as.data.frame, ...))
})
