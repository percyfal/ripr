##' AlignmentPairsList
##'
##' Convert list to AlignmentPairsList
##'
##' @export
##' @rdname AligmnentPairsList
##'
##' @importFrom methods new
##'
setMethod("AlignmentPairsList", "list",
          function(obj) new("AlignmentPairsList", listData = obj))



##' as.data.frame
##'
##' Convert AlignmentPairsList to data.frame.
##'
##'
##' @param x AlignmentPairsList object
##' @param sequences include sequences column or not
##' @param ...
##' @return data.frame
##' @author Per Unneberg
##'
##' @export
##'
setMethod("as.data.frame", signature = "AlignmentPairsList",
          function(x, ...) {
    do.call("rbind", lapply(x, as.data.frame, ...))
})
