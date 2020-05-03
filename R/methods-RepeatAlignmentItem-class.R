##' @section Constructor:
##'
##' \code{RepeatAlignmentItem(seqnames=NULL, ranges=NULL, strand=NULL, ...,
##'   seqlengths=NULL, seqinfo=NULL, bases=NULL, sequence=NULL)}:
##'   Creates an RepeatAligmentItem object.
##'
##' The constructor uses the same arguments as
##' \code{\link[GenomicRanges]{GRanges}} constructor and adds two new arguments:
##'
##' \code{bases} \code{NULL} or an integer vector with bases from
##'   alignment end to end of sequence
##'
##' \code{sequence} \code{NULL} or an XStringSet object containing the
##'   bases of the alignment item
##'
##'
##' @param seqnames sequence names
##' @param ranges IRanges
##' @param strand sequence strand
##' @param ... Additional parameters to pass to GRanges constructor
##' @param seqlengths sequence lengths
##' @param seqinfo Seqinfo object
##' @param bases bases from alignment end point to end of sequence
##' @param sequence XStringSet corresponding to repeat sequence
##' @param repeat_class repeat sequence identifier
##'
##' @examples
##'
##' RepeatAlignmentItem(
##'   S4Vectors::Rle("chr1"),
##'   ranges=IRanges::IRanges(names=c("chr1", "chr1"), start=c(10, 20), end=c(30, 22)),
##'   bases=as.integer(c(20,18)), strand=c("+", "+"),
##'   repeat_class=c("Gypsy", "LTR")
##' )
##'
##' @return RepeatAlignmentItem
##'
##' @export
##' @rdname RepeatAlignmentItem-class
##'
##' @seealso \code{\link[GenomicRanges]{GRanges}}
##'
RepeatAlignmentItem <- function(seqnames=NULL, ranges=NULL, strand=NULL,
                          ..., seqlengths=NULL, seqinfo=NULL, bases=NULL,
                          sequence = NULL, repeat_class = NULL) {
    ai <- AlignmentItem(seqnames = seqnames, ranges = ranges,
                        strand = strand, ..., seqlengths = seqlengths,
                        seqinfo = seqinfo, bases = bases, sequence = sequence)
    rai <- new("RepeatAlignmentItem", ai,
               repeat_class = as.factor(repeat_class))
    rai
}

##' Convert RepeatAlignmentItem to data.frame.
##'
##' @param x RepeatAlignmentItem object
##' @param ... additional arguments to as.data.frame
##'
##' @return data.frame
##'
##' @export
##'
setMethod("as.data.frame", "RepeatAlignmentItem",
          function(x, ...) {
    mcols_df <- as.data.frame(GRanges(x), ...)
    mcols_df[, "bases"] <- x@bases
    mcols_df[, "repeat_class"] <- x@repeat_class
    data.frame(mcols_df,
               stringsAsFactors = FALSE)
})



##' repeatClass
##'
##' @description get the repeat_class slot
##'
##' @param x a RepeatAlignmentItem
##'
##' @return factor
##'
##' @export
##' @rdname repeatClass
##'
setGeneric("repeatClass", signature = c("x"),
           function(x)
    standardGeneric("repeatClass"))


##'
##' @rdname repeatClass
##' @export
##'
setMethod("repeatClass", "RepeatAlignmentItem",
          function(x) {
    x@repeat_class
})


##' repeatClass<-
##'
##' @description set the repeat_class slot
##'
##' @param x a RepeatAlignmentItem
##' @param value a factor or character vector
##'
##' @return factor
##'
##' @export
##' @rdname repeatClass
##'
setGeneric("repeatClass<-",
           signature = c("x", "value"),
           function(x, value)
    standardGeneric("repeatClass<-"))



##'
##' @rdname repeatClass
##' @export
##'
##' @importFrom methods validObject slot<-
##'
setMethod("repeatClass<-",
                 signature=c("RepeatAlignmentItem", "factorOrcharacter"),
                 function(x, value) {
    if (is.character(value))
        value <- factor(value)
    slot(x, "repeat_class") <- value
    validObject(x)
    x
})
