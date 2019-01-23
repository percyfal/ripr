##' AlignmentPairs
##'
##' AlignmentPairs constructor
##'
##'
##' @export
##' @rdname AlignmentPairs
##'
setMethod("AlignmentPairs", signature = c("AlignmentItem", "AlignmentItem"),
          definition = function(query, subject, ...){
    if (!missing(...)) {
        elementMetadata <- DataFrame(...)
    } else {
        elementMetadata <- new("DataFrame", nrows = length(query))
    }
    elementMetadata$query <- query
    elementMetadata$subject <- subject
    new("AlignmentPairs", first = 1:length(query), second = 1:length(subject),
        elementMetadata = elementMetadata)
})

##############################
## Getters
##############################

##' Get query
##'
##' Get query of an AlignmentPairs object
##'
##' @return AlignmentItem
##'
##' @export
##' @rdname query
##'
setMethod("query", "AlignmentPairs", function(x) mcols(x)$query)

##' Get subject
##'
##' Get subject of an AlignmentPairs object
##'
##' @param x AlignmentPairs object
##' @return AlignmentItem
##'
##' @export
##' @rdname subject
##'
setMethod("subject", "AlignmentPairs", function(x) mcols(x)$subject)

##' score
##'
##' @param x AlignmentPairs object
##' @export
##' @rdname score
##'
setMethod("score", "AlignmentPairs", function(x) mcols(x)$score)

##' divergence
##'
##' @export
##' @rdname divergence
##'
setMethod("divergence", "AlignmentPairs", function(x) mcols(x)$divergence)

##' deletions
##'
##' @export
##' @rdname deletions
##'
setMethod("deletions", "AlignmentPairs", function(x) mcols(x)$deletions)

##' insertions
##'
##' @export
##' @rdname insertions
##'
setMethod("insertions", "AlignmentPairs", function(x) mcols(x)$insertions)

##' linkage_id
##'
##' @export
##' @rdname linkage_id
##'
setMethod("linkage_id", "AlignmentPairs", function(x) mcols(x)$linkage_id)

##' count
##'
##' @export
##' @rdname count
##'
##' @param width width of oligonucleotide
##' @param step window step
##' @param which which AlignmentItem to operate on (subject or query)
##' @param exclude Exclude characters from calculation
##'
setMethod("count", "AlignmentPairs",
          function(x, width = 2, step = 1, exclude = c("X", "-"), which = "query", ...) {
    which <- match.arg(which, c("query", "subject"))
    count(mcols(x)[[which]], width, step, exclude, ...)
})

##' subseqByRef
##'
##' Retrieve subsequences from reference
##'
##' @param x AlignmentPairs
##' @param ref DNAStringSet
##' @param which which AlignmentItem to operate on
##' @param ...
##' @return XStringSet
##' @author Per Unneberg
##'
##' @rdname subseqByRef
##' @export
##'
setMethod("subseqByRef", c("AlignmentPairs", "DNAStringSet"),
          function(x, ref, which = "query", ...) {
    which <- match.arg(which, c("query", "subject"))
    subseqByRef(mcols(x)[[which]], ref, ...)
})
