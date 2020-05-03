##' @title AlignmentPairs
##'
##' @export
##' @rdname AlignmentPairs-class
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

##'
##' @export
##' @rdname query
##'
setMethod("query", "AlignmentPairs", function(x) mcols(x)$query)

##'
##' @export
##' @rdname query
##'
##' @importFrom methods validObject
##'
setMethod("query<-",
          signature = c("AlignmentPairs", "AlignmentItem"),
          function(x, value) {
    mcols(x)$query <- value
    validObject(x)
    x
})

##'
##' @export
##' @rdname sbjct
##'
setMethod("sbjct", "AlignmentPairs", function(x) mcols(x)$subject)

##'
##' @export
##' @rdname sbjct
##'
##' @importFrom methods validObject
##'
setMethod("sbjct<-",
          signature = c("AlignmentPairs", "RepeatAlignmentItem"),
          function(x, value) {
    mcols(x)$subject <- value
    validObject(x)
    x
})

##'
##' @export
##' @rdname score
##'
setMethod("score", "AlignmentPairs", function(x) mcols(x)$score)

##'
##' @export
##' @rdname divergence
##'
setMethod("divergence", "AlignmentPairs", function(x) mcols(x)$divergence)

##'
##' @export
##' @rdname deletions
##'
setMethod("deletions", "AlignmentPairs", function(x) mcols(x)$deletions)

##'
##' @export
##' @rdname insertions
##'
setMethod("insertions", "AlignmentPairs", function(x) mcols(x)$insertions)

##'
##' @export
##' @rdname linkage_id
##'
setMethod("linkage_id", "AlignmentPairs", function(x) mcols(x)$linkage_id)

##'
##' @param which which AlignmentItem to operate on (subject or query)
##'
##' @export
##' @rdname count
##'
setMethod("count", "AlignmentPairs",
          function(x, width = 2, step = 1, exclude = c("X", "-"), which = "query", ...) {
    which <- match.arg(which, c("query", "subject"))
    count(mcols(x)[[which]], width, step, exclude, ...)
})


##' genome
##'
##' @importFrom GenomeInfoDb genome
##'
##' @rdname genome
##' @export
##'
##' @param x AlignmentPairs object
##'
setMethod("genome", "AlignmentPairs",
          function(x) {
    genome(query(x))
})


##'
##' @param which which AlignmentItem to operate on
##'
##' @rdname subseqByRef
##' @export
##'
setMethod("subseqByRef", c("AlignmentPairs", "DNAStringSet"),
          function(x, ref, which = "query", ...) {
    which <- match.arg(which, c("query", "subject"))
    subseqByRef(mcols(x)[[which]], ref, ...)
})



##' as.data.frame
##'
##' Convert AlignmentPairs to data.frame.
##'
##'
##' @param x AlignmentPairs object
##' @param ... additional arguments
##'
##' @return data.frame
##'
##' @export
##'
setMethod("as.data.frame", "AlignmentPairs",
          function(x, ...) {
    mcols_df <- as.data.frame(mcols(x), ...)
    cnames <- colnames(mcols_df)
    data.frame(mcols_df[, cnames],
               stringsAsFactors = FALSE)
})

##'
##' @export
##' @rdname calculateRIP
##'
##' @examples
##'
##' rm.alignment <- system.file("extdata", "repeatmasker_alignment.txt", package="ripr")
##' genome <- Biostrings::readDNAStringSet(system.file("extdata", "g5129s420.fasta", package="ripr"))
##' ap <- readRepeatMaskerAlignment(rm.alignment)
##' cr <- calculateRIP(ap, genome)
##'
setMethod("calculateRIP",
          signature = c("AlignmentPairs", "DNAStringSetOrMissing"),
          function(x, ref = NULL, ...) {
    if (is.null(ref)) {
        rip.product = RIPProductIndex(x, ...)
        rip.substrate = RIPSubstrateIndex(x, ...)
        rip.composite = RIPCompositeIndex(x, ...)
    } else {
        rip.product = RIPProductIndex(subseqByRef(x, ref), ...)
        rip.substrate = RIPSubstrateIndex(subseqByRef(x, ref), ...)
        rip.composite = RIPCompositeIndex(subseqByRef(x, ref), ...)
    }
    mcols(x)$rip.product <- rip.product
    mcols(x)$rip.substrate <- rip.substrate
    mcols(x)$rip.composite <- rip.composite
    x
})
