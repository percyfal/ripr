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



##' as.data.frame
##'
##' Convert AlignmentPairs to data.frame.
##'
##'
##' @param x AlignmentPairs object
##' @param sequences include sequences column or not
##' @param metadata include metadata or not
##' @param ...
##'
##' @return data.frame
##' @author Per Unneberg
##'
##' @export
##'
setMethod("as.data.frame", "AlignmentPairs",
          function(x, sequences = FALSE, metadata = FALSE, ...) {
    mcols_df <- as.data.frame(mcols(x), ...)
    cnames <- colnames(mcols_df)
    if (!sequences)
        cnames <- cnames[!grepl("sequence", cnames)]
    if (metadata) {
        md <- metadata(x)
        mcols_df[, names(md)] <- md
        cnames <- c(cnames, names(md))
    }
    data.frame(mcols_df[, cnames],
               stringsAsFactors = FALSE)
})

##' calculateRIP
##'
##' @param x
##' @param ref
##' @param ...
##'
##' @export
##' @rdname calculateRIP
##'
setMethod("calculateRIP", c("AlignmentPairs", "DNAStringSetOrMissing"),
          function(x, ref = NULL, sequence = FALSE, metadata = FALSE, ...) {
    if (is.null(ref))
        cbind(as.data.frame(x, sequence = sequence, metadata = metadata),
              rip.product = RIPProductIndex(x, ...),
              rip.substrate = RIPSubstrateIndex(x, ...),
              rip.composite = RIPCompositeIndex(x, ...))
    else
        cbind(as.data.frame(x, sequence = sequence, metadata = metadata),
              rip.product = RIPProductIndex(subseqByRef(x, ref), ...),
              rip.substrate = RIPSubstrateIndex(subseqByRef(x, ref), ...),
              rip.composite = RIPCompositeIndex(subseqByRef(x, ref), ...))
})
