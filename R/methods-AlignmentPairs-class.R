##' AlignmentPairs
##'
##' AlignmentPairs constructor
##'
##'
##' @export
##' @rdname AlignmentPairs
##'
setMethod("AlignmentPairs", signature=c("AlignmentItem", "AlignmentItem"),
          definition=function(query, subject, ...){
    if (!missing(...)) {
        elementMetadata <- DataFrame(...)
    } else {
        elementMetadata <- new("DataFrame", nrows=length(query))
    }
    elementMetadata$query <- query
    elementMetadata$subject <- subject
    new("AlignmentPairs", first=1:length(query), second=1:length(subject),
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
          function(x, width = 2, step = 1, which = "query", exclude = c("X", "-"), ...) {
    which <- match.arg(which, c("query", "subject"))
    seq <- gsub(paste(exclude, collapse = "|"), "", mcols(x)[[which]]$sequence)
    oligonucleotideFrequency(DNAStringSet(seq), width = width, step = step, ...)
})

##'
##' @export
##' @rdname RIPProductIndex
##'
##'
setMethod("RIPProductIndex", "AlignmentPairs",
          function(x, ...) {
    counts <- count(x)
    TpA <- counts[, "TA"]
    ApT <- counts[, "AT"]
    TpA / ApT
})

##'
##' @export
##' @rdname RIPSubstrateIndex
##'
##'
setMethod("RIPSubstrateIndex", "AlignmentPairs",
          function(x, ...) {
    counts <- count(x)
    CpA <- counts[, "CA"]
    TpG <- counts[, "TG"]
    ApC <- counts[, "AC"]
    GpT <- counts[, "GT"]
    (CpA + TpG) / (ApC + GpT)
})

##'
##' @export
##' @rdname RIPCompositeIndex
##'
##'
setMethod("RIPCompositeIndex", "AlignmentPairs",
          function(x, ...) {
    RIPProductIndex(x) - RIPSubstrateIndex(x)
})
