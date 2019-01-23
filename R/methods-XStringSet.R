##' count
##'
##' @export
##' @rdname count
##'
##' @param width width of oligonucleotide
##' @param step window step
##' @param exclude Exclude characters from calculation
##'
setMethod("count", "XStringSet",
          function(x, width = 2, step = 1, exclude = c("X", "-"), ...) {
    seq <- gsub(paste(exclude, collapse = "|"), "", x)
    oligonucleotideFrequency(DNAStringSet(seq), width = width, step = step, ...)
})

##'
##' @export
##' @rdname RIPProductIndex
##'
##'
setMethod("RIPProductIndex", "XStringSet",
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
setMethod("RIPSubstrateIndex", "XStringSet",
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
setMethod("RIPCompositeIndex", "AlignmentItem",
          function(x, ...) {
    RIPProductIndex(x) - RIPSubstrateIndex(x)
})
