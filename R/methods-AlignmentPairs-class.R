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
setReplaceMethod("query",
                 signature = c("AlignmentPairs", "AlignmentItem"),
                 function(x, value) {
    mcols(x)$subject <- value
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
setReplaceMethod("sbjct",
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


##'
##' @importFrom GenomeInfoDb genome
##'
##' @export
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

##'
##' @rdname windowScore
##' @export
##'
##' @param lambda list of functions that can be applied to the tuple
##'     x, ref, and the generated windows
##'
##'
setMethod("windowScore", c("AlignmentPairs", "DNAStringSet"),
          function(x, ref, window.size = 10000L, window.step = NULL,
                   which = c("rip", "repeat.content", "gc"), ...,
                   lambda = list()) {
    dots <- list(...)
    which <- match.arg(which, c("rip", "repeat.content", "gc.content"), several.ok = TRUE)
    if (is.null(window.step)) window.step <- window.size
    windows <- unlist(slidingWindows(ref, width = window.size, step = window.step))
    windows$window <- unlist(seq_along(seqnames(windows)))
    windows$window.size <- window.size
    windows$window.step <- window.step
    seqinfo(windows) <- seqinfo(ref)
    if ("rip" %in% which) {
        message("Calculating rip scores")
        arglist <- list(x = AlignmentItem(windows), ref = ref)
        arglist <- append(arglist, dots[which(names(dots) %in% names(formals(count)))])
        .rip <- do.call(calculateRIP, arglist)
        windows$rip.product <- .rip$rip.product
        windows$rip.substrate <- .rip$rip.substrate
        windows$rip.composite <- .rip$rip.composite
    }
    if ("repeat.content" %in% which) {
        message("Calculating repeat content")
        windows$repeats.observed <- NA
        ## Find overlaps between query annotations (i.e. repeats) and windows
        hits <- findOverlaps(windows, query(x), ignore.strand = TRUE)
        obs <- c(by(as.data.frame(hits), as.data.frame(hits)$queryHits,
                    function(h){
            sum(width(intersect(ranges(query(x)[h$subjectHits]),
                                reduce(ranges(windows[h$queryHits])))))}))
        windows[as.integer(names(obs))]$repeats.observed <- obs
        ## Calculate expected repeats based on 1. genome-wide estimate; proportion repeats times width of windows
        frac.repeats <- sum(windows$repeats.observed, na.rm = TRUE) / sum(width(windows))
        windows$repeats.expected <- frac.repeats * width(windows)
        ## Calculate expected repeats per chromosome
        frac.repeats.per.chr <- tapply(windows, seqnames(windows),
                                       function(x) {sum(x$repeats.observed, na.rm = TRUE) / sum(width(x))})
        windows$repeats.expected.chr <- as.numeric(frac.repeats.per.chr[as.integer(match(seqnames(windows), names(frac.repeats.per.chr)))] * width(windows))
    }
    if ("gc.content" %in% which) {
        message("Calculating gc content")
        windows$gc <- unlist(lapply(subseqByRef(windows, ref),
                                    function(x) {sum(alphabetFrequency(x)[c("G", "C")]) / length(x) }))
    }
    for (f in names(lambda)) {
        message("Applying ", f)
        windows <- lambda[[f]](x, ref, windows)
    }
    windows
})
