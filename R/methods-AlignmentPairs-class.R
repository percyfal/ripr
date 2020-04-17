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

##' windowScore
##'
##'
##'
##' @param x
##' @param ref
##' @param window.size
##' @param window.step
##' @param which
##' @param metadata.which
##' @param ... windows metadata
##'
##' @return GRanges
##'
##'
setMethod("windowScore", c("AlignmentPairs", "DNAStringSet"),
          function(x, ref, window.size = 10000L, window.step = NULL,
                   which = c("rip", "repeat.content", "gc"), metadata.which = NULL,
                   lambda = list(),
                   ...) {
    dots <- list(...)
    which <- match.arg(which, c("rip", "repeat.content", "gc.content"), several.ok = TRUE)
    if (is.null(window.step)) window.step <- window.size
    windows <- unlist(slidingWindows(ref, width = window.size, step = window.step))
    windows$window <- unlist(seq_along(seqnames(windows)))
    windows$window.size <- window.size
    windows$window.step <- window.step
    seqinfo(windows) <- seqinfo(ref)
    if (!is.null(metadata)) {
        metadata.which <- match.arg(metadata.which, names(metadata(x)), several.ok = TRUE)
        for (md in metadata.which) {mcols(windows)[[md]] <- metadata(x)[[md]]}
    }
    if ("rip" %in% which) {
        message("Calculating rip scores")
        arglist <- list(x = AlignmentItem(windows), ref = ref, sequence = FALSE,
                        metadata = FALSE)
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
