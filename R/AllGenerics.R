##' @section Constructor:
##'
##' \code{AlignmentPairs(query, subject, ...)}: Constructs an AlignmentPairs
##' object by aligning the vectors \code{query} and \code{subject}
##'
##' @param query AlignmentItem
##' @param subject AlignmentItem
##' @param ... Arguments to pass to constructor
##'
##' @return AlignmentPair object
##'
##' @export
##' @rdname AlignmentPairs-class
##'
setGeneric("AlignmentPairs", signature=c("query", "subject"),
           function(query, subject, ...)
    standardGeneric("AlignmentPairs"))


##' Retrieve score from object
##'
##' @description For any object with a score slot retrieve the score
##'     object
##'
##' @param x An AlignmentPairs object
##' @param ... Additional arguments
##'
##' @return An AlignmentItem object
##'
##' @export
##' @rdname score
##'
setGeneric("score", function(x, ...) standardGeneric("score"))

##' Retrieve query from object
##'
##' @description For any object with a query slot retrieve the query
##'     object
##'
##' @param x An AlignmentPairs object
##' @param ... Additional arguments
##'
##' @return An AlignmentItem object
##'
##' @export
##' @rdname query
##'
setGeneric("query", function(x, ...) standardGeneric("query"))

##' Set the query of an object
##'
##' @description For any object with a query slot set the
##'     query
##'
##' @param x An AlignmentPairs object
##' @param value An AlignmentItem object
##'
##' @return An AlignmentPairs object
##'
##' @export
##' @rdname query
##'
setGeneric("query<-", function(x, value) standardGeneric("query<-"))

##' Retrieve subject from object
##'
##' @description For any object with a subject slot retrieve the
##'     subject object
##'
##' @param x An AlignmentPairs object
##' @param ... Additional arguments
##'
##' @return An AlignmentItem object
##'
##' @export
##' @rdname sbjct
##'
setGeneric("sbjct", function(x, ...) standardGeneric("sbjct"))

##' Set the subject of an object
##'
##' @description For any object with a subject slot set the
##'     subject
##'
##' @param x An AlignmentPairs object
##' @param value A RepeatAlignmentItem
##'
##' @return An AlignmentPairs object
##'
##' @export
##' @rdname sbjct
##'
setGeneric("sbjct<-", function(x, value) standardGeneric("sbjct<-"))

##' divergence
##'
##' @description Get divergence from an alignment
##'
##' @export
##' @rdname divergence
##'
##' @param x An AlignmentPairs object
##' @param ... Additional arguments
##'
##' @return numeric
##'
setGeneric("divergence", function(x, ...) standardGeneric("divergence"))

##' deletions
##'
##' @description Get deletions from an alignment
##'
##' @export
##' @rdname deletions
##'
##' @param x An AlignmentPairs object
##' @param ... Additional arguments
##'
##' @return numeric
##'
setGeneric("deletions", function(x, ...) standardGeneric("deletions"))

##' insertions
##'
##' @description Get insertions from an alignment
##'
##' @export
##' @rdname insertions
##'
##' @param x An AlignmentPairs object
##' @param ... Additional arguments
##'
##' @return numeric
##'
setGeneric("insertions", function(x, ...) standardGeneric("insertions"))

##' linkage_id
##'
##' @description Get repeatmasker linkage_id from alignment
##'
##' @export
##' @rdname linkage_id
##'
##' @param x An AlignmentPairs object
##' @param ... Additional arguments
##'
##' @return character
##'
setGeneric("linkage_id", function(x, ...) standardGeneric("linkage_id"))


##' count
##'
##' @description Count nucleotides in sequence
##'
##' @export
##' @rdname count
##'
##' @param x Count nucleotides in object x
##' @param width width of oligonucleotide
##' @param step window step
##' @param exclude Exclude characters from calculation
##' @param ... additional arguments
##'
##' @return numeric
##'
setGeneric("count", signature = "x",
           function(x, width = 2, step = 1, exclude = c("X", "-"), ...)
    standardGeneric("count"))


##' subseqByRef
##'
##' @description Retrieve subsequence from a reference based on
##'     coordinates in x (typically an GRanges-derived object)
##'
##' @param x an AlignmentItem, AlignmentPairs, or GRanges object
##' @param ref DNAStringSet
##' @param ... additional parameters
##'
##' @return XStringSet
##'
##' @examples
##'
##' rm.alignment <- system.file("extdata", "repeatmasker_alignment.txt", package="ripr")
##' genome <- system.file("extdata", "g5129s420.fasta", package="ripr")
##' alnpair <- readRepeatMaskerAlignment(rm.alignment)
##' g <- Biostrings::readDNAStringSet(genome)
##'
##' ## extract AlignmentItem with query()
##' s1 <- subseqByRef(query(alnpair), g)
##' ## pass AlignmentPairs and choose query with 'which'
##' s2 <- subseqByRef(alnpair, g, which="query")
##' ## coerce AlignmentItem to GRanges
##' gr <- GenomicRanges::GRanges(query(alnpair))
##' s3 <- subseqByRef(gr, g)
##'
##'
##' @export
##'
##' @rdname subseqByRef
##'
setGeneric("subseqByRef", signature = c("x", "ref"),
           function(x, ref, ...)
    standardGeneric("subseqByRef"))



##' @section Constructor:
##'
##' \code{AlignmentPairsList(obj, ...)}: Constructs an
##' AlignmentPairsList object from the input
##'
##' @param obj list
##' @param ... ellipsis
##'
##' @export
##' @rdname AlignmentPairsList-class
##'
##' @return AlignmentPairsList object
##'
setGeneric("AlignmentPairsList", signature = c("obj"),
           function(obj, ...)
    standardGeneric("AlignmentPairsList"))



##' @section Constructor:
##'
##' \code{AlignmentItemList(obj, ...)}: Constructs an
##' AlignmentItemList object from the input
##'
##' @param obj list
##' @param ... ellipsis
##'
##' @export
##' @rdname AlignmentItemList-class
##'
##' @return AlignmentItemList object
##'
setGeneric("AlignmentItemList", signature = c("obj"),
           function(obj, ...)
    standardGeneric("AlignmentItemList"))


##' calculateRIP
##'
##' @description calculate multiple RIP scores
##'
##' @param x an AlignmentPairs, DNAStringSet object or NULL
##' @param ref a DNAStringSet object corresponding to an AlignmentItem
##'     in x, by default the query
##' @param ... additional parameters passed to any of the RIP
##'     functions \code{\link[ripr]{RIPProductIndex}},
##'     \code{\link[ripr]{RIPSubstrateIndex}}, and
##'     \code{\link[ripr]{RIPCompositeIndex}}
##'
##' @return AlignmentPairs with RIP score data columns
##'
##' @export
##' @rdname calculateRIP
##'
setGeneric("calculateRIP", signature = c("x", "ref"),
           function(x, ref=NULL, ...)
    standardGeneric("calculateRIP"))


##' windowScore
##'
##' @description calculate scores in windows
##'
##' @details Apply a sliding window to a sequence object and calculate
##'     different scores. Scores include RIP, repeat content, and GC
##'     content.
##'
##' @param windows ranges containing information on the windows on
##'     which to operate
##' @param x An AlignmentPairs object
##' @param ref A DNAStringSet object corresponding to the subject
##'     reference
##' @param which which score to calculate; any of rip, repeat.content,
##'     and gc
##' @param ... additional parameters
##'
##' @return GRanges - window ranges with values consisting of
##'     different scores
##'
##' @export
##' @rdname windowScore
##'
setGeneric("windowScore", signature = c("windows", "x", "ref"),
           function(windows, x, ref,
                    which = c("rip", "repeat.content", "gc"), ...)
    standardGeneric("windowScore"))

##' Calculate GC content in windows
##'
##' @description Calculate GC content in windows
##'
##' @param windows ranges on which to operate
##' @param ref A DNAStringSet object corresponding to the subject
##'     reference
##'
##' @return GRanges
##'
##' @export
##' @rdname windowGC
##'
setGeneric("windowGC", signature = c("windows", "ref"),
           function(windows, ref) standardGeneric("windowGC"))

##' Calculate RIP index in windows
##'
##' @description Calculate RIP index in windows
##'
##' @param windows ranges on which to operate
##' @param ref A DNAStringSet object corresponding to the subject
##'     reference
##' @param ... additional parameters for \code{count}
##'
##' @return GRanges
##'
##' @export
##' @rdname windowRIP
##'
setGeneric("windowRIP", signature = c("windows", "ref"),
           function(windows, ref, ...) standardGeneric("windowRIP"))

##' Calculate repeat content in windows
##'
##' @description Calculate repeat content in windows
##'
##' @param windows ranges on which to operate
##' @param obj ranges to intersect with windows
##' @param ref A DNAStringSet object corresponding to the subject
##'     reference
##'
##' @return GRanges
##'
##' @export
##' @rdname windowRepeatContent
##'
setGeneric("windowRepeatContent", signature = c("windows", "obj", "ref"),
           function(windows, obj, ref) standardGeneric("windowRepeatContent"))
