##' methods-repeatmasker
##'
##' Description Methods for dealing with RepeatMasker output
##'
##' Details
##'
##' @name RepeatMasker methods
##'
NULL

.header <- c("score", "divergence", "deletions",
             "insertions", "query_name", "query_start",
             "query_end", "query_bases", "complement",
             "matching_repeat", "repeat_class", "subject_start",
             "subject_end", "subject_bases", "linkage_id")
.aln_header <- c("score", "divergence", "deletions",
                 "insertions", "query_name", "query_start",
                 "query_end", "query_bases", "complement",
                 "repeat_label", "subject_start",
                 "subject_end", "subject_bases", "linkage_id")

.numeric_columns <- c("score", "divergence", "deletions", "insertions", "query_start",
                      "query_end", "subject_start", "subject_end")

## Scan a annotation line or alignment header line. Note that they are
## slightly different in that the summary file contains an extra
## column for repeat class/family, which often is excluded in
## alignments
.scan_line <- function(line, summary = TRUE) {
    d <- as.list(unlist(strsplit(line, "\\s+")))
    if (summary) {
        ## Remove trailing asterix if present
        if (length(d) == 16)
            d[16] <- NULL
        ## Linkage ID is missing
        if (length(d) == 14) d[[15]] <- NA
        ## Now assume 15 columns
        names(d) <- .header
        d$repeat_label <- paste(d$matching_repeat, d$repeat_class, sep="#")
    } else {
        ## matching_repeat is always(?) missing; linkage_id is always
        ## present
        if (length(d) == 13) {
            ## 13: complement missing
            names(d) <- .aln_header[-9]
            d$complement = "+"
        } else if (length(d) == 14) {
            ## 14: complement included
            names(d) <- .aln_header
        } else if (length(d) == 15) {
            ## Assume last column is irrelevant
            d <- d[1:14]
            names(d) <- .aln_header
        }
        d$matching_repeat <- unlist(strsplit(d$repeat_label, "#"))[1]
        d$repeat_class <- unlist(strsplit(d$repeat_label, "#"))[2]
    }
    if (d$complement == "C") {
        tmp <- d$subject_bases
        d$subject_bases <- d$subject_start
        d$subject_start <- tmp
    }
    unlist(d)
}

## Process header
.process_header <- function(data) {
    data[.numeric_columns] <- sapply(data[.numeric_columns], function(x) {as.numeric(as.character(x))})

    data$query_bases <- as.integer(gsub("(\\(|\\))", "", as.character(data$query_bases)))
    data$subject_bases <- as.integer(gsub("(\\(|\\))", "", as.character(data$subject_bases)))

    data$query_length <- data$query_end + data$query_bases
    data$subject_length <- data$subject_end + data$subject_bases

    data
}

.create_query_subject <- function(data, qinfo=NULL) {
    if (is.null(qinfo)) {
        message("manually inferring seqinfo for query")
        qinfo <- unique(data.frame(names = as.character(data$query_name),
                                   lengths = data$query_length))
        qinfo$names <- as.character(qinfo$names)
        if (any(duplicated(qinfo$names)))
            qinfo <- qinfo[!duplicated(qinfo$names),]
        qinfo <- GenomeInfoDb::Seqinfo(qinfo$names, qinfo$lengths)
    }
    stopifnot(inherits(qinfo, "Seqinfo"))

    strand <- rep("+", length(data$complement))
    strand[data$complement == "C"] <- "-"
    strand <- Rle(strand)

    query_seq <- NULL
    subject_seq <- NULL
    if ("query_seq" %in% colnames(data))
        query_seq <- BStringSet(data$query_seq)
    if ("subject_seq" %in% colnames(data))
        subject_seq <- BStringSet(data$subject_seq)
    i.na <- is.na(data$query_start) | is.na(data$query_end) | is.na(data$subject_start) | is.na(data$subject_end)
    if (any(i.na))
        sprintf("%s ", which(i.na))
    query <- AlignmentItem(
        seqnames = S4Vectors::Rle(data$query_name, rep(1, length(data$query_name))),
        ranges = IRanges::IRanges(start = data$query_start,
                                  end = data$query_end),
        strand = strand,
        bases = data$query_bases,
        sequence = query_seq,
        seqinfo = qinfo)

    subject <- RepeatAlignmentItem(
        seqnames = S4Vectors::Rle(data$matching_repeat, rep(1, length(data$matching_repeat))),
        ranges = IRanges::IRanges(start = data$subject_start,
                                  end = data$subject_end),
        strand = strand,
        bases = data$subject_bases,
        sequence = subject_seq,
        seqinfo = NULL,
        repeat_class = data$repeat_class
    )
    list(query = query, subject = subject)
}
##' Read RepeatMasker annotation output
##'
##' @description Parse RepeatMasker annotation output. The
##'     RepeatMasker annotation output is a simple tabular format
##'     where each line consists of a query-subject pair and
##'     additional statistics such as score and percentage divergence.
##'     \code{readRepeatMaskerAnnotation} converts each to an entry in
##'     an \code{\link[ripr]{AlignmentPairs}} object.
##'     \code{AlignmentPairs} has slot names \code{score},
##'     \code{divergence}, \code{deletions} and \code{insertions}
##'     corresponding to the first four columns in the tabular output,
##'     and a slot \code{linkage_id} that maps to the last column.
##'     There are also two slots \code{query} and \code{subject} that
##'     correspond to the query and subject columns. These slots are
##'     populated by \code{GenomicRanges::GRanges}-derived objects
##'     \code{\link[ripr]{AlignmentItem}} and
##'     \code{\link[ripr]{RepeatAlignmentItem}}.
##'
##'     A \code{GenomeInfoDb::Seqinfo} can be passed as an option to
##'     setup sequence information for the query (genome) sequence.
##'     Otherwise, sequence information is inferred from the alignment
##'     coordinates.
##'
##'     No attempt is made to infer the sequence information for the
##'     subject (repeat) for two reasons:
##'
##'     1. low-complexity repeats (e.g. polyA-sequences) have no
##'        well-defined sequence length
##'
##'     2. RepeatMasker makes adjustments to alignments which are
##'        recorded in the annotation output. Therefore the
##'        coordinates are sometimes modified which can result in
##'        inferred alignment ranges being longer than the repeat
##'        sequence itself
##'
##'
##'
##' @param filename File to parse
##' @param seqinfo Seqinfo object for query (genome) sequence
##' @param ... Arguments to pass to data access functions
##'
##' @export
##' @rdname readRepeatMaskerAnnotation
##'
setGeneric("readRepeatMaskerAnnotation", function(filename, seqinfo=NULL, ...)
    standardGeneric("readRepeatMaskerAnnotation"))

##'
##' @rdname readRepeatMaskerAnnotation
##' @export
##'
##' @examples
##' fn <- system.file("extdata", "repeatmasker_summary.txt",
##'                   package="ripr")
##' ap <- readRepeatMaskerAnnotation(fn)
##'
setMethod("readRepeatMaskerAnnotation", signature = c("character"),
          definition = function (filename, seqinfo=NULL, ...) {
    start_time <- Sys.time()
    con <- file(filename, "r")
    on.exit(close(con))
    message("Reading file ", filename)
    ## Discard header
    readLines(con, n = 3)
    ## Read rest of file
    lines <- lapply(readLines(con), trim)

    data <- as.data.frame(do.call(rbind, lapply(lines, .scan_line)))
    data <- .process_header(data)
    obj <- .create_query_subject(data, qinfo=seqinfo)

    message("Processed ", length(lines), " lines in ",
            format(Sys.time() - start_time, digits = 2))
    AlignmentPairs(obj$query, obj$subject, score = data$score, divergence = data$divergence,
                         deletions = data$deletions, insertions = data$insertions,
                         linkage_id = as.character(data$linkage_id))
})



## Get next alignment chunk
.nextAlignmentChunk <- function(con, buf, chunksize = 10, regex = "^\\d+",
                                comment.char = "^#") {
    while (TRUE) {
        ## Read more lines if there is no alignment in buffer
        if (sum(grepl(regex, buf)) <= 1)
            lines <- readLines(con, n = chunksize)
        buf <- c(buf, lines)
        pos <- grep(regex, buf)
        eof <- FALSE
        start <- pos[1]
        ## If no more lines were read we are at end
        if (length(lines) == 0) {
            end <- grep(comment.char, buf)[1] - 1
            if (is.na(end))
                end <- length(buf)
            eof <- TRUE
        } else {
            end <- pos[2] - 1
        }
        if (!(is.na(start) | is.na(end))) {
            aln <- buf[start:end]
            buf[start:end] <- NULL
            break
        }
    }
    res <- list(aln, buf, eof)
    names(res) <- c("aln", "buf", "eof")
    res
}


.repeatMaskerAlignment <- function(aln) {
    data <- as.data.frame(t(.scan_line(aln[[1]], summary = FALSE)))
    data[.numeric_columns] <- sapply(data[.numeric_columns], function(x) {as.numeric(as.character(x))})
    data <- .process_header(data)

    ## Simpler: just match alignment lines and assume odd are query,
    ## even are subject
    aln_regex <- "^C?\\s+([^ ]+)\\s+(\\d+)\\s+([^ ]+)\\s+(\\d+).*$"
    seqs <- gsub(aln_regex, "\\3", aln[grepl(aln_regex, aln)])
    if (length(seqs) == 0) {
        warning("No sequences in alignment chunk! Skipping!")
        message(aln)
        return (NULL)
    }
    data$query_seq <- paste0(seqs[seq(1, length(seqs), 2)], collapse = "")
    data$subject_seq <- paste0(seqs[seq(2, length(seqs), 2)], collapse = "")

    if (data$complement == "C")
        data$subject_seq <- paste0(
            rev(unlist(strsplit(data$subject_seq, ""))),
            collapse = "")

    data
}


##' Read repeatmasker alignment output
##'
##' @param filename File to parse
##' @param ... Arguments to pass to data access functions
##'
##' @export
##' @rdname readRepeatMaskerAlignment
##' @importFrom GenomeInfoDb Seqinfo
##'
setGeneric("readRepeatMaskerAlignment", function(filename, ...)
    standardGeneric("readRepeatMaskerAlignment"))

##'
##' @rdname readRepeatMaskerAlignment
##' @export
##'
##' @param num_aln number of alignments to read
##' @param update_freq print progress message every update_freq alignment
##'
setMethod("readRepeatMaskerAlignment", signature = "character",
          definition = function (filename, num_aln = NULL,
                                 update_freq = NULL, ...) {
    start_time <- Sys.time()
    old_time <- start_time
    con <- file(filename, "r")
    on.exit(close(con))
    message("Reading repeatmasker alignment file ", filename)
    buf <- list()
    nlines <- 0
    localenv <- new.env()
    localenv$counter <- 0
    localenv$size <- 1
    data <- data.frame()
    localenv$res <- list(NULL)
    ## https://stackoverflow.com/questions/17046336/here-we-go-again-append-an-element-to-a-list-in-r
    AddItemDoubling <- function(item) {
        if (localenv$counter == localenv$size )
        {
            length (localenv$res) <- localenv$size <- localenv$size * 2
        }
        localenv$counter <- localenv$counter + 1
        localenv$res[[localenv$counter]] <- item
    }

    while (TRUE) {
        chunk <- .nextAlignmentChunk(con, buf)
        nlines <- nlines + length(chunk$aln)
        ## Process alignment
        if (!is.null(chunk$aln)) {
            AddItemDoubling(chunk$aln)
            if (!is.null(update_freq)) {
                if (localenv$counter %% update_freq  == 0) {
                    now_time <- Sys.time()
                    message("Read ", localenv$counter, " alignments in ",
                            format(now_time - old_time, digits = 2))
                    old_time <- now_time
                }
            }
            if (!is.null(num_aln))
                if (localenv$counter %% num_aln == 0)
                    break
        }
        if (chunk$eof) break
        buf <- chunk$buf
    }
    localenv$res <- compact(localenv$res)
    message("Converting alignments to data frame objects...")
    data <- do.call("rbind", compact(lapply(localenv$res, .repeatMaskerAlignment)))
    obj <- .create_query_subject(data)
    message("Processed ", nlines, " lines in ",
            format(Sys.time() - start_time, digits = 2))
    AlignmentPairs(
        obj$query, obj$subject, score = data$score,
        divergence = data$divergence, deletions = data$deletions,
        insertions = data$insertions,
        linkage_id = as.character(data$linkage_id)
    )
})

##' asRepeatMaskerAnnotation
##'
##' Convert AlignmentPairs or AlignmentPairsList into RepeatMasker Annotation format
##'
##'
##' @param obj AlignmentPairs or AlignmentPairsList object
##'
##' @return data frame with columns corresponding to RepeatMasker Annotation format
##'
##' @export
##' @rdname asRepeatMaskerAnnotation
##'
asRepeatMaskerAnnotation <- function(obj) {
    cols <- c("score", "divergence", "deletions", "insertions",
              "query.seqnames", "query.start", "query.end",
              "query.bases", "query.strand", "subject.seqnames",
              "subject.repeat_class", "subject.start", "subject.end",
              "subject.bases", "linkage_id")
    x <- as.data.frame(obj)
    y <- x[, cols]

    ## Add parentheses to query.bases and subject.bases
    y$query.bases <- paste("(", y$query.bases, ")", sep = "")
    y$subject.bases <- paste("(", y$subject.bases, ")", sep = "")

    ## shift subject columns for inverse hits
    i <- y$query.strand == "-"
    tmp <- y[i, ]$subject.start
    y[i, ]$subject.start <- y[i, ]$subject.bases
    y[i, ]$subject.bases <- tmp

    y$query.strand <- as.character(y$query.strand)
    y[i, ]$query.strand <- "C"
    y
}
