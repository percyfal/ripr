.header <- c("score", "divergence", "deletions",
             "insertions", "query_name", "query_start",
             "query_end", "query_bases", "complement",
             "matching_repeat", "repeat_class", "subject_start",
             "subject_end", "subject_bases", "linkage_id")
.numeric_columns <- c("score", "divergence", "deletions", "insertions", "query_start",
                      "query_end", "subject_start", "subject_end")

## Scan a summary line or alignment header line
.scan_line <- function(line) {
    d <- as.list(unlist(strsplit(line, "\\s+")))
    ## Remove trailing asterix if present
    if (length(d) == 16)
        d[16] <- NULL
    ## Linkage ID is missing
    if (length(d) == 14) d[[15]] <- NA

    ## Assume 15 columns
    names(d) <- .header
    d$repeat_class <- paste(d$matching_repeat, d$repeat_class,
                            sep = "#")
    if (d$complement == "C") {
        tmp <- d$subject_bases
        d$subject_bases <- d$subject_start
        d$subject_start <- tmp
    }
    if (grepl("#", d$repeat_class)) {
        repeatinfo <- unlist(strsplit(d$repeat_class, "#"))
    } else {
        message("missing repeat name/family separator; setting both to ",
                d$repeat_class)
        repeatinfo <- c(d$repeat_class, d$repeat_class)
    }
    d$repeat_name <- repeatinfo[1]
    d$repeat_family <- repeatinfo[2]
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

.create_query_subject <- function(data) {
    qinfo <- unique(data.frame(names = as.character(data$query_name),
                               lengths = data$query_length))
    sinfo <- unique(data.frame(names = as.character(data$repeat_name),
                               lengths = data$subject_length))
    qinfo$names <- as.character(qinfo$names)
    sinfo$names <- as.character(sinfo$names)
    if (any(duplicated(qinfo$names)))
        qinfo <- qinfo[!duplicated(qinfo$names),]
    if (any(duplicated(sinfo$names)))
        sinfo <- sinfo[!duplicated(sinfo$names),]

    strand <- rep("+", length(data$complement))
    strand[data$complement == "C"] <- "-"
    strand <- Rle(strand)

    query_seq <- NULL
    subject_seq <- NULL
    if ("query_seq" %in% colnames(data))
        query_seq <- BStringSet(data$query_seq)
    if ("subject_seq" %in% colnames(data))
        subject_seq <- BStringSet(data$subject_seq)

    query <- AlignmentItem(
        seqnames = S4Vectors::Rle(data$query_name, rep(1, length(data$query_name))),
        ranges = IRanges::IRanges(start = data$query_start,
                                  end = data$query_end),
        strand = strand,
        bases = data$query_bases,
        sequence = query_seq,
        seqinfo = Seqinfo(qinfo$names, qinfo$lengths))

    subject <- AlignmentItem(
        seqnames = S4Vectors::Rle(data$repeat_name, rep(1, length(data$repeat_name))),
        ranges = IRanges::IRanges(start = data$subject_start,
                                  end = data$subject_end),
        strand = strand,
        bases = data$subject_bases,
        sequence = subject_seq,
        seqinfo = Seqinfo(sinfo$names, sinfo$lengths))
    list(query = query, subject = subject)
}
##' Read repeatmasker summary output
##'
##' @param filename File to parse
##' @param ... Arguments to pass to data access functions
##'
##' @export
##' @rdname readRepeatMaskerSummary
##'
setGeneric("readRepeatMaskerSummary", function(filename, ...)
    standardGeneric("readRepeatMaskerSummary"))

##'
##' @rdname readRepeatMaskerSummary
##' @export
##'
setMethod("readRepeatMaskerSummary", signature = "character", definition = function (filename, ...) {
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
    obj <- .create_query_subject(data)

    message("Processed ", length(lines), " lines in ", format(Sys.time() - start_time, digits=2))
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
            message("set end to ", end)
            eof <- TRUE
        } else {
            end <- pos[2] - 1
        }
        if (!is.na(end))
            message("start: ", start, " end: ", end, " pos: ", pos)
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
    data <- as.data.frame(t(.scan_line(aln[[1]])))
    data[.numeric_columns] <- sapply(data[.numeric_columns], function(x) {as.numeric(as.character(x))})

    data <- .process_header(data)

    query <- data$query_name
    subject <- data$matching_repeat

    query.regex <- paste0("^\\s+", query, "\\s+", "\\d+", "\\s+")
    i <- grepl(query.regex, aln)
    data$query_seq <- paste0(gsub("^\\s+.+\\s+\\d+\\s+(.+)\\s+\\d+\\s+", "\\1", aln[i]),
                             collapse = "")

    subject.regex <- paste0("^C?\\s+", subject, "\\s+", "\\d+", "\\s+")
    i <- grepl(subject.regex, aln)
    data$subject_seq <- paste0(gsub("^C?\\s+.+\\s+\\d+\\s+(.+)\\s+\\d+\\s+", "\\1", aln[i]),
                               collapse = "")

    if (data$complement == "C")
        data$subject_seq <- paste0(rev(unlist(strsplit(data$subject_seq, ""))), collapse="")

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
##'
setMethod("readRepeatMaskerAlignment", signature = "character", definition = function (filename, ...) {
    start_time <- Sys.time()
    con <- file(filename, "r")
    on.exit(close(con))
    message("Reading repeatmasker alignment file ", filename)
    buf <- list()
    nlines <- 0
    data <- data.frame()
    while(TRUE) {
        res <- .nextAlignmentChunk(con, buf)
        nlines <- nlines + length(res$aln)
        ## Process alignment
        if (!is.null(res$aln)) {
            data <- rbind(data, .repeatMaskerAlignment(res$aln))
        }
        if (res$eof) break
        buf <- res$buf
    }

    obj <- .create_query_subject(data)
    message("Processed ", nlines, " lines in ", format(Sys.time() - start_time, digits=2))
    AlignmentPairs(obj$query, obj$subject, score = data$score, divergence = data$divergence,
                   deletions = data$deletions, insertions = data$insertions,
                   linkage_id = as.character(data$linkage_id))

})
