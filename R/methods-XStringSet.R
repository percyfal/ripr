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
