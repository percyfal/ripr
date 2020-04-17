##' Calculate RIP Product Index
##'
##' @description Calculate RIP product index defined as TpA / ApT
##'
##' @param x an AligmnentPairs, AlignmentItem or XStringSet object
##' @param ... arguments to pass to count
##'
##' @return RIP product index
##'
##' @export
##'
RIPProductIndex <- function(x, ...) {
    counts <- count(x, width = 2, ...)
    TpA <- counts[, "TA"]
    ApT <- counts[, "AT"]
    TpA / ApT
}

##' Calculate RIP Substrate Index
##'
##' @description Calculate RIP Substrate Index defined as (CpA + TpG) / (ApC + GpT)
##'
##' @param x an AligmnentPairs, AlignmentItem or XStringSet object
##' @param ... arguments to pass to count
##'
##' @return RIP Substrate Index
##'
##' @export
##'
RIPSubstrateIndex <- function(x, ...) {
    counts <- count(x, width = 2, ...)
    CpA <- counts[, "CA"]
    TpG <- counts[, "TG"]
    ApC <- counts[, "AC"]
    GpT <- counts[, "GT"]
    (CpA + TpG) / (ApC + GpT)
}

##' Calculate RIP Composite Index
##'
##' @description Calculate RIP Composite Index defined as RIP Product
##'     Index minus RIP Substrate Index
##'
##' @param x  an AligmnentPairs, AlignmentItem or XStringSet object
##' @param ... additional arguments for RIPProductIndex and RIPSubstrateIndex
##'
##' @return RIP Composite Index
##'
##' @export
##'
RIPCompositeIndex <- function(x, ...) {
    RIPProductIndex(x, ...) - RIPSubstrateIndex(x, ...)
}
