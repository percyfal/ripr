#' ripr: A package for calculating rip statistics
#'
#' The ripr package provides functionality for calculating rip statistics
#'
#'
#' @section ripr functions:
#'
#' plotRip
#'
#' @docType package
#' @name ripr
#'
#' @import methods
#'
NULL


### ------------------------------
### RIP functions
###

##' Calculate RIP Product Index
##'
##' Calculate RIP product index defined as TpA / ApT
##'
##' @param x an AligmnentPairs, AlignmentItem or XStringSet object
##' @param ...
##' @return RIP product index
##' @author Per Unneberg
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
##' Calculate RIP Substrate Index defined as (CpA + TpG) / (ApC + GpT)
##'
##' @param x an AligmnentPairs, AlignmentItem or XStringSet object
##' @param ...
##' @return RIP Substrate Index
##' @author Per Unneberg
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
##' Calculate RIP Composite Index defined as RIP Product Index minus RIP Substrate Index
##'
##' @param x  an AligmnentPairs, AlignmentItem or XStringSet object
##' @param ...
##' @return RIP Composite Index
##' @author Per Unneberg
##' @export
##'
RIPCompositeIndex <- function(x, ...) {
    RIPProductIndex(x, ...) - RIPSubstrateIndex(x, ...)
}
