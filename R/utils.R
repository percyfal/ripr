## trim whitespace
trim <- function (x) gsub("^\\s+|\\s+$", "", x)

## Defaults for NULL values
`%||%` <- function(a, b) if (is.null(a)) b else a

## Remove NULLs from a list
compact <- function(x) {
  x[!vapply(x, is.null, logical(1))]
}

## Escape regex special characters
escape_special_chars <- function(x, special_chars = "([\\\\.$^?*|+)(}{])") {
    gsub("\\]", "\\\\]", gsub("\\[", "\\\\[", gsub(special_chars, "\\\\\\1", x)))
}
