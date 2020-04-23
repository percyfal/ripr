context("Test repeatmasker-related functions")

rm.summary <- system.file("extdata", "repeatmasker_summary.txt", package="ripr")
rm.alignment <- system.file("extdata", "repeatmasker_alignment.txt", package="ripr")

ap <- readRepeatMaskerAnnotation(rm.summary)
ap2 <- readRepeatMaskerAlignment(rm.alignment)

test_that("repeatmasker alignment file returns AlignmentPairs object", {
    expect_equal(length(ap), 8)
    expect_equal(sum(width(query(ap)@sequence)), 0)
    expect_equal(sum(width(sbjct(ap)@sequence)), 0)
})
