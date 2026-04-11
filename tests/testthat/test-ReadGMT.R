## Helper: write a temporary GMT file and return its path
write_gmt <- function(lines) {
  path <- tempfile(fileext = ".gmt")
  writeLines(lines, path)
  path
}

test_that("ReadGMT parses a basic GMT file correctly", {
  gmt <- write_gmt(c(
    "PathwayA\thttp://example.com\tGENE1\tGENE2\tGENE3",
    "PathwayB\thttp://example.com\tGENE4\tGENE5"
  ))
  out <- ReadGMT(gmt)
  expect_s3_class(out, "data.frame")
  expect_named(out, c("Gene", "PathwayID", "PathwayName"))
  expect_equal(nrow(out), 5)
  expect_setequal(out$Gene[out$PathwayID == "PathwayA"],
                  c("GENE1", "GENE2", "GENE3"))
  expect_setequal(out$Gene[out$PathwayID == "PathwayB"],
                  c("GENE4", "GENE5"))
  # PathwayID and PathwayName are both set to the pathway name field
  expect_equal(unique(out$PathwayName[out$PathwayID == "PathwayA"]), "PathwayA")
})

test_that("ReadGMT applies lower size limit", {
  gmt <- write_gmt(c(
    "Big\tNA\tG1\tG2\tG3\tG4\tG5",
    "Small\tNA\tG6"                   # only 1 gene — below limit=c(2,Inf)
  ))
  out <- ReadGMT(gmt, limit = c(2, Inf))
  expect_true(all(out$PathwayID == "Big"))
  expect_false("Small" %in% out$PathwayID)
})

test_that("ReadGMT applies upper size limit", {
  gmt <- write_gmt(c(
    "Huge\tNA\tG1\tG2\tG3\tG4\tG5\tG6",  # 6 genes — above limit
    "Fine\tNA\tG7\tG8\tG9"
  ))
  out <- ReadGMT(gmt, limit = c(1, 4))
  expect_false("Huge" %in% out$PathwayID)
  expect_true("Fine" %in% out$PathwayID)
})

test_that("ReadGMT drops lines with fewer than 3 fields", {
  gmt <- write_gmt(c(
    "PathwayA\thttp://example.com\tGENE1\tGENE2",
    "TooShort\tNA"                              # only 2 fields
  ))
  out <- ReadGMT(gmt)
  expect_false("TooShort" %in% out$PathwayID)
  expect_true("PathwayA" %in% out$PathwayID)
})

test_that("ReadGMT filters out empty gene strings", {
  gmt <- write_gmt(c(
    "PathwayA\tNA\tGENE1\t\tGENE2\t"   # two empty fields
  ))
  out <- ReadGMT(gmt)
  expect_equal(nrow(out), 2)
  expect_setequal(out$Gene, c("GENE1", "GENE2"))
})

test_that("ReadGMT returns empty data frame for empty file", {
  gmt <- write_gmt(character(0))
  out <- ReadGMT(gmt)
  expect_s3_class(out, "data.frame")
  expect_equal(nrow(out), 0)
})

test_that("ReadGMT default limit keeps all pathways", {
  gmt <- write_gmt(c(
    "P1\tNA\tG1",
    "P2\tNA\tG1\tG2\tG3\tG4\tG5\tG6\tG7\tG8\tG9\tG10"
  ))
  out <- ReadGMT(gmt)   # default limit = c(1, Inf)
  expect_setequal(unique(out$PathwayID), c("P1", "P2"))
})
