## Helper: build a minimal enrichment result data frame
make_enrich_df <- function(ids, geneIDs, pvalues = NULL, nes = NULL) {
  n <- length(ids)
  data.frame(
    ID          = ids,
    Description = ids,
    pvalue      = if (is.null(pvalues)) seq(0.001, by = 0.001, length.out = n) else pvalues,
    NES         = if (is.null(nes))     rep(1, n) else nes,
    geneID      = geneIDs,
    stringsAsFactors = FALSE
  )
}

test_that("EnrichedFilter returns input unchanged when fewer than 3 rows", {
  df1 <- make_enrich_df("A", "G1/G2/G3")
  expect_identical(EnrichedFilter(df1), df1)

  df2 <- make_enrich_df(c("A", "B"), c("G1/G2", "G3/G4"))
  expect_identical(EnrichedFilter(df2), df2)
})

test_that("EnrichedFilter keeps all rows when gene sets are non-overlapping", {
  df <- make_enrich_df(
    ids     = c("A", "B", "C"),
    geneIDs = c("G1/G2/G3", "G4/G5/G6", "G7/G8/G9")
  )
  out <- EnrichedFilter(df, cutoff = 0.8)
  expect_equal(nrow(out), 3)
  expect_setequal(out$ID, c("A", "B", "C"))
})

test_that("EnrichedFilter removes the lower-ranked of two identical gene sets", {
  # A and B share all genes (Jaccard = 1.0 > 0.8); A has lower p-value so B removed
  df <- make_enrich_df(
    ids     = c("A", "B", "C"),
    geneIDs = c("G1/G2/G3", "G1/G2/G3", "G4/G5/G6"),
    pvalues = c(0.001, 0.01, 0.05)
  )
  out <- EnrichedFilter(df, cutoff = 0.8)
  expect_true("A" %in% out$ID)
  expect_false("B" %in% out$ID)
  expect_true("C" %in% out$ID)
})

test_that("EnrichedFilter keeps both pathways when Jaccard is below cutoff", {
  # A and B share 1 of 4 genes: Jaccard = 1/(4+4-1) ≈ 0.14 — well below 0.8
  df <- make_enrich_df(
    ids     = c("A", "B", "C"),
    geneIDs = c("G1/G2/G3/G4", "G1/G5/G6/G7", "G8/G9/G10/G11")
  )
  out <- EnrichedFilter(df, cutoff = 0.8)
  expect_equal(nrow(out), 3)
})

test_that("EnrichedFilter sorts by pvalue before filtering", {
  # B has lower p-value than A; both share all genes — B should survive, A removed
  df <- make_enrich_df(
    ids     = c("A", "B", "C"),
    geneIDs = c("G1/G2/G3", "G1/G2/G3", "G4/G5/G6"),
    pvalues = c(0.05, 0.001, 0.01)
  )
  out <- EnrichedFilter(df, cutoff = 0.8)
  expect_true("B" %in% out$ID)
  expect_false("A" %in% out$ID)
})

test_that("EnrichedFilter returns a data frame", {
  df <- make_enrich_df(
    ids     = c("A", "B", "C"),
    geneIDs = c("G1/G2/G3", "G1/G2/G3", "G4/G5/G6")
  )
  out <- EnrichedFilter(df)
  expect_s3_class(out, "data.frame")
})

test_that("EnrichedFilter respects the cutoff parameter", {
  # Jaccard(A,B) = 2/(3+3-2) = 0.5
  df <- make_enrich_df(
    ids     = c("A", "B", "C"),
    geneIDs = c("G1/G2/G3", "G2/G3/G4", "G5/G6/G7")
  )
  # cutoff 0.8: Jaccard 0.5 < 0.8, both kept
  out_loose <- EnrichedFilter(df, cutoff = 0.8)
  expect_equal(nrow(out_loose), 3)

  # cutoff 0.4: Jaccard 0.5 > 0.4, lower-ranked removed
  out_strict <- EnrichedFilter(df, cutoff = 0.4)
  expect_equal(nrow(out_strict), 2)
})
