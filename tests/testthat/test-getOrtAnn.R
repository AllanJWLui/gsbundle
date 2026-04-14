## Tests for getOrtAnn() cache behaviour

test_that("getOrtAnn returns cached RDS without making network calls", {
  tmp <- tempfile()
  on.exit(unlink(tmp, recursive = TRUE))
  dir.create(tmp)

  ## Build a minimal ortholog data frame matching the expected column structure
  fake_ort <- data.frame(
    mmu_symbol = c("Trp53", "Brca1", "Myc"),
    mmu_entrez = c("22059", "12189", "17869"),
    hsa_symbol = c("TP53",  "BRCA1", "MYC"),
    hsa_entrez = c("7157",  "672",   "4609"),
    stringsAsFactors = FALSE
  )
  rds_path <- file.path(tmp, "GeneID_Annotation_mmu_hsa.rds")
  saveRDS(fake_ort, rds_path)

  withr::with_options(list(gsbundle.cache = tmp), {
    result <- getOrtAnn(fromOrg = "mmu", toOrg = "hsa")
  })

  expect_s3_class(result, "data.frame")
  expect_named(result, c("mmu_symbol", "mmu_entrez", "hsa_symbol", "hsa_entrez"))
  expect_equal(nrow(result), 3L)
  expect_equal(result$hsa_symbol, c("TP53", "BRCA1", "MYC"))
})
