test_that("getOrg returns correct result for all supported organisms", {
  expected <- list(
    list(input = "human",  org = "hsa", pkg = "org.Hs.eg.db"),
    list(input = "mouse",  org = "mmu", pkg = "org.Mm.eg.db"),
    list(input = "rat",    org = "rno", pkg = "org.Rn.eg.db"),
    list(input = "bovine", org = "bta", pkg = "org.Bt.eg.db"),
    list(input = "canine", org = "cfa", pkg = "org.Cf.eg.db"),
    list(input = "chimp",  org = "ptr", pkg = "org.Pt.eg.db"),
    list(input = "pig",    org = "ssc", pkg = "org.Ss.eg.db")
  )
  for (e in expected) {
    res <- getOrg(e$input)
    expect_equal(res$org, e$org,
                 info = paste("org mismatch for", e$input))
    expect_equal(res$pkg, e$pkg,
                 info = paste("pkg mismatch for", e$input))
  }
})

test_that("getOrg accepts KEGG codes as input", {
  expect_equal(getOrg("hsa")$org, "hsa")
  expect_equal(getOrg("mmu")$org, "mmu")
  expect_equal(getOrg("rno")$org, "rno")
  expect_equal(getOrg("bta")$org, "bta")
  expect_equal(getOrg("cfa")$org, "cfa")
  expect_equal(getOrg("ptr")$org, "ptr")
  expect_equal(getOrg("ssc")$org, "ssc")
})

test_that("getOrg is case-insensitive", {
  expect_equal(getOrg("Human")$org, "hsa")
  expect_equal(getOrg("HUMAN")$org, "hsa")
  expect_equal(getOrg("Mouse")$org, "mmu")
  expect_equal(getOrg("HSA")$org,   "hsa")
})

test_that("getOrg returns a list with org and pkg elements", {
  res <- getOrg("human")
  expect_type(res, "list")
  expect_true("org" %in% names(res))
  expect_true("pkg" %in% names(res))
})

test_that("getOrg errors on invalid organism", {
  expect_error(getOrg("zebra"))
  expect_error(getOrg(""))
  expect_error(getOrg("homo sapiens"))   # full Latin name not supported
})
