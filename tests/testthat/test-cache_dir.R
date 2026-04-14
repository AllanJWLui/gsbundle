## Tests for cache directory resolution in .gsbundle_cache() and gsGetter()

test_that(".gsbundle_cache() returns the system default when nothing is set", {
  withr::with_options(list(gsbundle.cache = NULL), {
    d <- gsbundle:::.gsbundle_cache()
    expect_equal(d, tools::R_user_dir("gsbundle", which = "cache"))
  })
})

test_that(".gsbundle_cache() respects options('gsbundle.cache')", {
  tmp <- tempfile()
  on.exit(unlink(tmp, recursive = TRUE))
  withr::with_options(list(gsbundle.cache = tmp), {
    d <- gsbundle:::.gsbundle_cache()
    expect_equal(d, tmp)
    expect_true(dir.exists(tmp))
  })
})

test_that(".gsbundle_cache(path) overrides the option", {
  opt_dir   <- tempfile()
  param_dir <- tempfile()
  on.exit({ unlink(opt_dir, recursive = TRUE); unlink(param_dir, recursive = TRUE) })
  withr::with_options(list(gsbundle.cache = opt_dir), {
    d <- gsbundle:::.gsbundle_cache(param_dir)
    expect_equal(d, param_dir)
    expect_true(dir.exists(param_dir))
    expect_false(dir.exists(opt_dir))
  })
})

test_that("gsGetter reads from cache.dir when pre-built RDS files are present", {
  tmp <- tempfile()
  on.exit(unlink(tmp, recursive = TRUE))
  dir.create(tmp)

  ## Build a minimal KEGG-format data frame and save it where gsGetter expects it
  fake_kegg <- data.frame(
    EntrezID    = c("1", "2", "3"),
    PathwayID   = c("KEGG_hsa00001", "KEGG_hsa00001", "KEGG_hsa00002"),
    PathwayName = c("Pathway A", "Pathway A", "Pathway B"),
    stringsAsFactors = FALSE
  )
  saveRDS(fake_kegg, file.path(tmp, "kegg.all.entrez.hsa.rds"))

  out <- gsGetter(type = "KEGG", cache.dir = tmp)

  expect_s3_class(out, "data.frame")
  expect_named(out, c("Gene", "PathwayID", "PathwayName"))
  expect_setequal(unique(out$PathwayID), c("KEGG_hsa00001", "KEGG_hsa00002"))
})
