test_that("format_gs_name handles _DN and _UP suffixes", {
  expect_equal(format_gs_name("KRAS_SIGNALING_DN"), "KRAS signaling down")
  expect_equal(format_gs_name("KRAS_SIGNALING_UP"), "KRAS signaling up")
  # _DN/_UP only matched at end of string
  expect_equal(format_gs_name("UPSTREAM_TARGETS"), "Upstream targets")
})

test_that("format_gs_name applies sentence case", {
  expect_equal(format_gs_name("ADIPOGENESIS"), "Adipogenesis")
  expect_equal(format_gs_name("UV_RESPONSE"), "UV response")  # UV restored by step 5
  expect_equal(format_gs_name("DNA_REPAIR"),  "DNA repair")   # DNA restored by step 5
})

test_that("format_gs_name applies multi-word substitutions (step 3)", {
  expect_equal(format_gs_name("TNFA_SIGNALING_VIA_NFKB"),
               "TNF-\u03b1 signaling via NF-\u03baB")
  expect_equal(format_gs_name("PI3K_AKT_MTOR_SIGNALING"),
               "PI3K-AKT-mTOR signaling")
  expect_equal(format_gs_name("IL6_JAK_STAT3_SIGNALING"),
               "IL-6-JAK-STAT3 signaling")
  expect_equal(format_gs_name("IL2_STAT5_SIGNALING"),
               "IL-2-STAT5 signaling")
  expect_equal(format_gs_name("INTERFERON_GAMMA_RESPONSE"),
               "IFN-\u03b3 response")
  expect_equal(format_gs_name("INTERFERON_ALPHA_RESPONSE"),
               "IFN-\u03b1 response")
  expect_equal(format_gs_name("TGF_BETA_SIGNALING"),
               "TGF-\u03b2 signaling")
  expect_equal(format_gs_name("WNT_BETA_CATENIN_SIGNALING"),
               "Wnt-\u03b2-catenin signaling")
  expect_equal(format_gs_name("EPITHELIAL_MESENCHYMAL_TRANSITION"),
               "Epithelial-mesenchymal transition")
})

test_that("format_gs_name applies single-word substitutions (step 5)", {
  expect_equal(format_gs_name("MTORC1_SIGNALING"),  "mTORC1 signaling")
  expect_equal(format_gs_name("MTOR_SIGNALING"),    "mTOR signaling")
  expect_equal(format_gs_name("P53_PATHWAY"),       "p53 pathway")
  expect_equal(format_gs_name("MYC_TARGETS_V1"),    "Myc targets v1")
  expect_equal(format_gs_name("E2F_TARGETS"),       "E2F targets")
  expect_equal(format_gs_name("G2M_CHECKPOINT"),    "G2/M checkpoint")
  expect_equal(format_gs_name("DNA_DAMAGE_RESPONSE"), "DNA damage response")
  expect_equal(format_gs_name("RNA_PROCESSING"),    "RNA processing")
})

test_that("format_gs_name handles IL with numeric suffix", {
  expect_equal(format_gs_name("IL2_SIGNALING"),  "IL-2 signaling")
  expect_equal(format_gs_name("IL6_SIGNALING"),  "IL-6 signaling")
  expect_equal(format_gs_name("IL10_SIGNALING"), "IL-10 signaling")
})

test_that("format_gs_name handles STAT with numeric suffix", {
  expect_equal(format_gs_name("STAT3_TARGETS"), "STAT3 targets")
  expect_equal(format_gs_name("STAT5_SIGNALING"), "STAT5 signaling")
})

test_that("format_gs_name is vectorised", {
  input <- c("DNA_REPAIR", "MTORC1_SIGNALING", "P53_PATHWAY")
  out   <- format_gs_name(input)
  expect_length(out, 3)
  expect_equal(out[1], "DNA repair")
  expect_equal(out[2], "mTORC1 signaling")
  expect_equal(out[3], "p53 pathway")
})

test_that("format_gs_name handles edge cases", {
  expect_equal(format_gs_name(""), "")
  expect_equal(format_gs_name(character(0)), character(0))
})
