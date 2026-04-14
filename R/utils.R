## Internal helper: return (and create) the gsbundle cache directory.
## Resolution order: path argument > options("gsbundle.cache") > system default.
.gsbundle_cache <- function(path = NULL) {
  d <- if (!is.null(path)) path else
       if (!is.null(getOption("gsbundle.cache"))) getOption("gsbundle.cache") else
       tools::R_user_dir("gsbundle", which = "cache")
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  d
}

## Internal helper: save database version metadata to versions.json
.save_db_meta <- function(db, organism, meta, cache.dir = NULL) {
  if (!requireNamespace("jsonlite", quietly = TRUE))
    stop("Package 'jsonlite' is required. Please install it.", call. = FALSE)
  path <- file.path(.gsbundle_cache(cache.dir), "versions.json")
  existing <- if (file.exists(path)) jsonlite::fromJSON(path, simplifyVector = FALSE) else list()
  key <- paste0(db, "_", organism)
  existing[[key]] <- meta
  writeLines(jsonlite::toJSON(existing, auto_unbox = TRUE, pretty = TRUE), path)
  invisible(NULL)
}

## Internal helper: fetch Last-Modified header via HTTP HEAD
.get_last_modified <- function(url) {
  if (requireNamespace("httr", quietly = TRUE)) {
    tryCatch({
      h <- httr::HEAD(url)
      lm <- httr::headers(h)[["last-modified"]]
      if (!is.null(lm) && nzchar(lm)) return(lm)
    }, error = function(e) NULL)
  }
  paste0("(exact version unavailable, fetched ", format(Sys.time(), "%Y-%m-%d"), ")")
}

#' Retrieve recorded database version metadata
#'
#' Returns a data frame summarising which version of each database was
#' downloaded and cached locally by \code{\link{retrieve_gs}} and
#' \code{\link{getGeneAnn}}.
#'
#' @param cache.dir Path to the cache directory to read \code{versions.json}
#'   from. Overrides \code{options("gsbundle.cache")} and the default system
#'   cache. Should match the \code{cache.dir} used when the data were
#'   downloaded.
#'
#' @return A data frame with columns \code{database}, \code{organism},
#'   \code{version}, \code{downloaded}, and \code{source_url}. Returns an
#'   empty data frame if no metadata has been recorded yet.
#'
#' @examples
#' \dontrun{
#'   gs_versions()
#' }
#'
#' @export
gs_versions <- function(cache.dir = NULL) {
  if (!requireNamespace("jsonlite", quietly = TRUE))
    stop("Package 'jsonlite' is required. Please install it.", call. = FALSE)
  empty <- data.frame(database = character(), organism = character(),
                      version = character(), downloaded = character(),
                      source_url = character(), stringsAsFactors = FALSE)
  path <- file.path(.gsbundle_cache(cache.dir), "versions.json")
  if (!file.exists(path)) return(empty)
  raw <- jsonlite::fromJSON(path, simplifyVector = FALSE)
  if (length(raw) == 0) return(empty)
  rows <- lapply(names(raw), function(key) {
    parts <- strsplit(key, "_", fixed = TRUE)[[1]]
    db  <- parts[1]
    org <- paste(parts[-1], collapse = "_")
    m   <- raw[[key]]
    data.frame(database    = db,
               organism    = org,
               version     = as.character(m$version),
               downloaded  = as.character(m$downloaded),
               source_url  = as.character(m$source_url),
               stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

#' Format gene set names to human-readable style
#'
#' Converts SCREAMING_SNAKE_CASE MSigDB gene set names (after collection-prefix
#' stripping) to readable sentence-case labels with proper Greek letters and
#' acronym capitalisation, matching the style of the curated MSigDB Hallmark
#' name table.
#'
#' The transformation pipeline is:
#' \enumerate{
#'   \item Replace \code{_DN} / \code{_UP} suffixes with \code{_down} /
#'         \code{_up}.
#'   \item Replace underscores with spaces and convert to lower case.
#'   \item Apply multi-word substitutions for known compound terms
#'         (e.g. \code{tgf beta} \eqn{\rightarrow} \code{TGF-\beta}).
#'   \item Capitalise the first character (sentence case).
#'   \item Apply single-word substitutions to restore acronyms and terms
#'         that should not be sentence-cased (e.g. \code{dna} \eqn{\rightarrow}
#'         \code{DNA}, \code{mtorc1} \eqn{\rightarrow} \code{mTORC1},
#'         \code{p53} stays \code{p53}).
#' }
#'
#' This is applied automatically to MSigDB pathway names in
#' \code{\link{gsGetter}}. KEGG, Reactome, GO and CORUM names are already in
#' natural-language form and are not modified.
#'
#' @param x Character vector of gene set names in SCREAMING_SNAKE_CASE (with
#'   the leading collection prefix already stripped).
#'
#' @return Character vector of the same length as \code{x}.
#'
#' @examples
#' format_gs_name(c("TNFA_SIGNALING_VIA_NFKB", "PI3K_AKT_MTOR_SIGNALING",
#'                  "MTORC1_SIGNALING", "P53_PATHWAY",
#'                  "KRAS_SIGNALING_DN", "PANCREAS_BETA_CELLS"))
#'
#' @export
format_gs_name <- function(x) {
  ## Step 1: normalise _UP / _DN suffixes before case folding
  x <- gsub("_DN$", "_down", x)
  x <- gsub("_UP$", "_up",   x)

  ## Step 2: underscores -> spaces, then full lower-case
  x <- tolower(gsub("_", " ", x))

  ## Step 3: multi-word substitutions (longest / most-specific first)
  ##   Greek letters: \\u03b1 = alpha, \\u03b2 = beta, \\u03b3 = gamma, \\u03ba = kappa
  x <- gsub("tnfa signaling via nfkb",
            "TNF-\u03b1 signaling via NF-\u03baB", x)
  x <- gsub("wnt beta catenin",           "Wnt-\u03b2-catenin",        x)
  x <- gsub("epithelial mesenchymal transition",
            "Epithelial-mesenchymal transition",                        x)
  x <- gsub("pi3k akt mtor",              "PI3K-AKT-mTOR",             x)
  x <- gsub("il6 jak stat3",              "IL-6-JAK-STAT3",            x)
  x <- gsub("il2 stat5",                  "IL-2-STAT5",                x)
  x <- gsub("interferon alpha",           "IFN-\u03b1",                x)
  x <- gsub("interferon gamma",           "IFN-\u03b3",                x)
  x <- gsub("tgf beta",                   "TGF-\u03b2",                x)

  ## Step 4: sentence case (capitalise first character only)
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))

  ## Step 5: single-word substitutions (ignore.case so they fire whether or
  ##   not the word was capitalised by step 4)
  x <- gsub("\\bmtorc1\\b",        "mTORC1",          x, ignore.case = TRUE)
  x <- gsub("\\bmtor\\b",          "mTOR",            x, ignore.case = TRUE)
  x <- gsub("\\bp53\\b",           "p53",             x, ignore.case = TRUE)
  x <- gsub("\\bwnt\\b",           "Wnt",             x, ignore.case = TRUE)
  x <- gsub("\\bmyc\\b",           "Myc",             x, ignore.case = TRUE)
  x <- gsub("\\bkras\\b",          "KRAS",            x, ignore.case = TRUE)
  x <- gsub("\\buv\\b",            "UV",              x, ignore.case = TRUE)
  x <- gsub("\\be2f\\b",           "E2F",             x, ignore.case = TRUE)
  x <- gsub("\\bg2m\\b",           "G2/M",            x, ignore.case = TRUE)
  x <- gsub("\\bdna\\b",           "DNA",             x, ignore.case = TRUE)
  x <- gsub("\\brna\\b",           "RNA",             x, ignore.case = TRUE)
  x <- gsub("\\bnfkb\\b",          "NF-\u03baB",      x, ignore.case = TRUE)
  x <- gsub("\\btnfa\\b",          "TNF-\u03b1",      x, ignore.case = TRUE)
  x <- gsub("\\bil([0-9]+)\\b",    "IL-\\1",          x, ignore.case = TRUE)
  x <- gsub("\\bpi3k\\b",          "PI3K",            x, ignore.case = TRUE)
  x <- gsub("\\bakt\\b",           "AKT",             x, ignore.case = TRUE)
  x <- gsub("\\bjak\\b",           "JAK",             x, ignore.case = TRUE)
  x <- gsub("\\bstat([0-9]+)\\b",  "STAT\\1",         x, ignore.case = TRUE)

  x
}

#' Parse a GMT file into a data frame
#'
#' @param gmtpath Path to a GMT-format file.
#' @param limit A two-length numeric vector giving the minimum and maximum
#'   gene-set size to retain (inclusive on both ends).
#' @return A data frame with columns \code{Gene}, \code{PathwayID}, and
#'   \code{PathwayName}.
#' @export
ReadGMT <- function(gmtpath, limit = c(1, Inf)) {
  lines <- readLines(gmtpath, warn = FALSE)
  rows <- lapply(lines, function(ln) {
    fields <- strsplit(ln, "\t")[[1]]
    if (length(fields) < 3) return(NULL)
    pathway_name <- fields[1]
    genes        <- fields[-(1:2)]
    genes        <- genes[nzchar(genes)]
    data.frame(Gene        = genes,
               PathwayID   = pathway_name,
               PathwayName = pathway_name,
               stringsAsFactors = FALSE)
  })
  rows <- rows[!sapply(rows, is.null)]
  if (length(rows) == 0)
    return(data.frame(Gene = character(), PathwayID = character(),
                      PathwayName = character(), stringsAsFactors = FALSE))
  gene2path <- do.call(rbind, rows)
  ## apply size filter
  count_gene <- table(gene2path$PathwayID)
  keep <- names(count_gene)[count_gene >= limit[1] & count_gene <= limit[2]]
  gene2path[gene2path$PathwayID %in% keep, ]
}
