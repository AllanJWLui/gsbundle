## Internal helper: return (and create) the gsbundle user cache directory
.gsbundle_cache <- function() {
  d <- tools::R_user_dir("gsbundle", which = "cache")
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  d
}

## Internal helper: save database version metadata to versions.json
.save_db_meta <- function(db, organism, meta) {
  if (!requireNamespace("jsonlite", quietly = TRUE))
    stop("Package 'jsonlite' is required. Please install it.", call. = FALSE)
  path <- file.path(.gsbundle_cache(), "versions.json")
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
gs_versions <- function() {
  if (!requireNamespace("jsonlite", quietly = TRUE))
    stop("Package 'jsonlite' is required. Please install it.", call. = FALSE)
  empty <- data.frame(database = character(), organism = character(),
                      version = character(), downloaded = character(),
                      source_url = character(), stringsAsFactors = FALSE)
  path <- file.path(.gsbundle_cache(), "versions.json")
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
  gene2path <- do.call(rbind, rows)
  ## apply size filter
  count_gene <- table(gene2path$PathwayID)
  keep <- names(count_gene)[count_gene >= limit[1] & count_gene <= limit[2]]
  gene2path[gene2path$PathwayID %in% keep, ]
}
