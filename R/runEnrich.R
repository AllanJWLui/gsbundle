#' Run gene set enrichment analysis
#'
#' A self-contained enrichment wrapper around \pkg{clusterProfiler} that
#' operates directly on pre-fetched gene sets from \code{\link{gsGetter}},
#' without requiring \pkg{MAGeCKFlute}.
#'
#' @param geneList A named numeric vector. Names are gene IDs (Entrez by
#'   default); values are scores (e.g. log-fold changes).
#' @param gene2path A three-column data frame with columns \code{Gene},
#'   \code{PathwayID}, and \code{PathwayName}, as returned by
#'   \code{\link{gsGetter}}.
#' @param keytype Character. ID type of the names of \code{geneList}.
#'   \code{"Entrez"} (default) or \code{"Symbol"}.
#' @param method One of \code{"HGT"} (hypergeometric test, default),
#'   \code{"ORT"} (over-representation test), or \code{"GSEA"} (gene set
#'   enrichment analysis).
#' @param pvalueCutoff Adjusted p-value cutoff (default \code{0.25}).
#' @param universe A character vector specifying the background gene list.
#'   \code{NULL} (default) uses all genes in \code{gene2path} plus the input
#'   genes.
#' @param cluster Logical. If \code{TRUE}, compute DBSCAN cluster labels and
#'   add a \code{cluster} column to the result without removing any rows
#'   (default \code{FALSE}). Useful for downstream grouping (e.g.
#'   \code{slice_head(n = 1, by = cluster)}).
#' @param filter Logical. If \code{TRUE}, compute DBSCAN clusters and keep
#'   only the best-ranked (lowest p-value) representative per cluster
#'   (default \code{FALSE}). Implies \code{cluster = TRUE}.
#' @param eps Numeric. DBSCAN neighbourhood radius for the similarity-based
#'   distance matrix (default \code{1.2}). Used when \code{cluster} or
#'   \code{filter} is \code{TRUE}. Increase to merge more pathways into the
#'   same cluster.
#' @param minpts Integer. DBSCAN minimum points per cluster (default \code{2}).
#'   Used when \code{cluster} or \code{filter} is \code{TRUE}.
#' @param organism \code{'hsa'} (default) or \code{'mmu'}. Used for gene ID
#'   conversion when \code{keytype != "Entrez"}.
#' @param ... Additional arguments passed to \code{clusterProfiler::enricher}
#'   or \code{clusterProfiler::GSEA}.
#'
#' @return An \code{enrichResult} object (HGT/ORT) or \code{gseaResult}
#'   object (GSEA) from \pkg{clusterProfiler}, with an added \code{geneName}
#'   column of gene symbols. When \code{cluster = TRUE} or
#'   \code{filter = TRUE}, a \code{cluster} column is also present showing
#'   the DBSCAN cluster assignment.
#'
#' @seealso \code{\link{gsGetter}}, \code{\link{cluster.dbscan}}
#'
#' @examples
#' \dontrun{
#'   gene2path <- gsGetter(type = "KEGG+GOBP", organism = "hsa")
#'   scores    <- setNames(rnorm(200), sample(gene2path$Gene, 200))
#'   res <- runEnrich(scores, gene2path, method = "HGT", pvalueCutoff = 0.05)
#'   head(as.data.frame(res))
#'
#'   # With DBSCAN redundancy filtering
#'   # Cluster only — add column, keep all rows
#'   res_clust <- runEnrich(scores, gene2path, method = "HGT",
#'                          pvalueCutoff = 0.05, cluster = TRUE,
#'                          eps = 1.2, minpts = 2)
#'
#'   # Filter — keep one representative per cluster
#'   res_filt <- runEnrich(scores, gene2path, method = "HGT",
#'                         pvalueCutoff = 0.05, filter = TRUE,
#'                         eps = 1.2, minpts = 2)
#' }
#'
#' @export
runEnrich <- function(geneList,
                      gene2path,
                      keytype      = "Entrez",
                      method       = "HGT",
                      pvalueCutoff = 0.25,
                      universe     = NULL,
                      cluster      = FALSE,
                      filter       = FALSE,
                      eps          = 1.2,
                      minpts       = 2,
                      organism     = "hsa",
                      ...) {
  if (!requireNamespace("clusterProfiler", quietly = TRUE))
    stop("Package 'clusterProfiler' is required. Please install it.", call. = FALSE)

  method <- toupper(method)
  if (!method %in% c("HGT", "ORT", "GSEA"))
    stop("method must be one of 'HGT', 'ORT', 'GSEA'", call. = FALSE)

  ## Gene ID conversion if keytype is not Entrez
  if (keytype != "Entrez") {
    allsymbol <- names(geneList)
    gene      <- TransGeneID(allsymbol, keytype, "Entrez", organism = organism)
    idx       <- duplicated(gene) | is.na(gene)
    geneList  <- geneList[!idx]
    names(geneList) <- gene[!idx]
  }

  ## Build TERM2GENE and TERM2NAME
  TERM2GENE <- gene2path[, c("PathwayID", "Gene")]
  idx_uniq  <- !duplicated(gene2path$PathwayID)
  TERM2NAME <- gene2path[idx_uniq, c("PathwayID", "PathwayName")]

  if (method == "GSEA") {
    enrichedRes <- .enrich_GSE(geneList, TERM2GENE, TERM2NAME,
                               pvalueCutoff, organism, ...)
  } else {
    enrichedRes <- .enrich_ORT(geneList, TERM2GENE, TERM2NAME,
                               pvalueCutoff, universe, organism, ...)
  }

  ## DBSCAN clustering and/or filtering
  if (!is.null(enrichedRes) && nrow(enrichedRes@result) > 0 && (cluster || filter)) {
    res         <- enrichedRes@result
    res$cluster <- cluster.dbscan(res, eps = eps, minpts = minpts)
    if (filter) {
      ## Results are sorted by pvalue; !duplicated keeps the lowest-p row per cluster
      res <- res[!duplicated(res$cluster), ]
    }
    enrichedRes@result <- res
  }

  enrichedRes
}


## ---------------------------------------------------------------------------
## Internal: GSEA (ported from enrich.GSE.R)
## ---------------------------------------------------------------------------
.enrich_GSE <- function(geneList, TERM2GENE, TERM2NAME,
                        pvalueCutoff, organism, ...) {
  geneList <- sort(geneList, decreasing = TRUE)

  len <- length(unique(intersect(names(geneList), TERM2GENE[[2]])))
  message("\t", len, " genes are mapped ...")

  enrichedRes <- clusterProfiler::GSEA(
    geneList     = geneList,
    pvalueCutoff = pvalueCutoff,
    TERM2GENE    = TERM2GENE,
    TERM2NAME    = TERM2NAME,
    verbose      = FALSE,
    ...
  )

  if (!is.null(enrichedRes) && nrow(enrichedRes@result) > 0) {
    ## Rename core_enrichment -> geneID for consistency
    cn <- colnames(enrichedRes@result)
    if ("core_enrichment" %in% cn)
      colnames(enrichedRes@result)[cn == "core_enrichment"] <- "geneID"

    geneID   <- strsplit(enrichedRes@result$geneID, "/")
    symbols  <- TransGeneID(names(geneList), "Entrez", "Symbol", organism = organism)
    geneName <- lapply(geneID, function(gid) paste(symbols[gid], collapse = "/"))
    enrichedRes@result$geneName <- unlist(geneName)
    enrichedRes@result$Count    <- lengths(geneID)

    keep_cols <- intersect(
      c("ID", "Description", "NES", "pvalue", "p.adjust",
        "geneID", "geneName", "Count"),
      colnames(enrichedRes@result)
    )
    enrichedRes@result <- enrichedRes@result[order(enrichedRes@result$pvalue), keep_cols]
  }
  enrichedRes
}


## ---------------------------------------------------------------------------
## Internal: HGT / ORT (ported from enrich.ORT.R and enrich.R)
## ---------------------------------------------------------------------------
.enrich_ORT <- function(geneList, TERM2GENE, TERM2NAME,
                        pvalueCutoff, universe, organism, ...) {
  gene <- names(geneList)

  if (!is.null(universe)) {
    universe <- universe[!is.na(universe)]
  } else {
    universe <- unique(c(gene, TERM2GENE[[2]]))
  }

  len <- length(unique(intersect(gene, TERM2GENE[[2]])))
  message("\t", len, " genes are mapped ...")

  enrichedRes <- clusterProfiler::enricher(
    gene,
    universe     = universe,
    minGSSize    = 0,
    maxGSSize    = Inf,
    TERM2GENE    = TERM2GENE,
    TERM2NAME    = TERM2NAME,
    pvalueCutoff = pvalueCutoff,
    ...
  )

  if (!is.null(enrichedRes) && nrow(enrichedRes@result) > 0) {
    res <- enrichedRes@result[enrichedRes@result$p.adjust <= pvalueCutoff, ]
    res <- res[order(res$pvalue), ]

    ## Gene symbols
    symbols  <- TransGeneID(gene, "Entrez", "Symbol", organism = organism)
    geneID   <- strsplit(res$geneID, split = "/")
    geneName <- lapply(geneID, function(gid) paste(symbols[gid], collapse = "/"))
    res$geneName <- unlist(geneName)

    ## NES and meanVal
    res$NES <- vapply(geneID, function(gid) {
      mean(geneList[gid], na.rm = TRUE) * length(gid)^0.6
    }, numeric(1))
    res$meanVal <- vapply(geneID, function(gid) {
      mean(geneList[gid], na.rm = TRUE)
    }, numeric(1))

    ## Enrichment strength
    observed     <- vapply(strsplit(res$GeneRatio, "/"),
                           function(x) as.numeric(x[1]) / as.numeric(x[2]), numeric(1))
    expected     <- vapply(strsplit(res$BgRatio, "/"),
                           function(x) as.numeric(x[1]) / as.numeric(x[2]), numeric(1))
    res$strength <- log10(observed / expected)

    keep_cols <- c("ID", "Description", "NES", "meanVal", "strength",
                   "pvalue", "p.adjust", "GeneRatio", "BgRatio",
                   "geneID", "geneName", "Count")
    enrichedRes@result <- res[, intersect(keep_cols, colnames(res))]
  }
  enrichedRes
}


## ---------------------------------------------------------------------------
## Similarity helpers (ported from Analysis_omics_integrated_ZNF703_exploratory.qmd)
## ---------------------------------------------------------------------------

## Jaccard similarity between two gene sets
.overlap_ratio <- function(x, y) {
  x <- unlist(x)
  y <- unlist(y)
  length(intersect(x, y)) / length(unique(c(x, y)))
}

## Build a pairwise Jaccard similarity matrix from a result data frame
.get_similarity_matrix <- function(result) {
  geneSets <- setNames(strsplit(as.character(result$geneID), "/", fixed = TRUE),
                       result$ID)
  id <- result$ID
  n  <- nrow(result)
  w  <- matrix(NA_real_, nrow = n, ncol = n,
               dimnames = list(result$Description, result$Description))
  for (i in seq_len(n - 1)) {
    for (j in (i + 1):n) {
      w[i, j] <- .overlap_ratio(geneSets[id[i]], geneSets[id[j]])
    }
  }
  w[lower.tri(w)] <- t(w)[lower.tri(w)]
  diag(w) <- 1
  w
}


#' DBSCAN-based pathway clustering
#'
#' Groups enrichment result rows by pairwise gene-set similarity using DBSCAN.
#' Pathways with Jaccard similarity above \code{1 - eps} are placed in the same
#' cluster. Noise points (singletons not assigned to any cluster) each receive
#' a unique small fractional label so they are all retained when
#' deduplicating by cluster.
#'
#' @param result A data frame with columns \code{ID}, \code{Description}, and
#'   \code{geneID} (slash-separated gene IDs), as in the \code{@result} slot of
#'   a \pkg{clusterProfiler} enrichment object.
#' @param eps Numeric. DBSCAN neighbourhood radius applied to the
#'   \code{1 - Jaccard} distance matrix (default \code{1.2}). Larger values
#'   merge more pathways into the same cluster.
#' @param minpts Integer. Minimum number of points to form a dense region
#'   (default \code{2}).
#'
#' @return A numeric vector of cluster labels, one per row of \code{result}.
#'   Noise points (originally cluster 0 in DBSCAN) are assigned unique
#'   fractional labels so they are not collapsed together.
#'
#' @seealso \code{\link{runEnrich}}
#'
#' @examples
#' \dontrun{
#'   gene2path <- gsGetter(type = "KEGG", organism = "hsa")
#'   scores    <- setNames(rnorm(200), sample(gene2path$Gene, 200))
#'   res       <- runEnrich(scores, gene2path, method = "HGT")
#'   cluster.dbscan(res@result, eps = 1.2, minpts = 2)
#' }
#'
#' @export
cluster.dbscan <- function(result, eps = 1.2, minpts = 2) {
  if (!requireNamespace("dbscan", quietly = TRUE))
    stop("Package 'dbscan' is required. Please install it.", call. = FALSE)

  sim_mat  <- .get_similarity_matrix(result)
  dist_mat <- as.matrix(as.dist(1 - sim_mat))
  db_res   <- dbscan::dbscan(dist_mat, eps = eps, minPts = minpts)

  clusters <- db_res$cluster
  noise    <- which(clusters == 0)
  if (length(noise) > 0) {
    step   <- 10 ^ -ceiling(log10(length(noise)))
    clusters[noise] <- seq(step, by = step, length.out = length(noise))
  }
  clusters
}


#' Remove redundant pathways from enrichment results (Jaccard-based)
#'
#' Eliminates redundant pathways by computing the Jaccard index between gene
#' sets and discarding the lower-ranked member of any pair exceeding the
#' similarity cutoff. Ported from \pkg{MAGeCKFlute}. For DBSCAN-based
#' filtering see \code{\link{cluster.dbscan}} and the \code{filter} argument
#' of \code{\link{runEnrich}}.
#'
#' @param enrichment An \code{enrichResult} or \code{gseaResult} object, or
#'   a data frame with columns \code{ID}, \code{pvalue}, \code{NES}, and
#'   \code{geneID}.
#' @param cutoff Numeric. Jaccard index cutoff (default \code{0.8}).
#'
#' @return A filtered data frame (the \code{@result} slot content).
#'
#' @export
EnrichedFilter <- function(enrichment, cutoff = 0.8) {
  if (methods::is(enrichment, "enrichResult")) enrichment <- enrichment@result
  if (methods::is(enrichment, "gseaResult"))  enrichment <- enrichment@result
  if (nrow(enrichment) < 3) return(enrichment)

  enrichment <- enrichment[order(enrichment$pvalue, -abs(enrichment$NES)), ]
  genelist   <- strsplit(enrichment$geneID, "/")
  names(genelist) <- enrichment$ID

  tmp1 <- crossprod(table(utils::stack(genelist)))
  tmp2 <- outer(lengths(genelist), lengths(genelist), "+")
  ijc  <- tmp1 / (tmp2 - tmp1)
  diag(ijc) <- 0
  idx  <- which(ijc > cutoff, arr.ind = TRUE)
  colnames(idx) <- c("row", "col")
  idx  <- unlist(apply(idx, 1, max))
  enrichment <- enrichment[setdiff(seq_len(nrow(enrichment)), idx), ]

  genelist <- strsplit(enrichment$geneID, "/")
  names(genelist) <- enrichment$ID
  tmp1  <- crossprod(table(utils::stack(genelist)))
  ijc2  <- tmp1 / lengths(genelist)
  diag(ijc2) <- 0
  idx   <- which(ijc2 > cutoff, arr.ind = TRUE)
  colnames(idx) <- c("row", "col")
  idx   <- unlist(apply(idx, 1, max))
  enrichment <- enrichment[setdiff(seq_len(nrow(enrichment)), idx), ]

  enrichment
}
