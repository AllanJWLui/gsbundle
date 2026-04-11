#' gsbundle: Gene Set Retrieval and Gene ID Conversion Utilities
#'
#' gsbundle provides a unified interface for downloading and caching gene sets
#' from public databases, converting gene identifiers across formats and
#' organisms, and running gene set enrichment analyses.
#'
#' @section Gene set retrieval:
#' Use [gsGetter()] to fetch gene–pathway mappings from KEGG, Reactome,
#' Gene Ontology, MSigDB, or CORUM. Results are cached locally in
#' `tools::R_user_dir("gsbundle", "cache")` and reused on subsequent calls.
#' Call [retrieve_gs()] to force a re-download, and [gs_versions()] to
#' inspect the version of each cached database.
#'
#' @section Gene ID conversion:
#' [TransGeneID()] converts a character vector of gene identifiers between
#' Entrez, Symbol, Ensembl, and UniProt formats, within a single organism or
#' across organisms (e.g. mouse → human). Cross-organism mapping uses
#' ortholog tables from MGI, NCBI HomoloGene, and Ensembl. Annotation data
#' is fetched via [getGeneAnn()] (within-organism) and [getOrtAnn()]
#' (cross-organism).
#'
#' @section Enrichment analysis:
#' [runEnrich()] wraps \pkg{clusterProfiler} to run hypergeometric
#' (`"HGT"`), over-representation (`"ORT"`), or gene set enrichment
#' (`"GSEA"`) tests against any gene set table returned by [gsGetter()].
#' Redundant pathways can be grouped or filtered using DBSCAN clustering
#' ([cluster.dbscan()]) or a pairwise Jaccard cutoff ([EnrichedFilter()]).
#'
#' @section Supported organisms:
#' Human (`hsa`), mouse (`mmu`), rat (`rno`), bovine (`bta`),
#' canine (`cfa`), chimpanzee (`ptr`), pig (`ssc`).
#'
#' @section Acknowledgements:
#' Several functions are adapted from
#' \href{https://github.com/WubingZhang/MAGeCKFlute}{MAGeCKFlute}
#' (Zhang et al., GPL >= 3) by Wubing Zhang.
#'
#' @name gsbundle-package
#' @aliases gsbundle
#' @keywords internal
"_PACKAGE"
