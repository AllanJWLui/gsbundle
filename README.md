# gsbundle

Gene set retrieval, gene ID conversion, and enrichment analysis utilities
for R.

gsbundle provides a unified interface for downloading and caching gene sets
from **KEGG**, **Reactome**, **Gene Ontology**, **MSigDB**, and **CORUM**,
converting gene identifiers across formats and organisms, and running
enrichment analyses (hypergeometric test, over-representation test, GSEA)
with optional DBSCAN-based redundancy filtering.

Annotation data is downloaded once from public databases and cached locally,
so subsequent calls are fast and work offline. The cache location is
configurable — see [Configuring the cache directory](#configuring-the-cache-directory).

## Installation

```r
# install.packages("remotes")
remotes::install_github("AllanJWLui/gsbundle")
```

### Bioconductor dependencies

Some functions require Bioconductor packages that are not installed
automatically. Install them if needed:

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# Required for cross-organism gene ID conversion via Ensembl
BiocManager::install("biomaRt")

# Required for enrichment analysis
BiocManager::install("clusterProfiler")
```

## Supported organisms

| Code  | Species              |
|-------|----------------------|
| `hsa` | Human                |
| `mmu` | Mouse                |
| `rno` | Rat                  |
| `bta` | Bovine               |
| `cfa` | Canine               |
| `ptr` | Chimpanzee           |
| `ssc` | Pig                  |

## Quick start

### Retrieve gene sets

```r
library(gsbundle)

# KEGG and Reactome pathways for human
gene2path <- gsGetter(type = "KEGG+REACTOME", organism = "hsa")
head(gene2path)
#>   Gene      PathwayID              PathwayName
#> 1 10327 KEGG_hsa00010  Glycolysis / Gluconeogenesis
#> ...

# Gene Ontology (Biological Process + Molecular Function)
go <- gsGetter(type = "GOBP+GOMF", organism = "hsa")

# All available gene sets
all_gs <- gsGetter(type = "All", organism = "hsa")

# Filter by gene set size (15–500 genes)
kegg <- gsGetter(type = "KEGG", organism = "hsa", limit = c(15, 500))

# Custom GMT file
custom <- gsGetter(gmtpath = "my_gene_sets.gmt")
```

Available `type` values and shorthands:

| Shorthand   | Expands to                                                                          |
|-------------|-------------------------------------------------------------------------------------|
| `All`       | Pathway + GO + Complex + MSigDB                                                     |
| `Pathway`   | KEGG, REACTOME, C2_CP_PID, C2_CP_BIOCARTA, C2_CP_WIKIPATHWAYS, C2_CP_KEGG_MEDICUS  |
| `GO`        | GOBP, GOCC, GOMF                                                                    |
| `Complex`   | CORUM                                                                               |
| `MSIGDB`    | C1, C2, C3, C4, C5, C6, C7, C8, H                                                  |

Combine any types with `+`, e.g. `"KEGG+GOBP+CORUM"`.

### Convert gene IDs

```r
# Symbol to Entrez (default)
TransGeneID("HLA-A", organism = "hsa")

# Symbol to UniProt
TransGeneID("HLA-A", toType = "uniprot", organism = "hsa")

# Cross-organism: mouse symbol to human symbol
TransGeneID("H2-K1", toType = "Symbol", fromOrg = "mmu", toOrg = "hsa")

# Batch conversion — returns a named vector
genes <- c("TP53", "BRCA1", "EGFR", "MYC")
TransGeneID(genes, fromType = "Symbol", toType = "Entrez", organism = "hsa")
```

Supported ID types (built-in, no biomaRt needed): `Symbol`, `Entrez`,
`Ensembl`, `UniProt`. Additional types (e.g. HGNC, RefSeq) are available
via biomaRt.

### Run enrichment analysis

Requires the `clusterProfiler` Bioconductor package.

```r
gene2path <- gsGetter(type = "KEGG+GOBP", organism = "hsa")

# Named numeric vector: names = Entrez IDs, values = scores (e.g. log2FC)
geneList <- setNames(rnorm(200), sample(gene2path$Gene, 200))

# Hypergeometric test (over-representation)
res <- runEnrich(geneList, gene2path, method = "HGT", pvalueCutoff = 0.05)
head(as.data.frame(res))

# Gene Set Enrichment Analysis
res_gsea <- runEnrich(geneList, gene2path, method = "GSEA", pvalueCutoff = 0.25)

# With DBSCAN clustering to group redundant pathways
res_clust <- runEnrich(geneList, gene2path, method = "HGT",
                       pvalueCutoff = 0.05, cluster = TRUE)

# Keep only one representative per cluster
res_filt <- runEnrich(geneList, gene2path, method = "HGT",
                      pvalueCutoff = 0.05, filter = TRUE)
```

Methods: `"HGT"` (hypergeometric test), `"ORT"` (over-representation
test), `"GSEA"` (gene set enrichment analysis).

## Key functions

| Function           | Description                                                 |
|--------------------|-------------------------------------------------------------|
| `gsGetter()`       | Retrieve gene sets from KEGG, Reactome, GO, MSigDB, CORUM  |
| `retrieve_gs()`    | Download/update gene set databases to local cache           |
| `TransGeneID()`    | Convert gene IDs across formats and organisms               |
| `getGeneAnn()`     | Fetch gene annotations from NCBI, Ensembl, and UniProt      |
| `getOrtAnn()`      | Fetch ortholog annotations (cross-organism mapping)         |
| `runEnrich()`      | Run enrichment analysis (HGT, ORT, or GSEA)                |
| `EnrichedFilter()` | Remove redundant pathways by Jaccard similarity             |
| `cluster.dbscan()` | DBSCAN clustering of enrichment results by gene overlap     |
| `ReadGMT()`        | Parse a GMT file into a data frame                          |
| `format_gs_name()` | Convert MSigDB SCREAMING_SNAKE_CASE names to readable form  |
| `gs_versions()`    | Show versions of locally cached databases                   |
| `getOrg()`         | Look up organism annotation package names                   |

## Caching and versioning

All downloaded annotations are stored in a local cache directory and reused
on subsequent calls. The cache location is resolved in this priority order:

1. **`cache.dir` argument** — passed directly to `gsGetter()`,
   `retrieve_gs()`, `getGeneAnn()`, `getOrtAnn()`, or `gs_versions()` for
   that call only
2. **`options("gsbundle.cache")`** — session-wide default; set once and all
   functions pick it up automatically
3. **System default** — `tools::R_user_dir("gsbundle", "cache")`, typically
   something like `~/.cache/R/gsbundle/` on Linux/macOS or
   `%LOCALAPPDATA%\R\Cache\R\gsbundle\` on Windows

### Configuring the cache directory

**Per-call** (highest priority):

```r
# Download to and read from a specific folder for this call only
gsGetter(type = "KEGG+GOBP", organism = "hsa",
         cache.dir = "/data/shared/gsbundle_cache")
getGeneAnn("hsa", cache.dir = "/data/shared/gsbundle_cache")
getOrtAnn("mmu", "hsa", cache.dir = "/data/shared/gsbundle_cache")
gs_versions(cache.dir = "/data/shared/gsbundle_cache")
```

**Session-wide option** (recommended for interactive use):

```r
# Set once, applies to all subsequent calls in this session
options(gsbundle.cache = "/data/shared/gsbundle_cache")

gene2path <- gsGetter(type = "KEGG+GOBP", organism = "hsa")
ann       <- getGeneAnn("hsa")
```

**Persistent default** (add to `~/.Rprofile` or a project `.Rprofile`):

```r
# ~/.Rprofile
options(gsbundle.cache = "/data/shared/gsbundle_cache")
```

### Shared / HPC cache

On a cluster or shared workstation, point all sessions at a common
pre-built cache directory so that data only needs to be downloaded once:

```r
# Populate the shared cache once (e.g. by a data manager)
retrieve_gs(type = c("KEGG", "REACTOME", "GO"), organism = "hsa",
            cache.dir = "/data/shared/gsbundle_cache")
getGeneAnn("hsa", cache.dir = "/data/shared/gsbundle_cache")

# Every user then reads from the same location
options(gsbundle.cache = "/data/shared/gsbundle_cache")
gene2path <- gsGetter(type = "KEGG+GOBP", organism = "hsa")
```

### Forcing a re-download

To refresh from source databases, pass `update = TRUE`:

```r
# Re-download KEGG and GO gene sets
gsGetter(type = "KEGG+GOBP", organism = "hsa", update = TRUE)

# Re-download gene annotations
getGeneAnn("hsa", update = TRUE)
```

### Checking cached database versions

```r
gs_versions()
#>   database organism       version  downloaded
#> 1     kegg      hsa  Release 113  2025-01-15
#> 2  ensembl      hsa          114  2025-01-15
#> 3       go      hsa  (fetched..) 2025-01-15
#> ...
```

## License

GPL-3
