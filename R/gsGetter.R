#' Extract pathway annotation from GMT file.
#'
#' @docType methods
#' @name gsGetter
#' @rdname gsGetter
#'
#' @param gmtpath The path to customized gmt file.
#' @param type Molecular signatures for testing, available datasets include
#' Pathway (KEGG, REACTOME, C2_CP_PID, C2_CP_BIOCARTA, C2_CP_WIKIPATHWAYS, C2_CP_KEGG_MEDICUS), GO (GOBP, GOCC, GOMF),
#' MSIGDB (C1, C2 (C2_CP (C2_CP_PID, C2_CP_BIOCARTA, C2_CP_WIKIPATHWAYS, C2_CP_KEGG_MEDICUS), C2_CGP),
#' C3 (C3_MIR_MIRDB, C3_TFT_GTRD), C4 (C4_CGN, C4_CM, c4_3CA), C5 (C5_GO_BP, C5_GO_CC, C5_GO_MF, C5_HPO), C6, C7 (C7_IMMUNESIGDB, C7_VAX), C8, H)
#' and Complex (CORUM). Any combination of them are also accessible
#' (e.g. 'GOBP+GOMF+KEGG+REACTOME').
#' @param limit A two-length vector, specifying the minimal and
#' maximal size of gene sets to load.
#' @param organism 'hsa' or 'mmu'.
#' @param update Boolean, indicating whether update the gene sets from source database.
#' @param msigdb.path Path to msigdb sqlite file. Must be provided if `update` = `TRUE` and "MSIGDB" is specified in `type`
#'
#' @return A three-column data frame.
#'
#' @author Wubing Zhang
#'
#' @examples
#' gene2path = gsGetter(type = "REACTOME+KEGG")
#' head(gene2path)
#'
#' @import dplyr
#' @export
#'
gsGetter <- function(gmtpath = NULL, type = "All", limit = c(0, Inf),
                     organism = 'hsa', update = FALSE, msigdb.path = NULL){
  ## Normalize type
  type = toupper(unlist(strsplit(type, "\\+")))
  if("ALL" %in% type) type = c("PATHWAY", "GO", "COMPLEX", "MSIGDB")
  if("MSIGDB" %in% type) type = c("C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8","H", type)
  # if("GO" %in% type) type = c("C5_GO_BP", "C5_GO_CC", "C5_GO_MF", type)
  if("GO" %in% type) type = c("GOBP", "GOCC", "GOMF", type)
  if("C2" %in% type) type = c("KEGG", "REACTOME", "C2_CP_PID", "C2_CP_BIOCARTA", "C2_CP_WIKIPATHWAYS", "C2_CP_KEGG_MEDICUS", "C2_CGP", type)
  if("PATHWAY" %in% type) type = c("KEGG", "REACTOME", "C2_CP_PID", "C2_CP_BIOCARTA", "C2_CP_WIKIPATHWAYS", "C2_CP_KEGG_MEDICUS", type)
  if("COMPLEX" %in% type) type = c("CORUM", type)
  type = setdiff(type, c("ALL", "MSIGDB", "C2", "GO", "PATHWAY", "COMPLEX"))
  # type[type=="GOBP"] <- "C5_GO_BP"
  # type[type=="GOCC"] <- "C5_GO_CC"
  # type[type=="GOMF"] <- "C5_GO_MF"

  ## Update genesets
  if(update) retrieve_gs(organism=organism, type = type, msigdb.path = msigdb.path)

  ## read GMT files
  if(!is.null(gmtpath)){
    gene2path = ReadGMT(gmtpath, limit = limit)
  }else{
    gene2path = data.frame()
    if("KEGG" %in% type){
      gsfile = file.path(.gsbundle_cache(),
                         paste0("kegg.all.entrez.", organism, ".rds"))
      if(!file.exists(gsfile)) retrieve_gs(type = "KEGG", organism=organism)
      tmp = readRDS(gsfile)
      colnames(tmp) = c("ENTREZID", "PathwayID", "PathwayName")
      gene2path = rbind(gene2path, tmp)
    }
    if("REACTOME" %in% type){
      gsfile = file.path(.gsbundle_cache(),
                         paste0("reactome.all.entrez.", organism, ".rds"))
      if(!file.exists(gsfile)) retrieve_gs(type = "REACTOME", organism=organism)
      tmp = readRDS(gsfile)
      colnames(tmp) = c("ENTREZID", "PathwayID", "PathwayName")
      gene2path = rbind(gene2path, tmp)
    }
    if("CORUM" %in% type){
      gsfile = file.path(.gsbundle_cache(),
                         paste0("corum.all.entrez.", organism, ".rds"))
      if(!file.exists(gsfile)) retrieve_gs(type = "CORUM", organism=organism)
      tmp = readRDS(gsfile)
      colnames(tmp) = c("ENTREZID", "PathwayID", "PathwayName")
      gene2path = rbind(gene2path, tmp)
    }
    if(any(grepl("^GO", type))){
      gsfile = file.path(.gsbundle_cache(),
                         paste0("go.all.entrez.", organism, ".rds"))
      if(!file.exists(gsfile)) retrieve_gs(type = "GO", organism=organism)
      go = readRDS(gsfile)
      go = go[go$Category%in%gsub("GO", "", type), 1:3]
      colnames(go) = c("ENTREZID", "PathwayID", "PathwayName")
      gene2path = rbind(gene2path, go)
    }
    othertypes = setdiff(type, c("KEGG", "CORUM", "REACTOME", "GOBP", "GOMF", "GOCC"))
    if(length(othertypes)>0){
      gsfile = file.path(.gsbundle_cache(),
                         paste0("msigdb.all.entrez.", organism, ".rds"))
      if(!file.exists(gsfile)) retrieve_gs(type = "MSIGDB", organism=organism, msigdb.path = msigdb.path)
      m = readRDS(gsfile)
      colnames(m) = c("Collection","ENTREZID", "PathwayID", "PathwayName")
      m <- m %>%
        dplyr::filter(Collection %in% othertypes) %>%
        dplyr::mutate(PathwayName = format_gs_name(stringr::str_replace_all(PathwayName, "^[[:alnum:]]*_", ""))) %>%
        dplyr::select(!Collection)
      gene2path = rbind(gene2path, m)
    }
    gene2path = na.omit(gene2path)
    count_gene = table(gene2path$PathwayID)
    pathways = names(count_gene)[count_gene>limit[1] & count_gene<=limit[2]]
    gene2path = gene2path[gene2path$PathwayID%in%pathways, ]
  }
  names(gene2path) = c("Gene","PathwayID", "PathwayName")
  # gene2path$PathwayName = tolower(gsub("_", " ", gene2path$PathwayName))
  gene2path$Gene = as.character(gene2path$Gene)
  return(gene2path)
}
#' Update genesets from source database
#'
#' @docType methods
#' @name retrieve_gs
#' @rdname retrieve_gs
#'
#' @param type A vector of databases, such as KEGG, REACTOME, CORUM, GO, MSIGDB
#' @param organism 'hsa' or 'mmu'.
#' @param msigdb.path Path to msigdb sqlite file. Must be provided if `MSIGDB` is specified in `type`
#'
#' @return save data to local library.
#'
#' @author Wubing Zhang
#'
#' @import dplyr
#' @export
#'
retrieve_gs <- function(type = c("KEGG", "REACTOME", "CORUM", "GO","MSIGDB"), organism = 'hsa', msigdb.path = NULL, release = NULL){
  options(stringsAsFactors = FALSE)

  if("KEGG" %in% type){ ## Process genesets from KEGG
    if (!is.null(release))
      message("release= is ignored for KEGG — only the current version is available")
    message(format(Sys.time(), " Downloading genesets from KEGG ..."))
    gene2path = read.table(paste0("https://rest.kegg.jp/link/pathway/", organism),
                           sep = "\t", stringsAsFactors = FALSE)
    names(gene2path) = c("EntrezID", "PathwayID")
    pathways = read.table(paste0("https://rest.kegg.jp/list/pathway/", organism),
                          sep = "\t", stringsAsFactors = FALSE)
    names(pathways) = c("PathwayID","PathwayName")
    gene2path$EntrezID=gsub(".*:", "", gene2path$EntrezID)
    gene2path$PathwayID=gsub(".*:","",gene2path$PathwayID)
    pathways$PathwayID=gsub(".*:","",pathways$PathwayID)
    pathways$PathwayName=gsub(" - .*", "", pathways$PathwayName)
    rownames(pathways) = pathways$PathwayID
    gene2path$PathwayName = pathways[gene2path$PathwayID, "PathwayName"]
    locfname = file.path(.gsbundle_cache(),
                         paste0("kegg.all.entrez.", organism, ".rds"))
    gene2path$PathwayID = paste0("KEGG_", gene2path$PathwayID)
    saveRDS(gene2path, locfname)
    kegg_version <- tryCatch({
      line1 <- readLines(paste0("https://rest.kegg.jp/info/pathway/", organism), n = 1, warn = FALSE)
      m <- regmatches(line1, regexpr("Release [^\\s]+", line1))
      if (length(m) > 0) m else line1
    }, error = function(e) "(version unavailable)")
    .save_db_meta("kegg", organism, list(
      version    = kegg_version,
      downloaded = format(Sys.time(), "%Y-%m-%d"),
      source_url = paste0("https://rest.kegg.jp/link/pathway/", organism)
    ))
  }
  if("CORUM" %in% type){ ## Process genesets from CORUM
    message(format(Sys.time(), " Downloading genesets from CORUM ..."))
    corum_release <- if (!is.null(release)) paste0("releases/", release) else "releases/current"
    base_url = paste0("https://mips.helmholtz-muenchen.de/corum/download/", corum_release, "/allComplexes.txt.zip")
    corum_version <- .get_last_modified(base_url)
    locfname = file.path(.gsbundle_cache(), "allComplexes.txt.zip")
    download.file(base_url, locfname, quiet = TRUE)
    corum <- read.table(unz(locfname, "allComplexes.txt"), sep = "\t",
                        header = TRUE, quote = "", stringsAsFactors = FALSE)
    file.remove(locfname)
    if(organism=="hsa"){
      corum = corum[corum$Organism=="Human", ]
    }else if(organism=="mmu"){
      corum = corum[corum$Organism=="Mouse", ]
    }
    genes = strsplit(corum$subunits.Gene.name., ";")
    nset = unlist(lapply(genes, length))
    gene2corum = data.frame(EntrezID = unlist(genes),
                            ComplexID = rep(corum$ComplexID, nset),
                            ComplexName = rep(corum$ComplexName, nset))
    gene2corum$ComplexID = paste0("CORUM_", gene2corum$ComplexID)
    gene2corum$EntrezID = TransGeneID(gene2corum$EntrezID, "Symbol",
                                      "Entrez", organism = organism)
    gene2corum = na.omit(gene2corum)
    locfname = file.path(.gsbundle_cache(),
                         paste0("corum.all.entrez.", organism, ".rds"))
    saveRDS(gene2corum, locfname)
    .save_db_meta("corum", organism, list(
      version    = corum_version,
      downloaded = format(Sys.time(), "%Y-%m-%d"),
      source_url = base_url
    ))
  }
  if("REACTOME" %in% type){ ## Process genesets from REACTOME
    if (!is.null(release))
      message("release= is ignored for REACTOME — only the current version is available")
    message(format(Sys.time(), " Downloading genesets from REACTOME ..."))
    reactome_url <- "https://reactome.org/download/current/NCBI2Reactome.txt"
    reactome_version <- .get_last_modified(reactome_url)
    gene2path = read.table(reactome_url,
                           sep = "\t", stringsAsFactors = FALSE, comment.char = "", quote = "")
    colnames(gene2path) = c("EntrezID", "PathwayID", "link", "PathwayName", "Evidence", "Organism")
    gene2path = gene2path[grepl(organism, gene2path$PathwayID, ignore.case = TRUE), ]
    gene2path = gene2path[, c(1,2,4)]
    gene2path$PathwayID = gsub(paste0("R-", toupper(organism), "-"), "REACTOME_", gene2path$PathwayID)
    locfname = file.path(.gsbundle_cache(),
                         paste0("reactome.all.entrez.", organism, ".rds"))
    saveRDS(gene2path, locfname)
    .save_db_meta("reactome", organism, list(
      version    = reactome_version,
      downloaded = format(Sys.time(), "%Y-%m-%d"),
      source_url = reactome_url
    ))
  }
  ## Process genesets from Gene ontology
  if(any(grepl("^GO", type))){
    if (!is.null(release))
      message("release= is ignored for GO — only the current version is available")
    message(format(Sys.time(), " Downloading genesets from Gene Ontology ..."))
    go_url <- "https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2go.gz"
    go_version <- .get_last_modified(go_url)
    tmpfile = file.path(.gsbundle_cache(), "gene2go.gz")
    download.file(go_url, destfile = tmpfile, quiet = TRUE)
    go <- read.table(gzfile(tmpfile), sep = "\t", header = TRUE,
                     stringsAsFactors = FALSE, comment.char = "", quote = "")
    colnames(go) <- c("tax_id","EntrezID","GO_ID","Evidence","Qualifier","GO_term", "PubMed","Category")
    taxid = c("hsa"=9606, "mmu"=10090)
    go = go[go$tax_id==taxid[organism], ]
    go = unique(go[, c(2,3,6,8)])
    go$Category[go$Category=="Process"] = "BP"
    go$Category[go$Category=="Component"] = "CC"
    go$Category[go$Category=="Function"] = "MF"
    go$EntrezID = as.character(go$EntrezID)
    locfname = file.path(.gsbundle_cache(),
                         paste0("go.all.entrez.", organism, ".rds"))
    saveRDS(go, locfname)
    .save_db_meta("go", organism, list(
      version    = go_version,
      downloaded = format(Sys.time(), "%Y-%m-%d"),
      source_url = go_url
    ))
  }
  if("MSIGDB" %in% type){ ## Process genesets from MSigDB
    if(is.null(msigdb.path)){
      stop("Please specify msigdb.path", call. = FALSE)
    }
    if (!(requireNamespace("RSQLite", quietly = TRUE) && requireNamespace("DBI", quietly = TRUE))) {
      stop("Packages \"RSQLite\" and \"DBI\" are required. Please install it.", call. = FALSE)
    }
    message(format(Sys.time(), " Loading MSIGDB file ..."))
    msigdb <- DBI::dbConnect(RSQLite::SQLite(), msigdb.path)
    gs <- DBI::dbReadTable(msigdb,"gene_set")
    gsym <- DBI::dbReadTable(msigdb,"gene_symbol")
    msigdb.df <- DBI::dbReadTable(msigdb,"gene_set_gene_symbol")
    msigdb.df <- inner_join(msigdb.df,gs, by = c("gene_set_id" = "id"))
    msigdb.df <- inner_join(msigdb.df,gsym, by = c("gene_symbol_id" = "id"))
    msigdb.df <- data.frame(Collection = msigdb.df$collection_name,
                            EntrezID = msigdb.df$NCBI_id,
                            GeneSetID = msigdb.df$standard_name,
                            GeneSetName = msigdb.df$standard_name)
    msigdb.df$Collection <- gsub(":","_",msigdb.df$Collection)
    locfname = file.path(.gsbundle_cache(),
                         paste0("msigdb.all.entrez.", organism, ".rds"))
    saveRDS(msigdb.df, locfname)
    .save_db_meta("msigdb", organism, list(
      version    = as.character(file.mtime(msigdb.path)),
      downloaded = format(Sys.time(), "%Y-%m-%d"),
      source_url = msigdb.path
    ))
  }


}
