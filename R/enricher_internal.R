#
# From DOSE:::enricher_internal
# https://github.com/YuLab-SMU/DOSE/blob/devel/R/enricher_internal.R
#
# Author: Guangchuang Yu
# Updated: 2024-06-13
# License: Artistic 2.0 (GPL compatible)
#

##' interal method for enrichment analysis
##'
##' using the hypergeometric model
##' @title enrich.internal
##' @param gene a vector of entrez gene id.
##' @param pvalueCutoff Cutoff value of pvalue.
##' @param pAdjustMethod one of "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"
##' @param universe background genes, default is the intersection of the 'universe' with genes that have annotations.
##' Users can set `options(enrichment_force_universe = TRUE)` to force the 'universe' untouched.
##' @param minGSSize minimal size of genes annotated by Ontology term for testing.
##' @param maxGSSize maximal size of each geneSet for analyzing
##' @param qvalueCutoff cutoff of qvalue
##' @param USER_DATA ontology information
##' @return  A \code{enrichResult} instance.
##' @importClassesFrom methods data.frame
##' @importFrom qvalue qvalue
##' @importFrom methods new
##' @importFrom stats phyper
##' @importFrom stats p.adjust
##' @keywords manip
##' @author Guangchuang Yu \url{https://yulab-smu.top}
enricher_internal <- function(gene,
                              pvalueCutoff,
                              pAdjustMethod="BH",
                              universe = NULL,
                              minGSSize=10,
                              maxGSSize=500,
                              qvalueCutoff=0.2,
                              USER_DATA){

  ## query external ID to Term ID
  gene <- as.character(unique(gene))
  qExtID2TermID <- EXTID2TERMID(gene, USER_DATA)
  qTermID <- unlist(qExtID2TermID)
  if (is.null(qTermID)) {
    message("--> No gene can be mapped....")
    if (inherits(USER_DATA, "environment")) {
      p2e <- get("PATHID2EXTID", envir=USER_DATA)
      sg <- unique(unlist(p2e[1:10]))
    } else {
      sg <- unique(USER_DATA@gsid2gene$gene[1:100])
    }
    sg <- sample(sg, min(length(sg), 6))
    message("--> Expected input gene ID: ", paste0(sg, collapse=','))

    message("--> return NULL...")
    return(NULL)
  }

  ## Term ID -- query external ID association list.
  qExtID2TermID.df <- data.frame(extID=rep(names(qExtID2TermID),
                                           times=lapply(qExtID2TermID, length)),
                                 termID=qTermID)
  qExtID2TermID.df <- unique(qExtID2TermID.df)

  qTermID2ExtID <- with(qExtID2TermID.df,
                        split(as.character(extID), as.character(termID)))

  extID <- ALLEXTID(USER_DATA)
  if (missing(universe))
    universe <- NULL
  if(!is.null(universe)) {
    if (is.character(universe)) {
      force_universe <- getOption("enrichment_force_universe", FALSE)
      if (force_universe) {
        extID <- universe
      } else {
        extID <- intersect(extID, universe)
      }
    } else {
      ## https://github.com/YuLab-SMU/clusterProfiler/issues/217
      message("`universe` is not in character and will be ignored...")
    }
  }

  qTermID2ExtID <- lapply(qTermID2ExtID, intersect, extID)

  ## Term ID annotate query external ID
  qTermID <- unique(names(qTermID2ExtID))


  termID2ExtID <- TERMID2EXTID(qTermID, USER_DATA)
  termID2ExtID <- lapply(termID2ExtID, intersect, extID)

  geneSets <- termID2ExtID

  idx <- get_geneSet_index(termID2ExtID, minGSSize, maxGSSize)

  if (sum(idx) == 0) {
    msg <- paste("No gene sets have size between", minGSSize, "and", maxGSSize, "...")
    message(msg)
    message("--> return NULL...")
    return (NULL)
  }

  termID2ExtID <- termID2ExtID[idx]
  qTermID2ExtID <- qTermID2ExtID[idx]
  qTermID <- unique(names(qTermID2ExtID))

  ## prepare parameter for hypergeometric test
  k <- sapply(qTermID2ExtID, length)
  k <- k[qTermID]
  M <- sapply(termID2ExtID, length)
  M <- M[qTermID]

  N <- rep(length(extID), length(M))
  ## n <- rep(length(gene), length(M)) ## those genes that have no annotation should drop.
  n <- rep(length(qExtID2TermID), length(M))
  args.df <- data.frame(numWdrawn=k-1, ## White balls drawn
                        numW=M,        ## White balls
                        numB=N-M,      ## Black balls
                        numDrawn=n)    ## balls drawn


  ## calcute pvalues based on hypergeometric model
  pvalues <- apply(args.df, 1, function(n)
    phyper(n[1], n[2], n[3], n[4], lower.tail=FALSE)
  )

  ## gene ratio and background ratio
  #GeneRatio <- apply(data.frame(a=k, b=n), 1, function(x)
  #                   paste(x[1], "/", x[2], sep="", collapse="")
  #                   )
  #BgRatio <- apply(data.frame(a=M, b=N), 1, function(x)
  #                 paste(x[1], "/", x[2], sep="", collapse="")
  #                 )

  GeneRatio <- sprintf("%s/%s", k, n)
  BgRatio <- sprintf("%s/%s", M, N)
  RichFactor <- k / M
  FoldEnrichment <- RichFactor * N / n

  # mu and sigma are the mean and standard deviation of the hypergeometric distribution
  ## https://en.wikipedia.org/wiki/Hypergeometric_distribution
  mu <- M * n / N
  sigma <- mu * (N - n) * (N - M) / N / (N-1)
  zScore <- (k - mu)/sqrt(sigma)
  Over <- data.frame(ID = as.character(qTermID),
                     GeneRatio = GeneRatio,
                     BgRatio = BgRatio,
                     RichFactor = RichFactor,
                     FoldEnrichment = FoldEnrichment,
                     zScore = zScore,
                     pvalue = pvalues,
                     stringsAsFactors = FALSE)

  p.adj <- p.adjust(Over$pvalue, method=pAdjustMethod)
  qobj <- tryCatch(qvalue(p=Over$pvalue, lambda=0.05, pi0.method="bootstrap"), error=function(e) NULL)

  # if (class(qobj) == "qvalue") {
  if (inherits(qobj, "qvalue")) {
    qvalues <- qobj$qvalues
  } else {
    qvalues <- NA
  }

  geneID <- sapply(qTermID2ExtID, function(i) paste(i, collapse="/"))
  geneID <- geneID[qTermID]
  Over <- data.frame(Over,
                     p.adjust = p.adj,
                     qvalue = qvalues,
                     geneID = geneID,
                     Count = k,
                     stringsAsFactors = FALSE)

  Description <- TERM2NAME(qTermID, USER_DATA)

  if (length(qTermID) != length(Description)) {
    idx <- qTermID %in% names(Description)
    Over <- Over[idx,]
  }
  Over$Description <- Description
  nc <- ncol(Over)
  Over <- Over[, c(1,nc, 2:(nc-1))]


  Over <- Over[order(pvalues),]


  Over$ID <- as.character(Over$ID)
  Over$Description <- as.character(Over$Description)

  row.names(Over) <- as.character(Over$ID)

  x <- new("enrichResult",
           result         = Over,
           pvalueCutoff   = pvalueCutoff,
           pAdjustMethod  = pAdjustMethod,
           qvalueCutoff   = qvalueCutoff,
           gene           = as.character(gene),
           universe       = extID,
           geneSets       = geneSets,
           organism       = "UNKNOWN",
           keytype        = "UNKNOWN",
           ontology       = "UNKNOWN",
           readable       = FALSE
  )
  if (inherits(USER_DATA, "GSON")) {
    if (!is.null(USER_DATA@keytype)) {
      x@keytype <- USER_DATA@keytype
    }
    if (!is.null(USER_DATA@species)) {
      x@organism <- USER_DATA@species
    }
    if (!is.null(USER_DATA@gsname)) {
      x@ontology <- gsub(".*;", "", USER_DATA@gsname)
    }
  }
  return (x)
}


get_enriched <- function(object) {

  Over <- object@result

  pvalueCutoff <- object@pvalueCutoff
  if (length(pvalueCutoff) != 0) {
    ## if groupGO result, numeric(0)
    Over <- Over[ Over$pvalue <= pvalueCutoff, ]
    Over <- Over[ Over$p.adjust <= pvalueCutoff, ]
  }

  qvalueCutoff <- object@qvalueCutoff
  if (length(qvalueCutoff) != 0) {
    if (! any(is.na(Over$qvalue))) {
      if (length(qvalueCutoff) > 0)
        Over <- Over[ Over$qvalue <= qvalueCutoff, ]
    }
  }

  object@result <- Over
  return(object)
}


EXTID2TERMID <- function(gene, USER_DATA) {
  if (inherits(USER_DATA, "environment")) {
    EXTID2PATHID <- get("EXTID2PATHID", envir = USER_DATA)

    qExtID2Path <- EXTID2PATHID[gene]
  } else if (inherits(USER_DATA, "GSON")) {
    gsid2gene <- USER_DATA@gsid2gene
    qExtID2Path <- setNames(lapply(gene, function(x) {
      subset(gsid2gene, gsid2gene$gene == x)[["gsid"]]
    }), gene)
  } else {
    stop("not supported")
  }

  len <- sapply(qExtID2Path, length)
  notZero.idx <- len != 0
  qExtID2Path <- qExtID2Path[notZero.idx]

  return(qExtID2Path)
}

ALLEXTID <- function(USER_DATA) {
  if (inherits(USER_DATA, "environment")) {
    PATHID2EXTID <- get("PATHID2EXTID", envir = USER_DATA)
    res <- unique(unlist(PATHID2EXTID))
  } else if (inherits(USER_DATA, "GSON")) {
    gsid2gene <- USER_DATA@gsid2gene
    res <- unique(gsid2gene$gene)
  } else {
    stop("not supported")
  }

  return(res)
}


TERMID2EXTID <- function(term, USER_DATA) {
  if (inherits(USER_DATA, "environment")) {
    PATHID2EXTID <- get("PATHID2EXTID", envir = USER_DATA)
    res <- PATHID2EXTID[term]
  } else if (inherits(USER_DATA, "GSON")) {
    gsid2gene <- USER_DATA@gsid2gene
    res <- setNames(lapply(term, function(x) {
      subset(gsid2gene, gsid2gene$gsid == x)[["gene"]]
    }), term)
  } else {
    stop("not supported")
  }

  return(res)
}

TERM2NAME <- function(term, USER_DATA) {
  if (inherits(USER_DATA, "environment")) {
    PATHID2NAME <- get("PATHID2NAME", envir = USER_DATA)
    #if (is.null(PATHID2NAME) || is.na(PATHID2NAME)) {
    if (is.null(PATHID2NAME) || all(is.na(PATHID2NAME))) {
      return(as.character(term))
    }
    res <- PATHID2NAME[term]
    i <-  is.na(res)
    res[i] <- term[i]
  } else if (inherits(USER_DATA, "GSON")) {
    gsid2name <- USER_DATA@gsid2name
    i <- match(term, gsid2name$gsid)
    j <- !is.na(i)
    res <- term
    res[j] <- gsid2name$name[i[j]]
  } else {
    res <- as.character(term)
  }

  names(res) <- term
  return(res)
}

get_geneSet_index <- function(geneSets, minGSSize, maxGSSize) {
  if (is.na(minGSSize) || is.null(minGSSize))
    minGSSize <- 1
  if (is.na(maxGSSize) || is.null(maxGSSize))
    maxGSSize <- Inf #.Machine$integer.max

  ## index of geneSets in used.
  ## logical
  geneSet_size <- sapply(geneSets, length)
  idx <-  minGSSize <= geneSet_size & geneSet_size <= maxGSSize
  return(idx)
}

#
# From DOSE:::enricher_internal
# https://github.com/YuLab-SMU/DOSE/blob/
# 34a63655d1c24c4e855a669e61880187a28a7a1a/R/build_Anno.R#L4
#
# Author: Guangchuang Yu
# Updated: 2024-06-13
# License: Artistic 2.0 (GPL compatible)
#

##' interal method for enrichment analysis
##'
##' @param path2gene Pathway[,c("term", "gene")]# term and MetaboliteID (MetaboliteID= gene as syntax required for enricher)
##' @param path2name Pathway[,c("term", "Description")]# term and description
##' @return  A \code{enrichResult} instance.
##' @author Guangchuang Yu \url{https://yulab-smu.top}
##' @noRd
build_Anno <- function(path2gene, path2name) {
  if (!exists(".Anno_clusterProfiler_Env", envir = .GlobalEnv)) {
    pos <- 1
    envir <- as.environment(pos)
    assign(".Anno_clusterProfiler_Env", new.env(), envir = envir)
  }
  Anno_clusterProfiler_Env <- get(".Anno_clusterProfiler_Env", envir= .GlobalEnv)

  # if(class(path2gene[[2]]) == 'list') {
  if (inherits(path2gene[[2]], "list")){
    ## to compatible with tibble
    path2gene <- cbind(rep(path2gene[[1]],
                           times = vapply(path2gene[[2]], length, numeric(1))),
                       unlist(path2gene[[2]]))
  }

  path2gene <- as.data.frame(path2gene)
  path2gene <- path2gene[!is.na(path2gene[,1]), ]
  path2gene <- path2gene[!is.na(path2gene[,2]), ]
  path2gene <- unique(path2gene)

  PATHID2EXTID <- split(as.character(path2gene[,2]), as.character(path2gene[,1]))
  EXTID2PATHID <- split(as.character(path2gene[,1]), as.character(path2gene[,2]))

  assign("PATHID2EXTID", PATHID2EXTID, envir = Anno_clusterProfiler_Env)
  assign("EXTID2PATHID", EXTID2PATHID, envir = Anno_clusterProfiler_Env)

  if ( missing(path2name) || is.null(path2name) || all(is.na(path2name))) {
    assign("PATHID2NAME", NULL, envir = Anno_clusterProfiler_Env)
  } else {
    path2name <- as.data.frame(path2name)
    path2name <- path2name[!is.na(path2name[,1]), ]
    path2name <- path2name[!is.na(path2name[,2]), ]
    path2name <- unique(path2name)
    PATH2NAME <- as.character(path2name[,2])
    names(PATH2NAME) <- as.character(path2name[,1])
    assign("PATHID2NAME", PATH2NAME, envir = Anno_clusterProfiler_Env)
  }
  return(Anno_clusterProfiler_Env)
}
