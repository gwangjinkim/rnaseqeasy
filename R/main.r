


##############################################
# data
##############################################

#' Gene lengths from gtf mm10.
#'
#' symbol to gene length gtf mm10.
#'
#' @format numeric named vector: numbers are gene length name is symbol. 
#' 
#' @source \url{}
"sym2len"



######################################################################
# helper functions
######################################################################

na_to_zero <- function(vec) {
  vec[is.na(vec)] <- 0
  vec
}

####################################################################
# required packages
####################################################################

# require(DESeq2)         # DE analysis
# require(edgeR)          # DE analysis
# require(biomaRt)        # for Annotations  # sudo apt install libssl-dev
# require(org.Mm.eg.db)
# require(org.Hs.eg.db)
# require(GO.db)          # for GO analysis
# require(GOstats)        # for GO analysis
# require(gage)           # for KEGG analysis
# require(gageData)
# require(KEGG.db)
# require(annotate)
# require(genefilter)
# require(vsn)
# require(ggplot2)        # for plots
# require(Glimma)         # for interactive plots
# require(SummarizedExperiment)
# require(reshape)               # for graphics
# require(gplots)                # for heatmap
# require(RColorBrewer)          # for heatmap
# require(pheatmap)              # for heatmap
# require(gskb)                  # for gskb 
# require(clusterProfiler)
# require(GenomicFeatures)
# require(openxlsx)
# require(xlsx2dfs)

options(java.parameters = "-Xmx3000m")
# https://stackoverflow.com/questions/21937640/handling-java-lang-outofmemoryerror-when-writing-to-excel-from-r

####################################################################
# pseudoannotation
####################################################################

pseudodf <- function(vec) {
  df <- data.frame(GeneID = as.character(vec), symbol = as.character(vec), stringsAsFactors=FALSE)
  rownames(df) <- as.character(vec)
  colnames(df) <- "symbol"
  df
}

####################################################################
# non-overlapping gene lengths from gtf
# (required for DGEList obj and normalizations)
# thanks to Irsan # https://www.biostars.org/p/83901/
####################################################################

gtf2gene.length <- function(gtf.path) {
  gc()
  txdb <- GenomicFeatures::makeTxDbFromGFF(gtf.path, format = "gtf")
  exons.list.per.gene <- GenomicFeatures::exonsBy(txdb, by="gene")
  exonic.gene.sizes <- as.data.frame(sum(SummarizedExperiment::width(reduce(exons.list.per.gene))))
  colnames(exonic.gene.sizes) <- "length"
  exonic.gene.sizes
}

#####################################################################
# list flattener
#####################################################################

once.flatten <- function(obj) unlist(obj, recursive = FALSE)
k.flatten <- function(obj, k) repeated(obj, k, once.flatten)
flatten <- function(obj) unlist(obj)

#####################################################################
# select a df by list of rownames -> df.list
#####################################################################

select.by.names.vec.list <- function(data.df, list.of.names.vecs) {
  lapply(list.of.names.vecs, function(names.vec) data.df[names.vec, ])
}

#####################################################################
# Return for a df, its sortings by l2FC, p-val and rownames
#####################################################################

df.sortings <- function(DE.res.df) {
  list(l2FC.srt = DE.res.df[order(DE.res.df$log2FoldChange, decreasing = TRUE), ],
       p.srt = DE.res.df[order(DE.res.df$pvalue, decreasing = FALSE), ],
       names.srt = DE.res.df[order(rownames(DE.res.df), decreasing = FALSE), ])
}

#####################################################################
# print DE sortings
#####################################################################

print.DE.sortings <- function(DE.list, fname, dir = ".") {
  DE.list.with.sortings <- lapply(once.flatten(lapply(DE.list, df.sortings)), as.data.frame)
  xlsx2dfs::dfs2xlsx(DE.list.with.sortings, file.path(dir, fname))
}

#####################################################################
# print cnts selections of DE 
#####################################################################

print.cnts.DE.sortings <- function(cnts, DE.list, fname, dir = ".") {
  DE.list.with.sortings <- lapply(once.flatten(lapply(DE.list, df.sortings)),
                                  as.data.frame)
  DE.list.with.sortings.names <- lapply(DE.list.with.sortings, rownames)
  DE.cnts.list.with.sortings <- select.by.names.vec.list(as.data.frame(cnts),
                                                         DE.list.with.sortings.names)
  xlsx2dfs::dfs2xlsx(DE.cnts.list.with.sortings, file.path(dir, fname))
}

####################################################################
# helper functions for timestamp
####################################################################

time.now <- function(time_now = Sys.time()) format(time_now, "%y%m%d%H%M%S")

#####################################################################
# helper functions for averaging tables
#####################################################################

counts.avg.build <- function(cnts.df, cols.list, groups){
  cnts.df <- as.data.frame(cnts.df)
  ncol_old <- length(colnames(cnts.df))
  ncol_limit <- length(groups) + ncol_old
  new_col_names <- c(colnames(cnts.df), groups)
  cnts.df <- cbind(cnts.df,
                   lapply(cols.list,
                          function(idxs) rowMeans(cnts.df[, idxs, drop = FALSE])))
  colnames(cnts.df) <- new_col_names
  cnts.df[, (ncol_old + 1):ncol_limit]
}

counts.std.build <- function(cnts.df, cols.list, groups){
  cnts.df <- as.data.frame(cnts.df)
  rowSds <- function(df) apply(df, 1, sd)
  ncol_old <- length(colnames(cnts.df))
  ncol_limit <- length(groups) + ncol_old
  new_col_names <- c(colnames(cnts.df), groups)
  cnts.df <- cbind(cnts.df,
                   lapply(cols.list,
                          function(idxs) rowSds(cnts.df[, idxs, drop = FALSE])))
  colnames(cnts.df) <- new_col_names
  cnts.df[, (ncol_old + 1):ncol_limit]
}

# ####################################################################
# Read-in meta.df and indirpath to count-matrix
# ####################################################################

# the point is that count files have different genes
# while I want all count files to have same rows.
# now one can unify all names, make them unique,
# then select from all data frames,
# but there will be all-NA rows
# all NA values should become 0!

read.tab <- function(fpath) {
  read.delim(fpath, sep = "\t", head = F, row.names = 1, stringsAsFactors = F)
}

read.dfs2table <- function(meta.df, indirpath) {
    # reads in all meta files
    fpaths <- file.path(indirpath, as.character(meta.df$fileName))
    dfs <- lapply(fpaths, read.tab)
    
    # collect names and make them unique and sort them
    names_list <- lapply(dfs, rownames)
    unified_names <- sort(unique(unlist(names_list)))

    # select from each df these unique sorted names and make NA to 0
    dfs <- lapply(dfs, function(df) {
        df <- df[unified_names, , drop=FALSE]      # this drop = FALSE is very mean!
        rownames(df) <- unified_names              # rename rows # NA NA.1 etc.
        df[is.na(df)] <- 0                         # make all NA values to 0
        df
    })
    
    # bind them to one data frame
    res_df <- Reduce(cbind, dfs)
    # and give them sample names
    names(res_df) <- meta.df$sampleName
    res_df
}


# ####################################################################
# Translate vector values
# ####################################################################

translate.vec <- function(values.vec, from.vec, to.vec) {
  tr.vec <- to.vec
  names(tr.vec) <- from.vec
  tr.vec[values.vec]
}

# ####################################################################
# cluster df sorting helper functions
# ####################################################################

# df to list and back ################################################

select_df <- function(df, val, col.selector) {
  df[ df[, col.selector] == val, ]
}

df2dflist <- function(df, col.selector) { # actually it is split()
  col.vals <- unique(df[, col.selector])
  dfl <- lapply(seq(col.vals),
                function(i) select_df(df,
                                      val = col.vals[i],
                                      col.selector))
  names(dfl) <- col.vals
  dfl
}

dflist2df <- function(dfl) {
  Reduce(rbind, dfl)
}

# ordering df lists  #################################################

order.list.by.col.mean <- function(dfl, col.selector) {
  dfl[order(unlist(lapply(dfl, function(df) mean(df[, col.selector]))))]
}

order.list.by.df.function <- function(dfl, df.function, ...) {
  dfl[order(unlist(lapply(dfl, function(df) df.function(df, ...))))]
}

list.order.by.col.mean <- function(dfl, col.selector) {
  order(unlist(lapply(dfl, function(df) mean(df[, col.selector]))))
}

list.order.manually.by.vec <- function(dfl, vec) {
  dfl[vec]
}

repeated.values.into.df.list <- function(dfl, col.selector, vals) {
  Map(function(df, val) {df[, col.selector] <- val; df}, dfl, vals)
}

list.order.by.col.mean.diffs <- function(dfl, col.selector.one, 
                                         col.selector.two, 
                                         decreasing = FALSE) {
    dfl[order(unlist(lapply(dfl, function(df) mean(df[, col.selector.one]))) -
                unlist(lapply(dfl, function(df) mean(df[, col.selector.two]))),
              decreasing = decreasing)]
}

#######################################################################
# create a  DESeq2 obj out of
# - path to counts
# - file names in meta.txt
# - (outputpath)
#######################################################################

meta2DESeq2.obj <- function(meta.df, indirpath, normalized=FALSE) {
  # require(DESeq2)
  count.matrix <- read.dfs2table(meta.df, indirpath)
  # if dfs have differing column values, then fill the other rows with zero!
  DESeq2.obj <- DESeq2::DESeqDataSetFromMatrix(
    countData = count.matrix,
    colData = meta.df,
    design = ~ 0 + condition)
  if (normalized) {
    DESeq2.obj <- DESeq2::estimateSizeFactors(DESeq2.obj)
    DESeq2.obj <- DESeq2::estimateDispersions(DESeq2.obj)
  }
  DESeq2.obj
}

#######################################################################
# create a raw count table out of
# - path to counts
# - file names in meta.txt
# - (outputpath)
# create a DESeq2-normalized table out of
# - path to counts
# - file names in meta.txt
# - (outputpath)
# create an averaged DESeq2-normalized table out of
# - path to counts
# - file names in meta.txt
# - (outputpath)
#######################################################################
denom <- function(meta.df) as.character(meta.df$condition[meta.df$testing == "denom"][1])
num   <- function(meta.df) as.character(meta.df$condition[meta.df$testing == "num"][1])
core.name <- function(meta.df) paste0(num(meta.df), "-vs-", denom(meta.df), collapse = '')

meta2cnts <- function(meta.df, DESeq2.obj, outdirpath=".",
                      dataname = dataname,
                      printp=FALSE, normalized=FALSE, averaged=FALSE,
                      sheetName) {
  cnts <- DESeq2::counts(DESeq2.obj, normalized=normalized)
  if (averaged) {
    cond <- unique(as.character(meta.df$condition))
    cond.cols <- lapply(cond,
                        function(cd) which(as.character(meta.df$condition) == cd))
    names(cond.cols) <- cond
    cnts.avg    <- counts.avg.build(cnts, cond.cols, cond)
    cnts.sd <- counts.std.build(cnts, cond.cols, cond)

    res <- xlsx2dfs::withNames(ifelse(normalized, "nrm-counts-avg", "raw-counts-avg"), cnts.avg,
                     ifelse(normalized, "nrm-counts-sd", "raw-counts-sd"), cnts.sd)
    if (!printp) {
      return(res)
    }
  }
  if (printp) {
    if (averaged) {
      filename <- paste0(ifelse(normalized, "nrm-counts-", "raw-counts-"),
                         "avg-",
                         dataname, "-", core.name(meta.df), ".xlsx")
      xlsx2dfs::dfs2xlsx(res,
                file.path(outdirpath, filename))
      return(res)
    }
    filename <- paste0(ifelse(normalized, "nrm-counts-", "raw-counts-"),
                       ifelse(averaged, "avg-", ""),
                       dataname, "-", core.name(meta.df), ".xlsx")
    xlsx2dfs::dfs2xlsx(xlsx2dfs::withNames(ifelse(normalized, "nrm-counts", "raw-counts"), cnts),
              file.path(outdirpath, filename))
  }
  cnts
}

#######################################################################
# create list of differentially expressed genes
# create list of differentially upregulated genes
# create list of differentially downregulated genes
# print them out
# out of
# path to counts
# file name in meta.txt - metapath
# also for single-rep RNAseq analysis! ("DESeq2")
#######################################################################

DEanalysis <- function(meta.df, DESeq2.obj.disp, outdirpath=".", dataname="",
                       printp=FALSE, prior.count=0,
                       alpha=0.05, lFC=1, filterp=FALSE) {
  dds <- DESeq2::DESeq(DESeq2.obj.disp) ## changed from DESeq2.obj at 2020.09.10 16:27
  res <- DESeq2::results(dds, contrast=c("condition", num(meta.df), denom(meta.df)),
                 cooksCutoff = Inf,
                 independentFiltering = FALSE) # those two avoid NA!
  if (filterp) {
    res <- subset(subset(res, padj < alpha), abs(log2FoldChange) > lFC)
  }
  
  up <- subset(res, res$log2FoldChange > 0)
  down <- subset(res, res$log2FoldChange < 0)
  res <- list(all=res, up=up, down=down)
  
  if (printp) {
    filecorename <- paste0("DE_", ifelse(filterp, "sig_", ""),
                           dataname, "_", core.name(meta.df), "_", time.now())
    if (filterp) { filecorename <- paste0(filecorename, "_", alpha, "_", lFC) }
    filename <- paste0(filecorename, ".xlsx")
    print.DE.sortings(res, fname = filename, dir = outdirpath)
  }
  res
}

#######################################################################
# create an averaged heatmap and group-pictures out of
# - k
# - (DESeq2 obj OR DEGList obj OR raw count table OR cpm-normalized table OR
#   DESeq2-normalized table)
#   metapath and indirpath
# - (outputpath)
#######################################################################

meta2heatmap <- function(meta.df, cnts.avg.nrm, resSig, outdirpath=".",
                         dataname = dataname, selected.genes = NULL, name.add = "", # for showing in name
                         k = k, printp=FALSE,
                         alpha = alpha, lFC= lFC, filterp=FALSE,
                         xlsxp=TRUE, csvp=FALSE, tsvp=FALSE) {
                         
    res.names <- rownames(resSig$all)

    ## if sth given in 'selected.genes': if "up" or "down" use them from resSig, else the vector given:
    if (length(selected.genes) > 0) {
        if (selected.genes == "up") {
            res.names <- rownames(resSig$up)
        } else if (selected.genes == "down") {
            res.names <- rownames(resSig$down)
        } else {
            res.names <- selected.genes
        }
    }

    cnts.res <- cnts.avg.nrm[res.names, ]

    # gene-wise normalization
    scaledata <- t(scale(t(cnts.res)))
    scaledata <- scaledata[complete.cases(scaledata), ]

    # k means clustering
    kClust <- kmeans(scaledata, centers = k, nstart = 1000, iter.max = 30)
    kClusters <- kClust$cluster

    # function to find centroid (cluster core) in cluster i
    clust.centroid <- function(i, dat, clusters) {
        ind = (clusters == i)
        colMeans(dat[ind, ])
    }
    kClustcentroids <- sapply(levels(factor(kClusters)),
                            clust.centroid, scaledata, kClusters) ## is a matrix

    # plot centroids
    Kmolten <- reshape::melt(kClustcentroids)
    colnames(Kmolten) <- c("sample", "cluster", "value")
    # ensure correct factorizing
    Kmolten$sample <- factor(Kmolten$sample,
                           levels = unique(meta.df$condition))
    {
        p1 <- ggplot2::ggplot(Kmolten, ggplot2::aes(x=factor(sample, levels = unique(meta.df$condition)),
                                  y=value, group = cluster,
                                  colour=as.factor(cluster))) +
              ggplot2::geom_point() +
              ggplot2::geom_line() +
              ggplot2::xlab("Time") +
              ggplot2::ylab("Expression") +
              ggplot2::labs(title = "Cluster Expression by Group", color = "Cluster") +
              ggplot2::theme(text=ggplot2::element_text(size=16)) # added becaus of fonts problem
        out_fpath <- paste0(outdirpath, "/k", k, "_", core.name(meta.df), paste0("-ClusterAll_", alpha, "_", lFC, "_", time.now(), ".svg"))
        ggplot2::ggsave(plot=p1, filename=out_fpath, dpi=300)
    }

    # check similarity of centroids
    print(cor(kClustcentroids))


    clusters_unique <- sort(unique(Kmolten$cluster))
    cores <- list()
    Kmolten.list <- list()
    scores <- list()
    corescores <- list()

    for (i in clusters_unique) {
      core_i <- Kmolten[Kmolten$cluster == i, ]
      core_i$sample <- factor(core_i$sample, levels = unique(meta.df$condition))
      
      # get clusters
      K_i <- scaledata[kClusters == i, ]
      
      # calculate correlation with core
      corescore_i <- function(x) cor(x, core_i$value)
      score_i <- apply(K_i, 1, corescore_i)
      
      # get df into long format for plotting
      K_i_molten <- reshape::melt(K_i)
      colnames(K_i_molten) <- c("gene", "sample", "value")
      
      # add the score
      K_i_molten <- merge(K_i_molten, score_i, by.x = "gene", by.y = "row.names", all.x = T)
      colnames(K_i_molten) <- c('gene', 'sample', 'value', 'score')
      
      # order
      K_i_molten <- K_i_molten[order(K_i_molten$score), ]
      
      # collect all results
      Kmolten.list[[i]] <- K_i_molten
      scores[[i]] <- score_i
      corescores[[i]] <- corescore_i
      cores[[i]] <- core_i
      
      # plot cluster in groups

      sp_i <- ggplot2::ggplot(K_i_molten, 
                             ggplot2::aes(x=factor(sample, levels=unique(meta.df$condition)), y=value)) + 
             ggplot2::geom_line(ggplot2::aes(colour=score, group=gene)) + 
             ggplot2::scale_color_gradientn(colours=c('blue1', 'red2')) + 
             ggplot2::geom_line(data=core_i,  
                                ggplot2::aes(sample,value,group=cluster), 
                                color='black', 
                                inherit.aes=FALSE) + # adding score
             ggplot2::xlab('Time') + 
             ggplot2::ylab('Expression') + 
             ggplot2::labs(title=paste0('Cluster ', i, ' Expression by Group'), color = 'Score')
      fpath <- file.path(outdirpath, paste0("/k", k, "_", core.name(meta.df), "-Cluster", i, "_", dataname, "_", alpha, "_", lFC, name.add, ".svg"))
      ggplot2::ggsave(plot=sp_i, filename=fpath, dpi=300)
    }
    
    # name collected values
    names(Kmolten.list) <- paste("K", 1:k, "molten", sep = '')
    names(cores) <- paste("core", 1:k, sep = '')
    names(scores) <- paste("score", 1:k, sep = '')

    # prepare heatmap
    colors.kr <- colorRampPalette(c("black", "red"))(100)
    #  colors.kgr <- colorRampPalette(c("black", "grey88", "red"))(100)

    # add cluster number and score for sorting the data
    scaledata.k <-cbind(scaledata,
                        cluster = kClust$cluster,
                        score = Reduce(`c`, scores)[rownames(scaledata)])
    scaledata.k <- scaledata.k[order(scaledata.k[, "cluster"], scaledata.k[, "score"]), ]

    scaledata_list <- df2dflist(scaledata.k, "cluster")
    names(scaledata_list) <- paste("cluster", 1:k, sep = "")
    # outpath
    outname <- file.path(outdirpath, paste0("k", k, "_khmap_", dataname, "_", time.now(), "_UO_", alpha, "_", lFC, name.add))

    # for gaps
    gaps.idxs <- cumsum(table(scaledata.k[, "cluster"])) # scaledata.k is one matrix! only column 'cluster'
    # dim(scaledata.k)
    # unordered clusters

    ##################
    # black to red
    ##################

    {
        svg(paste0(outname, ".svg"))
            pheatmap::pheatmap(scaledata.k[, 1:length(unique(meta.df$condition))],
                     cluster_rows = F,
                     cluster_cols = F,
                     cellwidth = 40,
                     col = colors.kr,
                     fontsize_row = 0.5,
                     border_color = NA,
                     gaps_row = gaps.idxs          # gap after each block
            )
        dev.off()
    }

    {
        setEPS()
        postscript(paste0(outname, ".eps"))
            pheatmap::pheatmap(scaledata.k[, 1:length(unique(meta.df$condition))],
                     cluster_rows = F,
                     cluster_cols = F,
                     cellwidth = 40,
                     col = colors.kr,
                     fontsize_row = 0.5,
                     border_color = NA,
                     gaps_row = gaps.idxs          # gap after each block
            )
        dev.off()
    }

  # print scaledata
  outfname <- paste0(outname, '.xlsx')
  xlsx2dfs::dfs2xlsx(scaledata_list, outfname)
  
  # save all inbetween results
  data <- list(scaledata_list, Kmolten.list, cores, scores)
  names(data) <- c("scaledata_list", "Kmolten.list", "core.list", "score.list")
  saveRDS(data, file = paste0(outdirpath, "/scaledata.Kmolten.core.score.list.", dataname, ".", time.now(), ".", name.add, ".rds"))
  data
}



# meta2heatmap <- function(meta.df, cnts.avg.nrm, resSig, outdirpath=".",
#                          dataname = dataname, selected.genes = NULL, name.add = "", # for showing in name
#                          k = k, printp=FALSE,
#                          alpha = alpha, lFC= lFC, filterp=FALSE,
#                          xlsxp=TRUE, csvp=FALSE, tsvp=FALSE) {
#   res.names <- rownames(resSig$all)
#   
#   ## if sth given in 'selected.genes': if "up" or "down" use them from resSig, else the vector given:
#   if (length(selected.genes) > 0) {
#     if (selected.genes == "up") {
#       res.names <- rownames(resSig$up)
#     } else if (selected.genes == "down") {
#       res.names <- rownames(resSig$down)
#     } else {
#       res.names <- selected.genes
#     }
#   }
#   
#   cnts.res <- cnts.avg.nrm[res.names, ]
#   
#   # gene-wise normalization
#   scaledata <- t(scale(t(cnts.res)))
#   scaledata <- scaledata[complete.cases(scaledata), ]
#   
#   # k means clustering
#   kClust <- kmeans(scaledata, centers = k, nstart = 1000, iter.max = 30)
#   kClusters <- kClust$cluster
#   
#   # function to find centroid (cluster core) in cluster i
#   clust.centroid <- function(i, dat, clusters) {
#     ind = (clusters == i)
#     colMeans(dat[ind, ])
#   }
#   kClustcentroids <- sapply(levels(factor(kClusters)),
#                             clust.centroid, scaledata, kClusters) ## is a matrix
#   
#   # plot centroids
#   Kmolten <- reshape::melt(kClustcentroids)
#   colnames(Kmolten) <- c("sample", "cluster", "value")
#   # ensure correct factorizing
#   Kmolten$sample <- factor(Kmolten$sample,
#                            levels = unique(meta.df$condition))
#   {
#     p1 <- ggplot2::ggplot(Kmolten, ggplot2::aes(x=factor(sample, levels = unique(meta.df$condition)),
#                               y=value, group = cluster,
#                               colour=as.factor(cluster))) +
#       ggplot2::geom_point() +
#       ggplot2::geom_line() +
#       ggplot2::xlab("Time") +
#       ggplot2::ylab("Expression") +
#       ggplot2::labs(title = "Cluster Expression by Group", color = "Cluster") +
#       ggplot2::theme(text=ggplot2::element_text(size=16))
#     out_fpath <- file.path(outdirpath, paste0("k", k, "_", core.name(meta.df), "-ClusterAll_", alpha, "_", lFC, "_", time.now(), ".png"))
#     ggplot2::ggsave(plot=p1, filename=out_fpath, dpi=300)
#   }
#   
#   # check similarity of centroids
#   print(cor(kClustcentroids))
#   
#   for (i in 1:k) {
#     # subset cores molten df to plot core
#     assign(paste0("core", i), Kmolten[Kmolten$cluster == i, ])
#     eval(parse(text = paste0("core", i, "$sample <- factor(core", i,
#                              "$sample, levels = unique(meta.df$condition))")))
#     # get clusters
#     assign(paste0("K", i), scaledata[kClusters == i, ])
#     
#     # calculate correlation with core
#     assign(paste0("corescore", i),
#            eval(parse(text = paste0("function(x) {cor(x, core", i, "$value)}"))))
#     assign(paste0("score", i),
#            eval(parse(text = paste0("apply(K", i, ", 1, corescore", i, ")"))))
#     
#     # get data frame into long format for plotting
#     assign(paste0("K", i, "molten"),
#            eval(parse(text = paste0("reshape::melt(K", i, ")"))))
#     eval(parse(text = paste0("colnames(K", i, "molten) <- c('gene', 'sample', 'value')")))
#     
#     # add the score
#     eval(parse(text = paste0("K", i, "molten <- merge(K", i, "molten, score", i, ", by.x = 'gene', by.y = 'row.names', all.x = T)")))
#     eval(parse(text = paste0("colnames(K", i, "molten) <- c('gene', 'sample', 'value', 'score')")))
#     
#     # order
#     eval(parse(text = paste0("K", i, "molten <- K", i, "molten[order(K", i, "molten$score), ]")))
#   }
#   
#   # plot cluster groups
#   for (i in 1:k) {
#     text = paste0("sp", i, " <- ggplot2::ggplot(K", i, "molten, ggplot2::aes(x=factor(sample,",
#                   " levels=unique(meta.df$condition)), y=value)) + ",
#                   "ggplot2::geom_line(ggplot2::aes(colour=score, group=gene)) + ",
#                   "ggplot2::scale_color_gradientn(colours=c('blue1', 'red2')) + ",
#                   # this adds the core
#                   "ggplot2::geom_line(data=core", i, ", ggplot2::aes(sample,value,group=cluster), ",
#                   "color='black', inherit.aes=FALSE) + ",
#                   "ggplot2::xlab('Time') + ",
#                   "ggplot2::ylab('Expression') + ",
#                   "ggplot2::labs(title='Cluster ", i, " Expression by Group', color = 'Score'); ",
#                   "fpath <- '", outdirpath, "/k", k, "_", core.name(meta.df), "-Cluster", i, "_", dataname, "_", alpha, "_", lFC, name.add, ".png'; ",
#                   "ggplot2::ggsave(plot=sp", i, ", filename=fpath, dpi=300)"
#     )
#     eval(parse(text = text))
#   }
#   
#   # prepare heatmap
#   colors.kr <- colorRampPalette(c("black", "red"))(100)
# #  colors.kgr <- colorRampPalette(c("black", "grey88", "red"))(100)
#   
#   eval(parse(text = paste0("scores <- c(", paste(paste("score", 1:k, sep = ''),
#                                                  collapse = ", "), ")")))
#   # add cluster number and score for sorting the data
#   scaledata.k <-cbind(scaledata,
#                       cluster = kClust$cluster,
#                       score = scores[rownames(scaledata)])
#   scaledata.k <- scaledata.k[order(scaledata.k[, "cluster"], scaledata.k[, "score"]), ]
#   
#   # outpath
#   outname <- file.path(outdirpath, paste0("k", k, "_khmap_", dataname, "_", time.now(), "_UO_", alpha, "_", lFC, name.add))
#   
#   # for gaps
#   gaps.idxs <- cumsum(table(scaledata.k[, "cluster"])) # scaledata.k is one matrix! only column 'cluster'
#   # dim(scaledata.k)
#   # unordered clusters
#   
#   ##################
#   # black to red
#   ##################
#   
#   {
#     svg(paste0(outname, ".svg"))
#     pheatmap::pheatmap(scaledata.k[, 1:length(unique(meta.df$condition))],
#              cluster_rows = F,
#              cluster_cols = F,
#              cellwidth = 40,
#              col = colors.kr,
#              fontsize_row = 0.5,
#              border_color = NA,
#              gaps_row = gaps.idxs          # gap after each block
#     )
#     dev.off()
#   }
#   
#   {
#     setEPS()
#     postscript(paste0(outname, ".eps"))
#     pheatmap::pheatmap(scaledata.k[, 1:length(unique(meta.df$condition))],
#              cluster_rows = F,
#              cluster_cols = F,
#              cellwidth = 40,
#              col = colors.kr,
#              fontsize_row = 0.5,
#              border_color = NA,
#              gaps_row = gaps.idxs          # gap after each block
#     )
#     dev.off()
#   }
#   
# 
# 
# 
# 
#   # print scaledata
#   scaledata_list <- df2dflist(scaledata.k, "cluster")
#   outfname <- paste0(outname, '.xlsx')
#   names(scaledata_list) <- paste("cluster", 1:k, sep = "")
#   xlsx2dfs::dfs2xlsx(scaledata_list, outfname)
#   
#   eval(parse(text = paste0("Kmolten.list <- list(", paste(paste("K", 1:k, "molten", sep = ''), collapse = ", "), ")")))
#   names(Kmolten.list) <- paste("K", 1:k, "molten", sep = '')
#   
#   eval(parse(text = paste0("core.list <- list(", paste(paste("core", 1:k, sep = ''), collapse = ", "), ")")))
#   names(Kmolten.list) <- paste("core", 1:k, sep = '')
#   
#   eval(parse(text = paste0("score.list <- list(", paste(paste("score", 1:k, sep = ''), collapse = ", "), ")")))
#   names(Kmolten.list) <- paste("score", 1:k, sep = '')
#   
#   data <- list(scaledata_list, Kmolten.list, core.list, score.list)
#   names(data) <- c("scaledata_list", "Kmolten.list", "core.list", "score.list")
#   saveRDS(data, file = paste0(outdirpath, "/scaledata.Kmolten.core.score.list.", dataname, ".", time.now(), ".", name.add, ".rds"))
#   data
# }
# 
# 

#######################################################################
# MDS plot, iMDS plot
# - metapath
# - indirpath
# - lFC, alpha
# - method
# - outdirpath
#######################################################################
# BiocManager::install("DEFormats")
# meta2iMDSplot <- function(meta.df, DESeq2.obj.disp, outdirpath=".",
#                           dataname = dataname, top=500,
#                           launch = TRUE) {
#   title <- paste0("MDS Plot ", dataname, " ", time.now())
#   filename <- paste0("iMDS_", dataname, "_", core.name(meta.df), "_", 
#                      time.now(), "_DESeq2", collapse = "_")
#   Glimma::glMDSPlot(DEFormats::as.DGEList(DESeq2.obj.disp), 
#             ## changed from DESeq2.obj at 2020.09.10 16:27
#             ## and added DEFormats::as.DGEList
#             top = top,
#             path = outdirpath,
#             main = title,
#             html = filename,
#             launch = launch
#   )
# }
# one could contribute and write Glimma::glMDSPlot.DESeqDataSet

# corrected verion for correct coloring
meta2iMDSplot <- function(meta.df, 
                          DESeq2.obj.disp, 
                          outdirpath=".",
                          dataname = dataname, 
                          top=500,
                          launch = TRUE) {
  title <- paste0("MDS Plot ", dataname, " ", time.now())
  filename <- paste0("iMDS_", dataname, "_", core.name(meta.df), "_",
                     time.now(), "_DESeq2", collapse = "_")
  labels <- meta.df$condition
  Glimma::glMDSPlot(DEFormats::as.DGEList(DESeq2.obj.disp),
            top = top,
            labels = labels, # this cures the labeling problem!
            groups = labels, # this cures the coloring problem!
            path = outdirpath,
            main = title,
            html = filename,
            launch = launch
  )
}


#######################################################################
# volcano plot, iVolcano Plot, MD plot, iMD plot
# - metapath
# - indirpath
# - lFC, alpha
# - method
# - outdirpath
#######################################################################

meta2iVolcano <- function(meta.df, DESeq2.obj, DESeq2.obj.disp, outdirpath=".",
                          dataname= dataname,
                          lFC=lFC, alpha=alpha, launch = TRUE) {
  wald.test <- DESeq2::nbinomWaldTest(DESeq2.obj.disp)
  res.DESeq2 <- DESeq2::results(wald.test, alpha=alpha, pAdjustMethod = "BH",
                        contrast = c("condition", num(meta.df), denom(meta.df)),
                        cooksCutoff = Inf,
                        independentFiltering = FALSE) # avoid NA in padj
  title <- paste0("Volcano Plot ", dataname, " ", time.now())
  filename <- paste0("iVolcano_", dataname, "_", core.name(meta.df), "_", time.now(), "_", alpha, "_", lFC, "_DESeq2")
  {
    Glimma::glXYPlot(x=res.DESeq2$log2FoldChange, y=-log10(res.DESeq2$pvalue),
             counts = DESeq2::counts(DESeq2.obj)[rownames(res.DESeq2), ],
             anno = pseudodf(rownames(DESeq2::counts(DESeq2.obj))),
             groups = meta.df$condition,
             samples = meta.df$sampleName,
             xlab = "log2FC",
             ylab = "log10padj",
             main = title,
             status = na_to_zero(ifelse(abs(res.DESeq2$log2FoldChange) > lFC &
                               res.DESeq2$padj < alpha, 1, 0)), # na_to_zero() had to be done
             side.main = "symbol",
             side.xlab = "Group",
             side.ylab = "Counts",
             path = outdirpath,
             html = filename,
             launch = launch)
  }
  title <- paste0("iMD Plot ", dataname, " ", time.now())
  filename <- paste0("iMD_", dataname, "_", core.name(meta.df), "_", time.now(), "_DESeq2")
  {
    Glimma::glMDPlot(x = res.DESeq2,
             counts = DESeq2::counts(DESeq2.obj)[rownames(res.DESeq2), ],
             anno = pseudodf(rownames(DESeq2::counts(DESeq2.obj))), # GeneID and symbol as col
             groups = factor(meta.df$condition, levels=unique(meta.df$condition)),
             samples = meta.df$sampleName,
             ylab = "log2FC",
             xlab = "Average log10 CPM",
             main = title,
             status = ifelse(abs(res.DESeq2$log2FoldChange) > lFC &
                               res.DESeq2$padj < alpha, 1, 0), # maybe in future na_to_zero() necessary?
             side.xlab = "Group",
             side.ylab = "Counts",
             side.main = "symbol",
             path = outdirpath,
             html = filename,
             launch = launch)
  }
}



#####################################################################
# GOI
#####################################################################

bra.down <- c("Msgn1", "Osr1", "Rspo3", "Fgf8", "Wnt3a")
eo.down.repr <- c("Dmdx1", "Dpf3", "Foxa1", "Hey1", "Hhex", "Tcf7l2", "Tle2")
eo.down <- c("Mesp1", "Foxa2", "Sox17", "Lhx1", "Cer1")
pp.up <- c("Lefty1", "Lefty2", "Nodal", "Wnt8a", "Fgf5", "Otx2", "Cldn6")
episc.up <- c("Nanog", "Pou5f1", "Sox2", "L1td1", "Utf1")
ne.up <- c("Sox1", "Sox3", "Olig3", "Gbx2", "Pou3f1", "Msx3", "Slc7a3", "Zic1", "Zic2", "Nkx1-2", "Epha2", "Efnb1", "Nkx6-2")
me.down <- c("Pdgfra", "Foxc1", "Foxc2", "Isl1", "Kdr", "Mixl1", "Hand1", "Myh6")

gois <- c(me.down, ne.up, episc.up, pp.up, eo.down, eo.down.repr, bra.down)

other.gois <- c("T", "Wnt3", "Wnt3a", "Nodal", "Fgf8", "Mixl1")

#####################################################################
# List Genes
#####################################################################

down.genes <- c("Foxc2", "Gsc", "Mixl1", "Foxc1", "Pdgfra", "Kdr", "Tbx1", 
                "Amot", "Eya1", "Prrx1", "Lhx1", "Tnnt1", "Isl1", "Tbx20", 
                "Myh7", "Myh6", "Eya2", "Tbx18", "Snai1", "Snai2", "Hand1", 
                "Hand2", "Mesp1", "Pax2", "Mesp2", "Col2a1", "Myocd", "Pax3", 
                "Wt1", "Dll3", "Prickle1", "Msgn1", "Rspo3", "Osr1", "Twist1", 
                "Twist2", "Tbx6", "Meox1", "Gata6", "Gata4", "Foxa2", "Sox17", 
                "Lama1", "Tgfa", "Fn1", "Cer1", "Tgfb2", "Bmper", "Chrd", "Lgr5", 
                "Bmp2", "Fzd7", "Wnt3a", "Bmp7", "Cfc1", "Dkk1", "Fgf1", 
                "Tdgf1", "Wnt3", "Tgfb1", "Nodal", "Dact1", "Bmp4", "Dll1", 
                "Hey1", "Hhex", "Id1", "Dmbx1", "Dpf3", "Sall1", "Foxa1", 
                "Tcf7l2", "Hesx1", "Tcf7l1", "Hdac7", "Otx2", "Fbn2", "Otx1")

up.genes <- c("Six3", "Fabp7", "Ntrk2", "Zic1", "Mab21l2", "Nefl", "Olig3", 
              "Pou4f2", "Nkx1-2", "Sox3", "Sema3c", "Epha2", "Zic3", "Efnb1", 
              "Radil", "Syt11", "Slc7a3", "Nes", "Zic5", "Snph", "Nfasc", 
              "Pou4f1", "Elavl3", "Nrcam", "Grin1", "Gfra3", "Phox2a", "Nkx6-2",
              "Pax6", "Nsg2", "Msx3", "Ephb1", "Ncan", "Nova2", "Zic2", "Grik3",
              "Epha1", "Bcl11a", "Hoxa2", "Tubb3", "Sox1", "Neurod1", "Neurog1", 
              "Stmn2", "Atxn1", "Cntn2", "Neurl1a", "Sema3c", "Gap43", "Fgf5", 
              "Tbx3", "Cldn6", "Pou3f1", "Lefty2", "Gbx2", "Nanog", "L1td1", 
              "Wnt8a", "Pou5f1", "Lefty1", "Utf1", "Esrrb", "Sox2", "Lin28a", 
              "Dnmt3a", "Cbx7", "Kdm2b", "Atf7ip2")

meta2gois <- function(meta.df, cnts.avg.nrm, cnts.nrm, res, gois, outdirpath=".",
                      dataname = dataname, title.leader,
                      alpha = alpha, lFC= lFC) {
  cnts.gois <- cnts.avg.nrm[gois, ]
  
  # gene-wise normalization
  scaledata <- t(scale(t(cnts.gois)))
  scaledata <- scaledata[complete.cases(scaledata), ]
  
  # prepare heatmap
  colfunc <- colorRampPalette(c("black", "red"))
  
  # outpath
  outname <- file.path(outdirpath, paste0(title.leader, dataname, "_", time.now(), "_UO_", alpha, "_", lFC))
  
  # unordered clusters
  {
    setEPS()
    postscript(paste0(outname, ".eps"))
    # svg(paste0(outname, ".svg"))
    pheatmap::pheatmap(scaledata,
             cluster_rows = F,
             cluster_cols = F,
             cellwidth = 40,
             col = colfunc(100),
             fontsize_row = 1,
             border_color = NA
    )
    dev.off()
  }
  
  # print counts
  gois.present <- gois[gois %in% rownames(cnts.nrm)]
  res.all.gois <- gois[gois %in% rownames(res$all)]
  res.up.gois <- gois[gois %in% rownames(res$up)]
  res.down.gois <- gois[gois %in% rownames(res$down)]
  
  res.gois <- list(cnts.nrm[gois.present, ], cnts.gois, res.all.gois, res.up.gois, res.down.gois)
  names(res.gois) <- c("goi.counts", "goi.counts.avg", "goi.DE.all", "goi.DE.up", "goi.DE.down")
  xlsx2dfs::dfs2xlsx(res.gois, paste0(outname, ".xlsx"))
}

#######################################################################
# MDS plot, iMDS plot for svg using plotly
# - metapath
# - indirpath
# - lFC, alpha
# - method
# - outdirpath
#######################################################################

meta2plyMDSplot <- function(meta.df, DESeq2.obj.disp, outdirpath=".",
                            dataname = dataname, top=500,
                            launch = TRUE,
                            DE = FALSE) {
  
  title <- paste0("PCA Plot ", dataname, if (DE) { "_DE" } else {""},  " ", time.now())
  filename <- paste0("PCA_", dataname, "_", core.name(meta.df), "_", time.now(), "_DESeq2", if (DE) { "_DE" } else {""}, collapse = "_")
  
  a1 <- Glimma::glMDSPlot(if (DE) {as.data.frame(DESeq2.obj)} else {DESeq2.obj},
                  top = top,
                  path = outdirpath,
                  main = title,
                  html = filename,
                  launch = launch)
  
  df <- as.data.frame(a1$points)
  
  p <- plotly::plot_ly(df, x = df$V1, y = df$V2, text = rownames(df),
               mode = "markers", color = meta.df$condition, marker = list(size=11))
  p <- plotly::layout(p, title = title,
              xaxis = list(title = "PC1"),
              yaxis = list(title = "PC2"))
  
  export_plotly2SVG(p, 
                    filename = paste0(filename, "PC1_PC2.svg"), 
                    parent_path = outdirpath)
  
  p1 <- plotly::plot_ly(df, x = df$V2, y = df$V3, text = rownames(df),
                mode = "markers", color = meta.df$condition, marker = list(size=11))
  p1 <- plotly::layout(p1, title = title,
               xaxis = list(title = "PC2"),
               yaxis = list(title = "PC3"))
  
  export_plotly2SVG(p1, 
                    filename = paste0(filename, "PC2_PC3.svg"),
                    parent_path = outdirpath)
}


#######################################################################
# glimma plots as svg plots (2020.09.07)
#######################################################################

read_json_to_df <- function(fpath) {
  line <- readLines(fpath)[[2]]
  data.frame(jsonlite::fromJSON(substr(line, 59, (nchar(line) - 3))))
}

plot_ivolcano_to_svg <- function(fpaths, goi=NULL, ...) {
    plot_ivolcano <- function(fpath) {
        js <- read_json_to_df(fpath)
        p <- ggplot2::ggplot(js, ggplot2::aes(x=log2FC, y=log10padj, color=cols)) + 
           ggplot2::geom_point(size=2, alpha=0.3, show.legend=FALSE) + 
           ggplot2::scale_color_manual(values=c("grey", "red")) +
           ggplot2::theme_classic()
        if (!is.null(goi)) {
            p <- p + ggplot2::geom_text(label=ifelse(symbol %in% goi, symbol, ""), color="black", nudge_x = 0.4, size=3) + 
                     ggplot2::geom_point(data= js[symbol %in% goi, ], size=1, shape=20, color="black", show.legend=FALSE)
        }
        out_fpath <- file.path(dirname(dirname(fpath)), gsub("\\.js$", ".svg", basename(fpath))) 
        ggplot2::ggsave(plot = p, filename = out_fpath, ...)
        # p
    }
    sapply(fpaths, plot_ivolcano)
}

plot_imds_to_svg <- function(fpath) {
    js <- read_json_to_df(fpath)
    p <- ggplot2::ggplot(js, ggplot2::aes(x=dim1, y=dim2, colour=groups, alpha=0.5)) + 
         ggplot2::geom_point(size=5) + 
         ggplot2::theme_classic()
    out_fpath <- file.path(dirname(dirname(fpath)), gsub("\\.js$", ".svg", basename(fpath))) 
    ggplot2::ggsave(plot = p, filename = out_fpath)
    # p
}

plot_imd_to_svg <- function(fpaths, goi=NULL, ...) {
    plot_imd <- function(fpath) {
        js <- read_json_to_df(fpath)
        out_fpath <- file.path(dirname(dirname(fpath)), gsub("\\.js$", ".svg", basename(fpath))) 
        significance <- ifelse(js$Adj.PValue>=0.5, "grey", "red")
        p <- ggplot2::ggplot(js, ggplot2::aes(x=logMean, y=logFC, colour=significance)) + 
             ggplot2::geom_point(size=2, alpha=0.3, show.legend=FALSE) +       
             ggplot2::scale_color_manual(values=c("grey", "red")) +  
             ggplot2::theme_classic()
        if (!is.null(goi)) {
            p <- p + ggplot2:geom_text(label=ifelse(symbol %in% goi, symbol, ""), color="black", nudge_x = 0.4, size=3)
                 ggplot2::geom_point(data= js[symbol %in% goi, ], size=1, shape=20, color="black", show.legend=FALSE)
        }
        ggplot2::ggsave(plot = p, filename = out_fpath, ...)
        # p   
    }
    sapply(fpaths, plot_imd)
}



###########################
# GO analysis (stand 20190809)
###########################


#######################################
# load packages
#######################################

# keytypes(org.Mm.eg.db)
# 
# # [1] "ACCNUM"       "ALIAS"        "ENSEMBL"     
# # [4] "ENSEMBLPROT"  "ENSEMBLTRANS" "ENTREZID"    
# # [7] "ENZYME"       "EVIDENCE"     "EVIDENCEALL" 
# # [10] "GENENAME"     "GO"           "GOALL"       
# # [13] "IPI"          "MGI"          "ONTOLOGY"    
# # [16] "ONTOLOGYALL"  "PATH"         "PFAM"        
# # [19] "PMID"         "PROSITE"      "REFSEQ"      
# # [22] "SYMBOL"       "UNIGENE"      "UNIPROT"  
# 



#######################################
# ALIAS to ENTREZID data frame rownames
#######################################

alias2entrezid.genenames <- function(genenames) {
  "Convert using biomart ALIAS the equivalents of ENTREZID."
  "It is not a 1:1 mapping! So number of output != number of genenames!"
  genes2eg <- clusterProfiler::bitr(genenames,
                   fromType = "ALIAS",
                   toType   = "ENTREZID",
                   OrgDb    = "org.Mm.eg.db")
  genes2eg.list <- split(genes2eg, genes2eg$ALIAS)
  genes2eg.list <- lapply(genes2eg.list, function(df) df$ENTREZID)
  # keys are ALIAS and values are vector of entrezids
  entrezids.list <- genes2eg.list[genenames]
  entrezids <- unlist(entrezids.list) # flattening
  not.found.genes <- sort(setdiff(genenames, names(genes2eg.list)))
  attr(entrezids, "not.found.genes") <- not.found.genes
  found.genes     <- names(genes2eg.list)
  attr(entrezids, "found.genes") <- found.genes
  entrezids
}

alias2entrezid.df <- function(DE.df) {
  # converts rownames from ALIAS to ENTREZID
  genes <- rownames(DE.df)
  genes2eg <- clusterProfiler::bitr(genes,
                   fromType = "ALIAS",
                   toType   = "ENTREZID",
                   OrgDb    = "org.Mm.eg.db")
  
  # remove duplicated rows of aliases
  genes2eg.unique <- genes2eg[!duplicated(genes2eg$ALIAS), ]
  
  # remove duplicated rows of entrezids
  genes2eg.unique <- genes2eg[!duplicated(genes2eg$ENTREZID), ]
  df.new <- DE.df[genes2eg.unique$ALIAS, ]
  
  # again remove rows with duplicated rownames
  df.new <- df.new[!duplicated(rownames(df.new)), ]
  rownames(df.new) <- genes2eg.unique$ENTREZID
  
  # save/hide not found genes in the object # sorted looks nicer!
  not.found.genes <- sort(setdiff(genes, genes2eg.unique$ALIAS))
  attr(df.new, "not.found.genes") <- not.found.genes
  
  # save found genes in the object
  found.genes <- genes2eg.unique$ALIAS # don't sort to match with df.new!
  attr(df.new, "found.genes") <- found.genes
  
  df.new
}




#######################################
# found and not founds
#######################################

# this works also with our list of genenames

ListNotFounds <- function(eg.x.list) {
  lapply(eg.x.list, function(x) data.frame(alias = attr(x, "not.found.genes")))
}
dfListNotFounds <- ListNotFounds

dfListFounds    <- function(eg.df.list) {
  names.list <- lapply(eg.df.list, function(df) attr(df, "found.genes"))
  Map(f = function(df, names.vec) {df$alias <- names.vec; df}, eg.df.list, names.list)
}
ListFounds <- function(eg.x.list) {
  lapply(eg.x.list, function(x) data.frame(alias = attr(x, "found.genes")))
}




#######################################
# GO overrepresentation
# simplify is discussed here:
#######################################
# https://github.com/GuangchuangYu/clusterProfiler/issues/28


EnrichGO <- function(x, ont, simplifyp=FALSE) {
  # Return enchrichment for a dataframe with entrezid rownames or a vector of entrezids
  res <- clusterProfiler::enrichGO(gene          = if (is.data.frame(x)) {
                                        rownames(x)
                                  } else if (is.vector(x)) {
                                   x
                                  },
                  OrgDb         = org.Mm.eg.db::org.Mm.eg.db,
                  keyType       = "ENTREZID",
                  ont           = ont,
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.05,
                  qvalueCutoff  = 0.05,
                  readable      = TRUE)
  if (simplifyp) {
      res <- clusterProfiler::simplify(x=res,
                      cutoff = 0.7,
                      by = "p.adjust",
                      select_fun = min)
  }
  res <- clusterProfiler::setReadable(res, org.Mm.eg.db::org.Mm.eg.db, keyType = "ENTREZID")
  res
}
dfEnrichGO <- EnrichGO # x is a df

ListEnrichGO <- function(eg.x.list, ont) {
  lapply(eg.x.list, function(x) dfEnrichGO(x, ont))
}
dfListEnrichGO <- ListEnrichGO # eg.df.list


#######################################
# GO GSEA
#######################################
# prepare own geneList
# https://github.com/GuangchuangYu/DOSE/wiki/how-to-prepare-your-own-geneList

#############################################
# GO GSEA
#############################################
# for GO GSEA, the log2FoldChange for fanking is needed.
# so for a pure gene list this is not suitable

dfGseGO <- function(eg.df, ont) {
  genes <- eg.df$log2FoldChange
  names(genes) <- rownames(eg.df)
  genes <- sort(genes, decreasing = TRUE) # ensure decreasing order
  res <- clusterProfiler::gseGO(gene          = genes,
               OrgDb         = org.Mm.eg.db::org.Mm.eg.db,
               keyType       = "ENTREZID",
               ont           = ont,
               nPerm         = 1000,
               minGSSize     = 10,
               maxGSSize     = 1000,
               pAdjustMethod = "BH",
               pvalueCutoff  = 0.05)
  if (!is.null(res)) {
    res <- clusterProfiler::setReadable(res, org.Mm.eg.db, keyType = "ENTREZID")
  }
  res
}

dfListGseGO <- function(eg.df.list, ont) {
  lapply(eg.df.list, function(df) dfGseGO(df, ont))
}








#######################################################
# GO BP,CC,MF entire analysis with or without selection
# CNTS removed, since not used in GO
#######################################################

do_GO_enrichment <- function(fpathDExlsx,  # fpathCNTSxlsx, 
                             outBase, 
                             ont,
                             fpathCLUSTERxlsx = "",
                             select_clusters = c(),
                             DE.list.mode.p=T, # take only c(1, 4, 7) of xlsx sheets
                             gsea.p=T) { # do gsea? log2FC needed!
  
  fname <- basename(fpathDExlsx)
  
  if (length(select_clusters > 0)) {
    clust.numbers <- paste0("_", paste(unique(select_clusters), collapse="-"), "_")
  } else {
    clust.numbers <- "_"
  }
  
  decorateFname <- function(textPiece) gsub(".xlsx", 
                                            paste0(clust.numbers, textPiece, ont, ".xlsx"), 
                                            fname)
  overrepFname_ont <- decorateFname("over_GO_")
  gseaFname_ont <- decorateFname("gsea_GO_")
  notFoundFnameDE <- decorateFname("not_found_")
  FoundFnameDE <- decorateFname("found_")
  
  if (!dir.exists(outBase)) {
    dir.create(outBase, recursive = TRUE)
  }
  
  # read-in data
  if (DE.list.mode.p) {
    xlsxDE_dfList   <- xlsx2dfs::xlsx2dfs(fpathDExlsx)[c(1, 4, 7)]
  } else {
    xlsxDE_dfList   <- xlsx2dfs::xlsx2dfs(fpathDExlsx)
  }
  
  if (fpathCLUSTERxlsx != "" && length(select_clusters) > 0) {
    # read-in cluster tables
    cluster.dfs <- xlsx2dfs::xlsx2dfs(fpathCLUSTERxlsx)
    
    # extract gene names from clusters
    clusternames <- lapply(cluster.dfs, rownames)
    
    # extract genes from the clusters
    select_cluster_names <- unique(unlist(clusternames[select_clusters]))
    xlsxDE_dfList <- lapply(xlsxDE_dfList, function(df) df[intersect(rownames(df), select_cluster_names), ])
    
  }
  
  # translate to eg
  xlsxDE_dfListEg <- lapply(xlsxDE_dfList, alias2entrezid.df)
  
  # not founds
  not_founds_list_DE       <- dfListNotFounds(xlsxDE_dfListEg)
  
  # founds
  founds_list_DE           <- dfListFounds(xlsxDE_dfListEg)
  
  # do overrep ont
  xlsxDE_dfListEg_GO_over_ont <- dfListEnrichGO(xlsxDE_dfListEg, ont) 
  
  # do gsea ont
  if (gsea.p) {
    xlsxDE_dfListEg_GO_gsea_ont <- dfListGseGO(xlsxDE_dfListEg, ont) 
  }
  
  # write results
  xlsx2dfs::dfs2xlsx(not_founds_list_DE, file.path(outBase, notFoundFnameDE))
  xlsx2dfs::dfs2xlsx(founds_list_DE, file.path(outBase, FoundFnameDE))
  xlsx2dfs::dfs2xlsx(xlsxDE_dfListEg_GO_over_ont, file.path(outBase, overrepFname_ont))
  if (gsea.p) {
  xlsx2dfs::dfs2xlsx(xlsxDE_dfListEg_GO_gsea_ont, file.path(outBase, gseaFname_ont))
  }
}

do_all_GO_enrichment <- function(fpathDExlsx, 
                                 outBase,
                                 fpathCLUSTERxlsx = "",
                                 select_clusters = c(),
                                 DE.list.mode.p=T,
                                 gsea.p=T)  {
  do_GO_enrichment(fpathDExlsx, outBase, "BP", fpathCLUSTERxlsx, select_clusters, DE.list.mode.p, gsea.p)
  do_GO_enrichment(fpathDExlsx, outBase, "CC", fpathCLUSTERxlsx, select_clusters, DE.list.mode.p, gsea.p)
  do_GO_enrichment(fpathDExlsx, outBase, "MF", fpathCLUSTERxlsx, select_clusters, DE.list.mode.p, gsea.p)
}




################################################################################
# GO analysis only with symbol names
################################################################################

#######################################
# load packages
#######################################

# keyTypes(org.Mm.eg.db)
# 
# # [1] "ACCNUM"       "ALIAS"        "ENSEMBL"     
# # [4] "ENSEMBLPROT"  "ENSEMBLTRANS" "ENTREZID"    
# # [7] "ENZYME"       "EVIDENCE"     "EVIDENCEALL" 
# # [10] "GENENAME"     "GO"           "GOALL"       
# # [13] "IPI"          "MGI"          "ONTOLOGY"    
# # [16] "ONTOLOGYALL"  "PATH"         "PFAM"        
# # [19] "PMID"         "PROSITE"      "REFSEQ"      
# # [22] "SYMBOL"       "UNIGENE"      "UNIPROT"  
# 


#######################################
# ALIAS to ENTREZID data frame rownames
#######################################

alias2entrezid.genenames <- function(genenames) {
  "Convert using biomart ALIAS the equivalents of ENTREZID."
  "It is not a 1:1 mapping! So number of output != number of genenames!"
  genes2eg <- clusterProfiler::bitr(genenames,
                   fromType = "ALIAS",
                   toType   = "ENTREZID",
                   OrgDb    = "org.Mm.eg.db")
  genes2eg.list <- split(genes2eg, genes2eg$ALIAS)
  genes2eg.list <- lapply(genes2eg.list, function(df) df$ENTREZID)
  # keys are ALIAS and values are vector of entrezids
  entrezids.list <- genes2eg.list[genenames]
  entrezids <- unlist(entrezids.list) # flattening
  not.found.genes <- sort(setdiff(genenames, names(genes2eg.list)))
  attr(entrezids, "not.found.genes") <- not.found.genes
  found.genes     <- names(genes2eg.list)
  attr(entrezids, "found.genes") <- found.genes
  entrezids
}


#######################################
# found and not founds
#######################################

# this works also with our list of genenames

get.founds <- function(egs) {
  data.frame(alias=attr(egs, "found.genes"))
}

get.not.founds <- function(egs) {
  data.frame(alias=attr(egs, "not.found.genes"))
}

#######################################
# GO overrepresentation
# simplify is discussed here:
#######################################
# https://github.com/GuangchuangYu/clusterProfiler/issues/28


EnrichGO.genenames <- function(genenames, ont, simplifyp=FALSE) {
  # Return enchrichment for a dataframe with entrezid rownames or a vector of entrezids
  res <- clusterProfiler::enrichGO(gene          = genenames,
                  OrgDb         = org.Mm.eg.db::org.Mm.eg.db,
                  keyType       = "ALIAS",
                  ont           = ont,
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.05,
                  qvalueCutoff  = 0.05,
                  readable      = TRUE)
  if (simplifyp) {
      res <- clusterProfiler::simplify(x=res,
                      cutoff = 0.7,
                      by = "p.adjust",
                      select_fun = min)
  }
  res
}
# dfEnrichGO <- EnrichGO # x is a df




#######################################################
# GO BP,CC,MF entire analysis with or without selection
# CNTS removed, since not used in GO
#######################################################


do_GO_enrichment.genenames <- function(genenames,
                                       outfpath,
                                       simplifyp=FALSE) { # when genenames, then no GSEA
  outBase <- dirname(outfpath)
  fileName   <- basename(outfpath)
  decorateFname <- function(textPiece, fname=fileName) gsub(".xlsx", 
                                                   paste0(textPiece, "_GO.xlsx"), 
                                                   fname)
  overrepFname_ont <- decorateFname("over_GO_")
  gseaFname_ont <- decorateFname("gsea_GO_")
  notFoundFnameDE <- decorateFname("not_found_")
  FoundFnameDE <- decorateFname("found_")
  
  if (!dir.exists(outBase)) {
    dir.create(outBase, recursive = TRUE)
  }
  
  # Go enrichment from genenames
  df.bp <- EnrichGO.genenames(genenames, "BP", simplifyp=simplifyp)
  df.cc <- EnrichGO.genenames(genenames, "CC", simplifyp=simplifyp)
  df.mf <- EnrichGO.genenames(genenames, "MF", simplifyp=simplifyp)
  
  # collect the founds
  founds.df <- get.founds(genenames)
  not.founds.df <- get.not.founds(genenames)
  xlsx2dfs::dfs2xlsx(xlsx2dfs::withNames("BP", df.bp,
                              "CC", df.cc,
                              "MF", df.mf),
           outfpath)
#                      "found.genes", founds.df,
#                      "not.fount.genes", not.founds.df

}

# test: do_GO_enrichment.genenames(up.genes, "~/test.out.xlsx")

do.GO.xlsx <- function(xlsx) {
  genes <- unique(xlsx2dfs::xlsx2dfs(xlsx)[[1]]$gene)
  do_GO_enrichment.genenames(genes,
                       gsub(".xlsx", "_GO.xlsx", xlsx),
                       simplifyp=simplifyp)
  gc()
}


do.go.analyses <- function(Kmolten.rds.fpath, parallelize.p = F, simplifyp = F) {
  
  ####################################################################
  # inferred paths
  ####################################################################
  path <- dirname(dirname((Kmolten.rds.fpath)))
  jitter.xlsx.paths <- dir(path, pattern="_jitter.xlsx", recursive=T, full.names=T)
  
  if (parallelize.p) {
    bplapply(jitter.xlsx.paths, FUN=do.GO.xlsx,
      BPPARAM=mcp)
  } else {
    for (xlsx in jitter.xlsx.paths) {
      genes <- unique(xlsx2dfs::xlsx2dfs(xlsx)[[1]]$gene)
      do_GO_enrichment.genenames(genes,
                                 gsub(".xlsx", "_GO.xlsx", xlsx),
                                 simplifyp=simplifyp)
    }
  }
}

analyze.go <- function(xlsx.fpath, parallelize.p = F, simplifyp = T) {
  
  ####################################################################
  # inferred paths
  ####################################################################
  genes <- unique(xlsx2dfs::xlsx2dfs(xlsx.fpath)[[1]]$gene)
  do_GO_enrichment.genenames(genes,
                             gsub(".xlsx", if (simplifyp) {"_GO.xlsx"} else {"_GO_extended.xlsx"}, xlsx),
                             simplifyp=simplifyp)
}



################################################################################

#########################
# KEGG
#########################



#######################################
# KEGG overrepresentation 
#######################################


dfEnrichKEGG <- function(eg.df) {
  res <- clusterProfiler::enrichKEGG(gene          = rownames(eg.df),
                    keyType       = "kegg",
                    pAdjustMethod = 'BH',
                    organism      = 'mmu',
                    pvalueCutoff  =  0.05,
                    qvalueCutoff  =  0.05)
  if (!is.null(res)) {
    res <- clusterProfiler::setReadable(res, org.Mm.eg.db::org.Mm.eg.db, keyType = "ENTREZID")
  }
  res
}

dfListEnrichKEGG <- function(eg.df.list) {
  lapply(eg.df.list, function(df) dfEnrichKEGG(df))
}

#######################################
# KEGG GSEA
#######################################


dfGseKEGG <- function(eg.df, pCutOff = 0.05, org = 'mmu') {
  genes <- eg.df$log2FoldChange
  names(genes) <- rownames(eg.df)
  genes <- sort(genes, decreasing = TRUE)
  res <- clusterProfiler::gseKEGG(gene          = genes,
                 keyType       = "kegg",
                 nPerm         = 1000,
                 minGSSize     = 10,
                 organism      = org,
                 pvalueCutoff  = pCutOff)
  res <- clusterProfiler::setReadable(res, org.Mm.eg.db::org.Mm.eg.db, keyType = "ENTREZID")
  res
}

dfListGseKEGG <- function(eg.df.list) {
  lapply(eg.df.list, function(df) dfGseKEGG(df))
}


#############################
# MKEGG
#############################


#######################################
# KEGG Module overrepresentation test
#######################################


dfEnrichMKEGG <- function(eg.df, pCutOff = 0.05, qCutOff = 0.05, org = 'mmu') {
  res <- clusterProfiler::enrichMKEGG(gene          = rownames(eg.df),
                     keyType       = "kegg",
                     pAdjustMethod = 'BH',
                     organism      = org,
                     pvalueCutoff  = pCutOff,
                     qvalueCutoff  = qCutOff)
  if (!is.null(res)) {
    res <- clusterProfiler::setReadable(res, org.Mm.eg.db::org.Mm.eg.db, keyType = "ENTREZID")
  }
  res
}

dfListEnrichMKEGG <- function(eg.df.list) {
  lapply(eg.df.list, function(df) dfEnrichMKEGG(df))
}

#######################################
# KEGG Module GSEA test
#######################################


dfGseMKEGG <- function(eg.df, pCutOff = 0.05, org = 'mmu') {
  genes <- eg.df$log2FoldChange
  names(genes) <- rownames(eg.df)
  genes <- sort(genes, decreasing = TRUE)
  res <- clusterProfiler::gseMKEGG(gene          = genes,
                  keyType       = "kegg",
                  pAdjustMethod = 'BH',
                  organism      = org,
                  pvalueCutoff  = pCutOff)
  res <- clusterProfiler::setReadable(res, org.Mm.eg.db::org.Mm.eg.db, keyType = "ENTREZID")
  res
}

dfListGseMKEGG <- function(eg.df.list, pCutOff = 0.05, org = 'mmu') {
  lapply(eg.df.list, function(df) dfGseMKEGG(df, pCutOff = pCutOff, org = org))
}

do_KEGG_enrichment <- function(fpathDExlsx,
                               outBase, 
                               fpathCLUSTERxlsx = "",
                               select_clusters = c(),
                               DE.list.mode.p=T, # take only c(1, 4, 7) of xlsx sheets
                               gsea.p=T) { #gsea? log2FC column needed!
  
  fname <- basename(fpathDExlsx)
  
  if (length(select_clusters > 0)) {
    clustNumbers <- paste0("_", paste(unique(select_clusters), collapse="-"), "_")
  } else {
    clustNumbers <- "_"
  }
  
  if (!dir.exists(outBase)) {
    dir.create(outBase, recursive = TRUE)
  }
  
  decorateFname <- function(textPiece, clustNumbers_ = clustNumbers, fname_ = fname) {
    gsub(".xlsx", paste0(clustNumbers_, textPiece, ".xlsx"), fname_)
  }
  overKEGGFname <- decorateFname("over_KEGG")
  gseaKEGGFname <- decorateFname("gsea_KEGG")
  overMKEGGFname <- decorateFname("over_MKEGG")
  gseaMKEGGFname <- decorateFname("gsea_MKEGG")
  notFoundFnameDE <- decorateFname("not_found_")
  FoundFnameDE <- decorateFname("found_")
  
  # read-in data
  if (DE.list.mode.p) {
    xlsxDE_dfList   <- xlsx2dfs::xlsx2dfs(fpathDExlsx)[c(1, 4, 7)]
  } else {
    xlsxDE_dfList   <- xlsx2dfs::xlsx2dfs(fpathDExlsx)
  }
  
  if (fpathCLUSTERxlsx != "" && length(select_clusters) > 0) {
    # read-in cluster tables
    cluster.dfs <- xlsx2dfs::xlsx2dfs(fpathCLUSTERxlsx)
    
    # extract gene names from clusters
    clusternames <- lapply(cluster.dfs, rownames)
    
    # extract genes from the clusters
    select_cluster_names <- unique(unlist(clusternames[select_clusters]))
    xlsxDE_dfList <- lapply(xlsxDE_dfList, function(df) df[intersect(rownames(df), select_cluster_names), ])
    
  }
  
  # translate to eg
  xlsxDE_dfListEg <- lapply(xlsxDE_dfList, alias2entrezid.df)
  
  # not founds
  not_founds_list_DE       <- dfListNotFounds(xlsxDE_dfListEg)
  
  # founds
  founds_list_DE           <- dfListFounds(xlsxDE_dfListEg)
  
  # do overrep KEGG
  xlsxDE_dfListEg_er_KEGG <- dfListEnrichKEGG(xlsxDE_dfListEg) 
  
  # do gsea KEGG
  if (gsea.p) {
    xlsxDE_dfListEg_gs_KEGG <- dfListGseKEGG(xlsxDE_dfListEg)
  }
  
  # do overrep MKEGG
  xlsxDE_dfListEg_er_MKEGG <- dfListEnrichMKEGG(xlsxDE_dfListEg) 
  
  # do gsea MKEGG
  if (gsea.p) {
    xlsxDE_dfListEg_gs_MKEGG <- dfListGseMKEGG(xlsxDE_dfListEg)
  }
  
  # write results
  xlsx2dfs::dfs2xlsx(not_founds_list_DE, file.path(outBase, notFoundFnameDE))
  xlsx2dfs::dfs2xlsx(founds_list_DE, file.path(outBase, FoundFnameDE))
  xlsx2dfs::dfs2xlsx(xlsxDE_dfListEg_er_KEGG, file.path(outBase, overKEGGFname))
  xlsx2dfs::dfs2xlsx(xlsxDE_dfListEg_er_MKEGG, file.path(outBase, overMKEGGFname))
  if (gsea.p) {
    xlsx2dfs::dfs2xlsx(xlsxDE_dfListEg_gs_KEGG, file.path(outBase, gseaKEGGFname))
    xlsx2dfs::dfs2xlsx(xlsxDE_dfListEg_gs_MKEGG, file.path(outBase, gseaMKEGGFname))
  }
}




################################
# gskb enrichment
################################

all_ez_grouped <- mgskb::mgskb


names(all_ez_grouped)
# [1] "MousePath_Co-expression_eg.gmt" "MousePath_GO_eg.gmt"            "MousePath_Location_eg.gmt"      "MousePath_Metabolic_eg.gmt"     "MousePath_miRNA_eg.gmt"        
# [6] "MousePath_Other_eg.gmt"         "MousePath_Pathway_eg.gmt"       "MousePath_TF_eg.gmt"  



#######################################
# MSigDB/GSKB overrepresentation test
#######################################

dfEnrichDB <- function(df, gmt, pCutOff = 0.05, pAdj = "BH") {
  if (typeof(gmt) == "character") {
    term2gene <- clusterProfiler::read.gmt(gmt)
  } else {
    term2gene <- gmt
  }
  
  res <- clusterProfiler::enricher(rownames(df), TERM2GENE=term2gene, pvalueCutoff = pCutOff, pAdjustMethod = pAdj)
  if (!is.null(res)) {
    res <- clusterProfiler::setReadable(res, org.Mm.eg.db::org.Mm.eg.db, keyType = "ENTREZID")
  }
  res
}

dfListEnrichDB <- function(df.list, gmt, pCutOff = 0.05, pAdj = "BH") {
  if (typeof(gmt) == "character") {
    term2gene <- clusterProfiler::read.gmt(gmt)
  } else {
    term2gene <- gmt
  }
  lapply(df.list, function(df) dfEnrichDB(df, gmt, pCutOff = pCutOff, pAdj = pAdj))
}

#######################################
# MSigDB/GSKB GSEA
#######################################

dfGseDB <- function(df, gmt, pCutOff = 0.05, pAdj = "BH") {
  if (typeof(gmt) == "character") {
    term2gene <- clusterProfiler::read.gmt(gmt)
  } else {
    term2gene <- gmt
  }
  genes <- df$log2FoldChange
  names(genes) <- rownames(df)
  genes <- sort(genes, decreasing = TRUE)
  res <- clusterProfiler::GSEA(genes, 
              TERM2GENE=term2gene, 
              nPerm = 1000,
              minGSSize = 10,
              maxGSSize = 1000,
              pvalueCutoff = pCutOff, 
              pAdjustMethod = pAdj)
  res <- clusterProfiler::setReadable(res, org.Mm.eg.db::org.Mm.eg.db, keyType = "ENTREZID")
  res
}

dfListGseDB <- function(df.list, gmt, pCutOff = 0.05, pAdj = "BH") {
  if (typeof(gmt) == "character") {
    term2gene <- clusterProfiler::read.gmt(gmt)
  } else {
    term2gene <- gmt
  }
  lapply(df.list, function(df) dfGseDB(df, gmt, pCutOff = pCutOff, pAdj = pAdj))
}



#######################
# gskb general analysis
#######################

do_gskb_enrichment <- function(fpathDExlsx,
                               gmt,
                               gmtname,
                               outBase, 
                               fpathCLUSTERxlsx = "",
                               select_clusters = c(),
                               DE.list.mode.p=T, # take only c(1, 4, 7) of xlsx sheets
                               gsea.p=T) { #gsea? log2FC needed!
  
  fname <- basename(fpathDExlsx)
  outBase <- file.path(outBase, gmtname)
  
  if (length(select_clusters > 0)) {
    clustNumbers <- paste0("_", paste(unique(select_clusters), collapse="-"), "_")
  } else {
    clustNumbers <- "_"
  }
  
  if (!dir.exists(outBase)) {
    dir.create(outBase, recursive = TRUE)
  }
  
  decorateFname <- function(textPiece, gmtname, clustNumbers_ = clustNumbers, fname_ = fname) {
    gsub(".xlsx", paste0(clustNumbers_, textPiece, ".xlsx"), fname_)
  }
  overDBFname <- decorateFname(paste0("over_gskb_", gmtname))
  gseaDBFname <- decorateFname(paste0("gsea_gskb_", gmtname))
  notFoundFnameDE <- decorateFname(paste0("not_found_", gmtname))
  FoundFnameDE <- decorateFname(paste0("found_", gmtname))
  
  # read-in data
  if (DE.list.mode.p) {
    xlsxDE_dfList   <- xlsx2dfs::xlsx2dfs(fpathDExlsx)[c(1, 4, 7)]
  } else {
    xlsxDE_dfList   <- xlsx2dfs::xlsx2dfs(fpathDExlsx)
  }
  
  if (fpathCLUSTERxlsx != "" && length(select_clusters) > 0) {
    # read-in cluster tables
    cluster.dfs <- xlsx2dfs::xlsx2dfs(fpathCLUSTERxlsx)
    
    # extract gene names from clusters
    clusternames <- lapply(cluster.dfs, rownames)
    
    # extract genes from the clusters
    select_cluster_names <- unique(unlist(clusternames[select_clusters]))
    xlsxDE_dfList <- lapply(xlsxDE_dfList, function(df) df[intersect(rownames(df), select_cluster_names), ])
    
  }
  
  # translate to eg
  xlsxDE_dfListEg <- lapply(xlsxDE_dfList, alias2entrezid.df)
  
  # not founds
  not_founds_list_DE       <- dfListNotFounds(xlsxDE_dfListEg)
  
  # founds
  founds_list_DE           <- dfListFounds(xlsxDE_dfListEg)
  
  # do overrep gskb
  xlsxDE_dfListEg_er_DB <- dfListEnrichDB(xlsxDE_dfListEg, gmt = gmt) 
  
  # do gsea gskb
  if (gsea.p) {
    xlsxDE_dfListEg_gs_DB <- dfListGseDB(xlsxDE_dfListEg, gmt = gmt)
  }
  
  # write results
  xlsx2dfs::dfs2xlsx(not_founds_list_DE, file.path(outBase, notFoundFnameDE))
  xlsx2dfs::dfs2xlsx(founds_list_DE, file.path(outBase, FoundFnameDE))
  xlsx2dfs::dfs2xlsx(xlsxDE_dfListEg_er_DB, file.path(outBase, overDBFname))
  if (gsea.p) {
    xlsx2dfs::dfs2xlsx(xlsxDE_dfListEg_gs_DB, file.path(outBase, gseaDBFname))
  }
}



do_all_gskb_enrichments <- function(DEfpath, outBase, all_ez_grouped, xlsx.path = "", select_group = c(), DE.list.mode.p=T, gsea.p=T) {
  do_gskb_enrichment(DEfpath, 
                     all_ez_grouped$MousePath_Pathway_eg.gmt,
                     "mpathway",
                     outBase,
                     fpathCLUSTERxlsx = xlsx.path,
                     select_clusters = select_group,
                     DE.list.mode.p = DE.list.mode.p,
                     gsea.p = gsea.p)
  
  do_gskb_enrichment(DEfpath,
                     all_ez_grouped$MousePath_Metabolic_eg.gmt,
                     "mmetabolic",
                     outBase,
                     fpathCLUSTERxlsx = xlsx.path,
                     select_clusters = select_group,
                     DE.list.mode.p = DE.list.mode.p,
                     gsea.p = gsea.p)
  
  do_gskb_enrichment(DEfpath,
                     all_ez_grouped$MousePath_TF_eg.gmt,
                     "mTF",
                     outBase,
                     fpathCLUSTERxlsx = xlsx.path,
                     select_clusters = select_group,
                     DE.list.mode.p = DE.list.mode.p,
                     gsea.p = gsea.p)
  
  do_gskb_enrichment(DEfpath,
                     all_ez_grouped$MousePath_GO_eg.gmt,
                     "mGO",
                     outBase,
                     fpathCLUSTERxlsx = xlsx.path,
                     select_clusters = select_group,
                     DE.list.mode.p = DE.list.mode.p,
                     gsea.p = gsea.p)
  
  do_gskb_enrichment(DEfpath,
                     all_ez_grouped$MousePath_Other_eg.gmt,
                     "mOther",
                     outBase,
                     fpathCLUSTERxlsx = xlsx.path,
                     select_clusters = select_group,
                     DE.list.mode.p = DE.list.mode.p,
                     gsea.p = gsea.p)
  
  do_gskb_enrichment(DEfpath,
                     all_ez_grouped$MousePath_miRNA_eg.gmt,
                     "mmiRNA",
                     outBase,
                     fpathCLUSTERxlsx = xlsx.path,
                     select_clusters = select_group,
                     DE.list.mode.p = DE.list.mode.p,
                     gsea.p = gsea.p)
  
  do_gskb_enrichment(DEfpath,
                     all_ez_grouped$`MousePath_Co-expression_eg.gmt`,
                     "mCoexpr",
                     outBase,
                     fpathCLUSTERxlsx = xlsx.path,
                     select_clusters = select_group,
                     DE.list.mode.p = DE.list.mode.p,
                     gsea.p = gsea.p)
  
  do_gskb_enrichment(DEfpath,
                     all_ez_grouped$`MousePath_Location_eg.gmt`,
                     "mlocation",
                     outBase,
                     fpathCLUSTERxlsx = xlsx.path,
                     select_clusters = select_group,
                     DE.list.mode.p = DE.list.mode.p,
                     gsea.p = gsea.p)
}


################################################################################

##########################################
# for scaling
# and adding normed averaged values with dots
##########################################


# require(xlsx2dfs)
# if (!require(rlang)) {
#  install.packages("rlang")
#   require(tidyverse)
# }
#if (!require(tidyverse)) {
#  install.packages("tidyverse")
#  require(tidyverse)
# }
# require(reshape2)

#####################################################################
# helper functions for averaging tables
#####################################################################

counts.avg.build <- function(cnts.df, cols.list, groups){
  cnts.df <- as.data.frame(cnts.df)
  ncol_old <- length(colnames(cnts.df))
  ncol_limit <- length(groups) + ncol_old
  new_col_names <- c(colnames(cnts.df), groups)
  cnts.df <- cbind(cnts.df,
                   lapply(cols.list,
                          function(idxs) rowMeans(cnts.df[, idxs, drop = FALSE])))
  colnames(cnts.df) <- new_col_names
  cnts.df[, (ncol_old + 1):ncol_limit]
}

counts.std.build <- function(cnts.df, cols.list, groups){
  cnts.df <- as.data.frame(cnts.df)
  rowSds <- function(df) apply(df, 1, sd)
  ncol_old <- length(colnames(cnts.df))
  ncol_limit <- length(groups) + ncol_old
  new_col_names <- c(colnames(cnts.df), groups)
  cnts.df <- cbind(cnts.df,
                   lapply(cols.list,
                          function(idxs) rowSds(cnts.df[, idxs, drop = FALSE])))
  colnames(cnts.df) <- new_col_names
  cnts.df[, (ncol_old + 1):ncol_limit]
}

scale.raw.counts.with.SD <- function(cnts.DE.sig.fpath, meta.fpath, out.fpath = "") {
  cnts.DE.sig <- xlsx2dfs::xlsx2dfs(cnts.DE.sig.fpath)[["all.names.srt"]]
  meta.df <- read.table(meta.fpath, sep = '\t', header = T, stringsAsFactors = F)
  meta.df$condition <- factor(meta.df$condition,
                              levels = unique(meta.df$condition))
  cond <- unique(as.character(meta.df$condition))
  cond.cols <- lapply(cond,
                      function(cd) which(as.character(meta.df$condition) == cd))
  names(cond.cols) <- cond

  cnts.avg <- counts.avg.build(cnts.DE.sig, cond.cols, cond)
  std.avg  <- counts.std.build(cnts.DE.sig, cond.cols, cond)
  
  scaledata <- t(scale(t(cnts.avg)))
  
  scaled.sds <- std.avg/apply(cnts.avg, 1, sd)
  
  upper.sds.values <- scaledata + scaled.sds

  lower.sds.values <- scaledata - scaled.sds
  
  res <- list(scaled_values = scaledata,
              scaled_StdDevs = scaled.sds,
              upper_SD_values = upper.sds.values,
              lower_SD_values = lower.sds.values)
  res <- lapply(res, as.data.frame)
  if (out.fpath != "") {
    xlsx2dfs::dfs2xlsx(res, out.fpath)
  }
  res
}


scale.raw.counts.with.SD.with.orig.values <- function(cnts.DE.sig.fpath, meta.fpath, out.fpath = "") {
  cnts.DE.sig <- xlsx2dfs::xlsx2dfs(cnts.DE.sig.fpath)[["all.names.srt"]]
  meta.df <- read.table(meta.fpath, sep = '\t', header = T, stringsAsFactors = F)
  meta.df$condition <- factor(meta.df$condition,
                              levels = unique(meta.df$condition))
  cond <- unique(as.character(meta.df$condition))
  cond.cols <- lapply(cond,
                      function(cd) which(as.character(meta.df$condition) == cd))
  names(cond.cols) <- cond

  cnts.avg <- counts.avg.build(cnts.DE.sig, cond.cols, cond)
  std.avg  <- counts.std.build(cnts.DE.sig, cond.cols, cond)
  
  scaledata <- t(scale(t(cnts.avg)))
  
  scaled.sds <- std.avg/apply(cnts.avg, 1, sd)
  upper.sds.values <- scaledata + scaled.sds
  lower.sds.values <- scaledata - scaled.sds
  
  ## scale the counts
  avgs.avgs <- apply(cnts.avg, MARGIN=1, FUN=mean)
  sds.avgs  <- apply(cnts.avg, MARGIN=1, FUN=sd)
  
  # scaledata <- (cnts.avg/sds.avgs - avgs.avgs/sds.avgs) # exact scaledata
  scaled.cnts.DE.sig <- (cnts.DE.sig/sds.avgs - avgs.avgs/sds.avgs)
  
  res <- list(scaled_values = scaledata,
              scaled_StdDevs = scaled.sds,
              upper_SD_values = upper.sds.values,
              lower_SD_values = lower.sds.values,
              scaled_counts = scaled.cnts.DE.sig)
  res <- lapply(res, as.data.frame)
  if (out.fpath != "") {
    xlsx2dfs::dfs2xlsx(res, out.fpath)
  }
  res
}





############################################
# SEM
############################################

sem <- function(x) sd(x)/sqrt(length(x))

counts.sem.build <- function(cnts.df, cols.list, groups){
  cnts.df <- as.data.frame(cnts.df)
  rowSem <- function(df) apply(df, 1, sem)
  ncol_old <- length(colnames(cnts.df))
  ncol_limit <- length(groups) + ncol_old
  new_col_names <- c(colnames(cnts.df), groups)
  cnts.df <- cbind(cnts.df,
                   lapply(cols.list,
                          function(idxs) rowSem(cnts.df[, idxs, drop = FALSE])))
  colnames(cnts.df) <- new_col_names
  cnts.df[, (ncol_old + 1):ncol_limit]
}


scale.with.SD.SEM <- function(cnts.DE.sig.fpath, meta.fpath, out.fpath = "", tabname="all.names.srt") {
  cnts.DE.sig <- xlsx2dfs::xlsx2dfs(cnts.DE.sig.fpath)[[tabname]]
  meta.df <- read.table(meta.fpath, sep = '\t', header = T, stringsAsFactors = F)
  meta.df$condition <- factor(meta.df$condition,
                              levels = unique(meta.df$condition))
  cond <- unique(as.character(meta.df$condition))
  cond.cols <- lapply(cond,
                      function(cd) which(as.character(meta.df$condition) == cd))
  names(cond.cols) <- cond

  cnts.avg <- counts.avg.build(cnts.DE.sig, cond.cols, cond)
  std.avg  <- counts.std.build(cnts.DE.sig, cond.cols, cond)
  sem.avg  <- counts.sem.build(cnts.DE.sig, cond.cols, cond)
  
  scaledata <- t(scale(t(cnts.avg)))
  
  scaled.sds <- std.avg/apply(cnts.avg, 1, sd)
  scaled.sem <- sem.avg/apply(cnts.avg, 1, sem)
  
  res <- list(scaled_values = scaledata,
              scaled_StdDevs = scaled.sds,
              scaled_Sems = scaled.sem)
  res <- lapply(res, as.data.frame)
  if (out.fpath != "") {
    xlsx2dfs::dfs2xlsx(res, out.fpath)
  }
  res
}

scale.centering.with.SD.SEM <- function(cnts.DE.sig.fpath, meta.fpath, out.fpath = "", tabname="all.names.srt") {
  cnts.DE.sig <- xlsx2dfs::xlsx2dfs(cnts.DE.sig.fpath)[[tabname]]
  meta.df <- read.table(meta.fpath, sep = '\t', header = T, stringsAsFactors = F)
  meta.df$condition <- factor(meta.df$condition,
                              levels = unique(meta.df$condition))
  cond <- unique(as.character(meta.df$condition))
  cond.cols <- lapply(cond,
                      function(cd) which(as.character(meta.df$condition) == cd))
  names(cond.cols) <- cond

  cnts.avg <- counts.avg.build(cnts.DE.sig, cond.cols, cond)
  std.avg  <- counts.std.build(cnts.DE.sig, cond.cols, cond)
  sem.avg  <- counts.sem.build(cnts.DE.sig, cond.cols, cond)
  
  cnts.means.wt.dko <- rowMeans(cnts.avg[, c(1, 2)])
  
  scaledata <- t(scale(t(cnts.avg), center=cnts.means.wt.dko))
  
  scaled.sds <- std.avg/apply(cnts.avg, 1, sd)
  scaled.sem <- sem.avg/apply(cnts.avg, 1, sem)
  
  res <- list(scaled_values = scaledata,
              scaled_StdDevs = scaled.sds,
              scaled_Sems = scaled.sem)
  res <- lapply(res, as.data.frame)
  if (out.fpath != "") {
    xlsx2dfs::dfs2xlsx(res, out.fpath)
  }
  res
}

scale.centering.with.wt.dko.SD.SEM <- function(cnts.DE.sig.fpath, meta.fpath, out.fpath = "", tabname=1) {
  cnts.DE.sig <- xlsx2dfs::xlsx2dfs(cnts.DE.sig.fpath)[[tabname]]
  meta.df <- read.table(meta.fpath, sep = '\t', header = T, stringsAsFactors = F)
  meta.df$condition <- factor(meta.df$condition,
                              levels = unique(meta.df$condition))
  cond <- unique(as.character(meta.df$condition))
  cond.cols <- lapply(cond,
                      function(cd) which(as.character(meta.df$condition) == cd))
  names(cond.cols) <- cond

  cnts.avg <- counts.avg.build(cnts.DE.sig, cond.cols, cond)
  std.avg  <- counts.std.build(cnts.DE.sig, cond.cols, cond)
  sem.avg  <- counts.sem.build(cnts.DE.sig, cond.cols, cond)

  stds <- apply(cnts.DE.sig[, 1:6], 1, sd)
  
  cnts.means.wt.dko <- rowMeans(cnts.avg[, c(1, 2)])
  
  scaledata <- t(scale(t(cnts.avg), center=cnts.means.wt.dko, scale=stds))
  
  scaled.sds <- std.avg/apply(cnts.avg, 1, sd)
  scaled.sem <- sem.avg/apply(cnts.avg, 1, sem)
  res <- list(scaled_values = scaledata,
              scaled_StdDevs = scaled.sds,
              scaled_Sems = scaled.sem)
  res <- lapply(res, as.data.frame)
  if (out.fpath != "") {
    xlsx2dfs::dfs2xlsx(res, out.fpath)
  }
  res
}


scale.raw.counts.with.SD.SEM.with.orig.values <- function(cnts.DE.sig.fpath, meta.fpath, out.fpath = "") {
  cnts.DE.sig <- xlsx2dfs::xlsx2dfs(cnts.DE.sig.fpath)[["all.names.srt"]]
  meta.df <- read.table(meta.fpath, sep = '\t', header = T, stringsAsFactors = F)
  meta.df$condition <- factor(meta.df$condition,
                              levels = unique(meta.df$condition))
  cond <- unique(as.character(meta.df$condition))
  cond.cols <- lapply(cond,
                      function(cd) which(as.character(meta.df$condition) == cd))
  names(cond.cols) <- cond

  cnts.avg <- counts.avg.build(cnts.DE.sig, cond.cols, cond)
  std.avg  <- counts.std.build(cnts.DE.sig, cond.cols, cond)
  sem.avg  <- counts.sem.build(cnts.DE.sig, cond.cols, cond)
  
  scaledata <- t(scale(t(cnts.avg)))
  
  scaled.sds <- std.avg/apply(cnts.avg, 1, sd)
  scaled.sem <- sem.avg/apply(cnts.avg, 1, sd)
  upper.sds.values <- scaledata + scaled.sds
  lower.sds.values <- scaledata - scaled.sds
  upper.sem.values <- scaledata + scaled.sem
  lower.sem.values <- scaledata - scaled.sem
  
  ## scale the counts
  avgs.avgs <- apply(cnts.avg, MARGIN=1, FUN=mean)
  sds.avgs  <- apply(cnts.avg, MARGIN=1, FUN=sd)
  
  # scaledata <- (cnts.avg/sds.avgs - avgs.avgs/sds.avgs) # exact scaledata
  scaled.cnts.DE.sig <- (cnts.DE.sig/sds.avgs - avgs.avgs/sds.avgs)
  
  res <- list(scaled_values = scaledata,
              scaled_StdDevs = scaled.sds,
              scaled_SEM = scaled.sem,
              upper_SD_values = upper.sds.values,
              lower_SD_values = lower.sds.values,
              upper_SEM_values = upper.sem.values,
              lower_SEM_values = lower.sem.values,
              scaled_counts = scaled.cnts.DE.sig)
  res <- lapply(res, as.data.frame)
  if (out.fpath != "") {
    xlsx2dfs::dfs2xlsx(res, out.fpath)
  }
  res
}

###########################################
# plotting functions (they work!)
###########################################
# require(ggbeeswarm)
plot.scaled.avg <- function(scale.svg.xlsx.fpath, 
                            genes, 
                            error.type="sd", 
                            add.scatter=FALSE, 
                            out.svg.fpath="") {
  dfs <- xlsx2dfs::xlsx2dfs(scale.svg.xlsx.fpath)
  
  meta.fpath <- dir(file.path(dirname(dirname(scale.svg.xlsx.fpath)), "meta"), pattern=".txt", full.names=T)
  meta.df <- read.table(meta.fpath, sep = '\t', header = T, stringsAsFactors = F)
  meta.df$condition <- factor(meta.df$condition,
                              levels = unique(meta.df$condition))
  cond <- unique(as.character(meta.df$condition))
  cond.cols <- lapply(cond,
                      function(cd) which(as.character(meta.df$condition) == cd))
  names(cond.cols) <- cond

  colname2cond <- as.character(meta.df$condition)
  names(colname2cond) <- meta.df$sampleName
  
  cnts.melted <- reshape::melt(dfs[["scaled_values"]], variable.name="group")
  sd.melted   <- reshape::melt(dfs[["scaled_StdDevs"]], variable.name="group")
  sem.melted  <- reshape::melt(dfs[["scaled_SEM"]], variable.name="group")
  cnts.dots.melted <- reshape::melt(dfs[["scaled_counts"]], variable.name="group")
  cnts.melted$symbol <- rownames(dfs[["scaled_values"]])
  sd.melted$symbol   <- rownames(dfs[["scaled_StdDevs"]])
  sem.melted$symbol  <- rownames(dfs[["scaled_SEM"]])
  cnts.dots.melted$symbol <- rownames(dfs[["scaled_counts"]])
  cnts.melted <- cnts.melted[, c("symbol", "group", "value")]
  sd.melted   <- sd.melted[, c("symbol", "group", "value")]
  sem.melted  <- sem.melted[, c("symbol", "group", "value")]
  cnts.dots.melted <- cnts.dots.melted[, c("symbol", "group", "value")]
  cnts.dots.melted$group <- colname2cond[cnts.dots.melted$group]
  cnts.melted.sel <- cnts.melted[cnts.melted$symbol %in% genes, ]
  sd.melted.sel   <- sd.melted[sd.melted$symbol %in% genes, ]
  sem.melted.sel  <- sem.melted[sem.melted$symbol %in% genes, ]
  cnts.dots.melted.sel <- cnts.dots.melted[cnts.dots.melted$symbol %in% genes, ]
  cnts.sd.sem.melted.sel <- cbind(cnts.melted.sel, 
                                       sd=sd.melted.sel$value, 
                                       sem=sem.melted.sel$value)
  ## err <- if (error.type == "sd") {sd} else if (error.type == "sem") {sem} # doesn't work!
  
  ## that also simply doesn't work!
  p <- ggplot2::ggplot(cnts.sd.sem.melted.sel, ggplot2::aes(x = group, 
                                          y = value, 
                                          fill=factor(symbol, 
                                                      levels=genes))) +
       ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
               panel.grid.minor = ggplot2::element_blank(),
               panel.background = ggplot2::element_blank(),
               axis.line = ggplot2::element_line(colour = "black")) + # remove background and grid
       ggplot2::geom_bar(stat="identity", 
                width=.75,
                color="Black",
                position = ggplot2::position_dodge())
  if (error.type == "sd") {
    p <- p +
       ggplot2::geom_errorbar(ggplot2::aes(ymin=value-sd, 
                         ymax=value+sd), 
                     width=.6, 
                     size=.9,
                     color="Black",
                     position=ggplot2::position_dodge(.75))
  } else if (error.type == "sem") {
    p <- p +
       ggplot2::geom_errorbar(ggplot2::aes(ymin=value-sem, 
                         ymax=value+sem), 
                     width=.6, 
                     size=.9,
                     color="Black",
                     position=ggplot2::position_dodge(.75))
  }
  if (add.scatter) {
    p <- p + ggplot2::geom_jitter(data = cnts.dots.melted.sel,
                         mapping = ggplot2::aes(x = group,
                                       y = value),
                         size = 1,
                         # width=.75,
                         position = ggplot2::position_dodge(0.75))
  }
  if (out.svg.fpath != "") {
    ggplot2::ggsave(filename=out.svg.fpath, plot=p)
  }
} ## works!


create.scaled.values.sem.dots.dir.or.fpath <- function(dir.or.fpath, genes=NULL, add = NULL, error.type = "sem") {
  if (endsWith(dir.or.fpath, ".xlsx") || endsWith(dir.or.fpath, ".txt")) {
    dir.path <- dirname(dirname(dir.or.fpath))
  } else {
    dir.path <- file.path(dir.or.fpath, "k2-vs-wtn")
  }
  meta.fpath <- dir(file.path(dir.path, "meta"), pattern = ".txt", full.names = TRUE)
  cnts.DE.sig.fpath <- dir(file.path(dir.path, "DE-table"), pattern = "DE-cnts-sig-", full.names = TRUE)
  out.sig.fpath  <- gsub("-cnts-", "-scaled-avg-dots-sem-", cnts.DE.sig.fpath)
  result.1 <- scale.raw.counts.with.SD.SEM.with.orig.values(cnts.DE.sig.fpath, meta.fpath, out.sig.fpath)
  result.2 <- scale.raw.counts.with.SD.SEM.with.orig.values(cnts.DE.sig.fpath, meta.fpath, file.path(out.revision.dir, basename(out.sig.fpath)))

  if (!is.null(genes) && !is.null(add)) {
    out.svg.fpath  <- gsub(".txt", ".svg", gsub(".xlsx", ".svg", out.sig.fpath))
    out.svg.fpath  <- gsub("-sem-", paste0("-sem-", add, "-"), out.svg.fpath)
    out.svg.fpath  <- gsub("-sem-",
                           paste0("-", error.type, "-"),
                           out.svg.fpath)
    plot.scaled.avg(out.sig.fpath,
                    genes = genes,
                    error.type=error.type,
                    add.scatter=TRUE,
                    out.svg.fpath=out.svg.fpath)
    plot.scaled.avg(out.sig.fpath,
                    genes = genes,
                    error.type=error.type,
                    add.scatter=TRUE,
                    out.svg.fpath=file.path(out.revision.dir, basename(out.svg.fpath)))
  }
}


robustscale.counts.with.MAD <- function(cnts.DE.sig.fpath, meta.fpath, out.fpath = "", tabname="all.names.srt") {
  cnts.DE.sig <- xlsx2dfs::xlsx2dfs(cnts.DE.sig.fpath)[[tabname]]
  meta.df <- read.table(meta.fpath, sep = '\t', header = T, stringsAsFactors = F)
  meta.df$condition <- factor(meta.df$condition,
                              levels = unique(meta.df$condition))
  cond <- unique(as.character(meta.df$condition))
  cond.cols <- lapply(cond,
                      function(cd) which(as.character(meta.df$condition) == cd))
  names(cond.cols) <- cond

  cnts.avg <- counts.avg.build(cnts.DE.sig, cond.cols, cond)
  std.avg  <- counts.std.build(cnts.DE.sig, cond.cols, cond)
  
  robustscaled   <- quantable::robustscale(cnts.avg, dim=1, center=TRUE, scale=T, preserveScale=F)
  robustscaledata <- robustscaled$data
  
  robustscaled.mads <- robustscaled$mads
  
  upper.mads.values <- robustscaledata + robustscaled.mads

  lower.mads.values <- robustscaledata - robustscaled.mads
  
  res <- list(robustscaled_values = robustscaledata,
              robustscaled_MADs = robustscaled.mads, # mean absolute deviation
              upper_MAD_values = upper.mads.values,
              lower_MAD_values = lower.mads.values)
  res <- lapply(res, as.data.frame)
  if (out.fpath != "") {
    xlsx2dfs::dfs2xlsx(res, out.fpath)
  }
  res
}

## doesnt' work rownames problem
robustscale.counts.with.SD <- function(cnts.DE.sig.fpath, meta.fpath, out.fpath = "", tabname="all.names.srt") {
  cnts.DE.sig <- xlsx2dfs::xlsx2dfs(cnts.DE.sig.fpath)[[tabname]]
  meta.df <- read.table(meta.fpath, sep = '\t', header = T, stringsAsFactors = F)
  meta.df$condition <- factor(meta.df$condition,
                              levels = unique(meta.df$condition))
  cond <- unique(as.character(meta.df$condition))
  cond.cols <- lapply(cond,
                      function(cd) which(as.character(meta.df$condition) == cd))
  names(cond.cols) <- cond

  cnts.avg <- counts.avg.build(cnts.DE.sig, cond.cols, cond)
  std.avg  <- counts.std.build(cnts.DE.sig, cond.cols, cond)
  
  robustscaled   <- quantable::robustscale(cnts.avg, dim=1, center=TRUE, scale=F, preserveScale=F)
  robustscaledata <- robustscaled$data
  
  robustscaled.mads <- robustscaled$mads
  
  upper.mads.values <- robustscaledata + robustscaled.mads

  lower.mads.values <- robustscaledata - robustscaled.mads
  
  res <- list(robustscaled_values = robustscaledata,
              robustscaled_MADs = robustscaled.mads, # mean absolute deviation
              upper_MAD_values = upper.mads.values,
              lower_MAD_values = lower.mads.values)
  res <- lapply(res, as.data.frame)
  if (out.fpath != "") {
    xlsx2dfs::dfs2xlsx(res, out.fpath)
  }
  res
}


#######################################
# rpkm
#######################################

# require(xlsx2dfs)
# require(scater)

load("data/sym2len.rda")
print.rpkm <- Vectorize(function(cnts.fpath) {
  raw.count.sheet <- "raw-counts"
  out.fpath <- ""
  cnts <- xlsx2dfs::xlsx2dfs(cnts.fpath)[[raw.count.sheet]]
  sce <- SingleCellExperiment::SingleCellExperiment(assays=list(counts=as.matrix(cnts)))
  sce <- scater::normalize(sce)
  
  # load("../data/sym2len.rda")
  cnts.rpkm.len <- scater::calculateFPKM(object = sce, 
                                         effective_length = sym2len[rownames(cnts)], 
                                         use_size_factors=FALSE)
    
  # sym2len is a clojure!
  cnts.rpkm.len.cl <- na.omit(cnts.rpkm.len)
  if (out.fpath == "") {
    out.fpath <-gsub("raw-counts", "rpkm-raw", cnts.fpath)
  }
  xlsx2dfs::dfs2xlsx(xlsx2dfs::withNames("rpkm_raw_counts", cnts.rpkm.len.cl),
           out.fpath)
  cnts.rpkm.len.cl
})


# Error in file(description = xlsxFile) : invalid 'description' argument
# this is when several xlsxFiles are found for cnts.fpath

###################################################
# db handling (for automatic meta file generation)
###################################################

load_db <- function(db_fpath = NULL) {
    if (is.null(db_fpath)) {
        fpath <- "../inst/extdata/db.json"
    }
    jsonlite::fromJSON(fpath)
}
# it is a named list with vectors of fpaths, 
# however, when empty they are empty lists.

#' Add key/id file_path pair to db.json.
#'
#' Add key/id file_path pair to db.json.
#'
#' @param db_fpath file path to database - "../inst/extdata/db.json"
#' @param id unique key/id for group of samples
#' @param fpaths vector of file paths (biological replicates)
#' @param sort sort the keys/ids alphabetically?
#' @return nothing
#' @export
add_to_db <- function(db_fpath, id, fpaths, sort=FALSE) {
  db <- load_db()
  l <- list(c(fpaths))
  names(l) <- id
  db <- append(db, l)
  if (sort) {
    db <- db[sort(names(db))]
  }
  writeLines(jsonlite::toJSON(db, pretty=TRUE), db_fpath)
}

#' Remove entry (key/id) from db.json.
#'
#' Remove entry (key/id) from db.json.
#'
#' @param id unique key/id for group of samples
#' @param db_fpath file path to database - default is "../inst/extdata/db.json"
#' @return nothing
#' @export
remove_from_db <- function(id, db_fpath = NULL) {
  if (is.null(db_fpath)) db_fpath <- "../inst/extdata/db.json"
  db <- load_db()
  db <- db[names(db) != id]
  writeLines(jsonlite::toJSON(db, pretty=TRUE), db_fpath)
}

# setwd("~/local/src-r/rnaseqeasy/R")
# db <- load_db()
# add_to_db("../inst/extdata/db.json",
#           "test-id", c("a/b/c", "d/e/f"))
# remove_from_db("test-id")
# 

#######################################################
# automated meta file generation
#######################################################

meta_lines <- function(key, db, name=NULL, testing="") {
  # Returns meta lines df.
  # sampleName fileName	fileName	condition	testing
  # test = ["", "denom", "num"]
  if (is.null(name)) {
    name <- key
  }
  fileName <- db[[key]]
  k <- length(fileName)
  sampleName <- sapply(1:k, function(i) paste0(name, "-", i))
  condition <- rep(name, k)
  testing <- rep(testing, k)
  data.frame(sampleName = sampleName,
             fileName = fileName,
             condition = condition,
             testing = testing)
}

#' Generate entire meta file out of information.
#'
#' Automatic creation of meta file.
#'
#' @param keys list of keys for db. Order determines testing!: First element -> "denom", Second element: "num", rest is ""
#' @param db contains key and list of paths (vector of paths) to key files
#' @param meta_dir the directory to which meta file should be written to. Name is generated automatically using project_name and source
#' @param project_name inserted as signifier into meta filename
#' @param source inserted as signifier into meta filename
#' @param sep separator in meta file - tab
#' @param names_ names to be used instead of keys
#' @return list("fpath" = fpath, "core_name" = core_name)
#' @export
create_meta <- function(keys, db, meta_dir, project_name, source, sep="\t", names_=NULL) {
        if (is.null(names_)) {
                names_ <- keys
        }
        # automatic file_name/path generation:
        core_name <- paste0(names_[2], "-vs-", names_[1])
        if (length(names_) > 2) core_name <- paste0(core_name, paste(names_[-c(1, 2)], collapse=""))
        fname <- paste0("meta_",
                        project_name, "_",
                        source, "_",
                        core_name,
                        ".txt")
    fpath <- file.path(meta_dir, fname)

    # automatic file content generation:
    testings <- c("denom", "num", rep("", length(keys) - 2))
    res <- lapply(1:length(names_), function(i) meta_lines(keys[i], db, names_[i], testings[i]))
    df <- Reduce(rbind, res)
    names(df) <- c("sampleName", "fileName", "condition", "testing")
    print(df)
    write.table(x = df, file = fpath, col.names = TRUE, row.names = FALSE, quote=FALSE, sep=sep)
    list("fpath" = fpath, "core_name" = core_name)
}




#' Perform automated analysis.
#'
#' Perform automated analysis.
#'
#' @param keys list of keys for db. Order determines testing!: First element -> "denom", Second element: "num", rest is ""
#' @param db contains key and list of paths (vector of paths) to key files
#' @param meta_dir the directory to which meta file should be written to. Name is generated automatically using project_name and source
#' @param project_name inserted as signifier into meta filename
#' @param source inserted as signifier into meta filename
#' @param sep separator in meta file - tab
#' @param names names to be used instead of keys
#' @param out_dir output directory for the analysis
#' @param alpha p-value for DESeq2 analysis
#' @param lFC log2FoldChange threshold for significance
#' @param ks number of clusters
#' @param go "skip" or don't skip "yes"
#' @param goi gene of interest for iVolcano and iMDS plot svg
#' @return nothing
#' @export
analyze <- function(keys, 
                    db, 
                    meta_dir, 
                    project_name, 
                    source, 
                    sep, 
                    names, 
                    out_dir,
                    alpha = 0.05,
                    lFC = 1,
                    ks = "12",
                    go = "skip",
                    goi = NULL) {
    tmp <- create_meta(keys = keys,
                     db = db,
                     meta_dir = meta_dir,
                     project_name = project_name,
                     source = source,
                     sep = sep,
                     names = names)
    meta_fpath <- tmp[["fpath"]]
    dirname    <- tmp[["core_name"]]
    data_name <- paste0(project_name, "_", source)
    ks <- as.numeric(strsplit(ks, "\\s+")[[1]])
    indirpath <- ""
    dataname   <- data_name
    outdirpath <- out_dir
    dir.create(outdirpath, recursive=TRUE, showWarnings=FALSE)
    metapath   <- meta_fpath
    meta.df <- read.table(metapath, sep = '\t',
                        header = TRUE,
                        stringsAsFactors = FALSE)
    meta.df$condition <- factor(meta.df$condition, # ensures preferred order
                              levels = unique(meta.df$condition))
    denom <- function(meta.df) as.character(meta.df$condition[meta.df$testing == "denom"][1])
    num   <- function(meta.df) as.character(meta.df$condition[meta.df$testing == "num"][1])
    core.name <- function(meta.df) paste0(num(meta.df), "-vs-", denom(meta.df), collapse = '')

    outdirpath <- file.path(outdirpath, project_name, core.name(meta.df))
    meta.outdir <- file.path(outdirpath, "meta")
    cnts.outdir <- file.path(outdirpath, "count-table")
    DE.outdir <- file.path(outdirpath, "DE-table")
    plot.outdir <- file.path(outdirpath, "glimma-plots")
    hm.outdir <- file.path(outdirpath, "heatmap")
    hm.all.outdir <- file.path(hm.outdir, "all")
    GO.outdir <- file.path(outdirpath, "GO")
    pathway.outdir <- file.path(outdirpath, "pathway")
    gskb.outdir <- file.path(outdirpath, "gskb")
    PCA.outdirpath <- file.path(outdirpath, "pca")

    dir.create(path=outdirpath, recursive = TRUE, showWarnings = FALSE)
    dir.create(path=meta.outdir, recursive = TRUE, showWarnings = FALSE)
    dir.create(path=cnts.outdir, recursive = TRUE, showWarnings = FALSE)
    dir.create(path=DE.outdir, recursive = TRUE, showWarnings = FALSE)
    dir.create(path=plot.outdir, recursive = TRUE, showWarnings = FALSE)
    dir.create(path=hm.all.outdir, recursive = TRUE, showWarnings = FALSE)
    dir.create(path=GO.outdir, recursive = TRUE, showWarnings = FALSE)
    dir.create(path=pathway.outdir, recursive = TRUE, showWarnings = FALSE)
    dir.create(path=gskb.outdir, recursive = TRUE, showWarnings = FALSE)
    dir.create(path=PCA.outdirpath, recursive = TRUE, showWarnings = FALSE)

    # document meta table
    write.table(meta.df,
              file = file.path(meta.outdir, basename(metapath)),
              sep = "\t", row.names = FALSE, quote = FALSE)

    DESeq2.obj <- meta2DESeq2.obj(meta.df, indirpath)
    DESeq2.obj.disp <- meta2DESeq2.obj(meta.df, indirpath, normalized = TRUE)

    # create count tables
    cnts.raw <- meta2cnts(meta.df, DESeq2.obj, outdirpath = cnts.outdir,
                        dataname = dataname,
                        printp = TRUE, normalized = FALSE, averaged = FALSE,
                        sheetName = "raw.all")
    cnts.nrm <- meta2cnts(meta.df, DESeq2.obj.disp, outdirpath = cnts.outdir,
                        dataname = dataname,
                        printp = TRUE, normalized = TRUE, averaged = FALSE,
                        sheetName = "normalized.all")
    cnts.avg.nrm <- meta2cnts(meta.df, DESeq2.obj.disp, outdirpath = cnts.outdir,
                            dataname = dataname,
                            printp = TRUE, normalized = TRUE, averaged = TRUE,
                            sheetName = "avg.normalized.all")

    # create DE tables
    res <- DEanalysis(meta.df, DESeq2.obj.disp, outdirpath = DE.outdir,
                    dataname = dataname,
                    printp = TRUE)
    resSig <- DEanalysis(meta.df, DESeq2.obj.disp, outdirpath = DE.outdir,
                       dataname = dataname,
                       printp = TRUE,
                       filterp = TRUE, alpha = alpha, lFC = lFC)
    print.cnts.DE.sortings(cnts.nrm, 
                         resSig, 
                         file.path(paste0("DE-cnts-sig-", dataname, "_", 
                                          core.name(meta.df), "_",
                                          time.now(), "-", alpha, "-", lFC, ".xlsx")), 
                         dir = DE.outdir)

    print.cnts.DE.sortings(cnts.avg.nrm$`nrm-counts-avg`, 
                         resSig, 
                         file.path(paste0("DE-cnts-avg-sig-", dataname, "_", 
                                          core.name(meta.df), "_",
                                          time.now(), "-", alpha, "-", lFC, ".xlsx")), 
                         dir = DE.outdir)

    # scaled normalized counts
    DE.cnts.fpath <- dir(DE.outdir, pattern=paste0("DE-cnts-sig-", data_name , "_", 
                                                   core.name(meta.df)), full.names=T)[1]
    out.fpath <- gsub("-cnts-", "-scaled-avg-", DE.cnts.fpath)
    result <- scale.raw.counts.with.SD(DE.cnts.fpath, metapath, out.fpath)

    # robust scaled normalized counts
    # out.fpath.r <- gsub("-cnts-", "-robustscaled-mad-avg-", DE.cnts.fpath)
    # result.robust <- robustscale.counts.with.MAD(DE.cnts.fpath, metapath, out.fpath.r)
    #### due to scater error inactiate this
    # # rpkm counts
    # raw.fpath   <- dir(cnts.outdir, pattern="raw-counts-", full.names=T)
    # rpkm <- print.rpkm(raw.fpath)
    
    # glimma plots
    meta2iMDSplot(meta.df, DESeq2.obj.disp, outdirpath,
                  dataname, top=300, launch = FALSE)

    meta2iVolcano(meta.df, DESeq2.obj, DESeq2.obj.disp, outdirpath,
                  dataname, alpha = alpha, lFC = lFC,
                  launch = FALSE)
    
    # glima plots as svg
    ivolcano_fpath <- dir(file.path(outdirpath, "glimma-plots", "js"), 
                          pattern="iVolcano\\_", 
                          full.names=TRUE)
    imds_fpath <- dir(file.path(outdirpath, "glimma-plots", "js"), 
                          pattern="iMDS\\_", 
                          full.names=TRUE)
    imd_fpath <- dir(file.path(outdirpath, "glimma-plots", "js"), 
                          pattern="iMD\\_", 
                          full.names=TRUE)
    if (length(ivolcano_fpath) > 1) ivolcano_fpath <- ivolcano_fpath[length(ivolcano_fpath)]
    if (length(imds_fpath) > 1) imds_fpath <- imds_fpath[length(imds_fpath)]
    if (length(imd_fpath) > 1) imd_fpath <- imd_fpath[length(imd_fpath)]   
    
    try({
      plot_ivolcano_to_svg(ivolcano_fpath, goi=goi)
      plot_imds_to_svg(imds_fpath)
      plot_imd_to_svg(imd_fpath, goi=goi)
    })






    # heatmaps reproducible
    for (k in ks) {
        set.seed(123)
    try(scaledata.Kmolten.core.list <- meta2heatmap(meta.df, cnts.avg.nrm$`nrm-counts-avg`, resSig,
                                                    outdirpath = hm.all.outdir, 
                                                    selected.genes = NULL, 
                                                    name.add = "all",
                                                    dataname = dataname,
                                                    k = k, 
                                                    printp=TRUE,
                                                    alpha=alpha, 
                                                    lFC=lFC, 
                                                    filterp=TRUE,
                                                    xlsxp=TRUE, 
                                                    csvp=FALSE, 
                                                    tsvp=FALSE))    
    }



    if (go != "skip") {
    
        last <- function(l) l[length(l)]
        DEfpath <- last(dir(DE.outdir, pattern = "DE_sig", full.names = TRUE))


        try(do_all_GO_enrichment(DEfpath, GO.outdir))
        try(do_KEGG_enrichment(DEfpath, pathway.outdir))
        try(do_all_gskb_enrichments(DEfpath,
                              gskb.outdir,
                              all_ez_grouped))
    }
    
    fpath <- file.path(outdirpath, paste0("run", time.now(), ".RData"))
    save.image(file = fpath)

    sink(paste0(fpath, ".log"), split=TRUE)
        sessionInfo()
    sink()

    gc()
}








##########################################
# preparing pairings
##########################################

# create from this list the common denominator - also giving the connectors

common.string <- function(s1, s2, acc="") {
  # Return common beginning string
  first.s1 = substr(s1, 1, 1)
  first.s2 = substr(s2, 1, 1)
  if (s1 == "" && s2 == "" || first.s1 != first.s2) {
    acc
  } else  { # s1[1] == s2[1]
    common.string(substr(s1, 2, nchar(s1)), substr(s2, 2, nchar(s2)), paste0(acc, first.s1))
  }
} # works!

cdr <- function(vec) if (length(vec) < 2) c() else vec[2:length(vec)]
car <- first <- function(vec) vec[1]
cadr <- second <- function(vec) vec[2]
cons <- function(x, vec) c(x, vec)
nreverse <- rev
null <- function(x) length(x) == 0

grouper <- function(vec, acc=c(), last.common=c(), prev="") {
  if (null(vec)) {
    nreverse(acc)
  } else if (length(vec) == 1) {
    grouper(cdr(vec), cons(last.common, acc))
  } else {
    next.common <- common.string(first(vec), second(vec))
    if (null(last.common)) {
      grouper(cdr(vec), cons(next.common, acc), next.common, car(vec))
    } else {
      prev.common <- common.string(prev, car(vec))
      if (last.common == prev.common) {
        grouper(cdr(vec), cons(prev.common, acc), prev.common, car(vec))
      } else {
        grouper(cdr(vec), cons(next.common, acc), next.common, car(vec))
      }
    }
  }
} # works!

split_by_group <- function(vec, indexes=FALSE) {
  dfs <- split(data.frame(x= if (indexes) 1:length(vec) else vec,
                          stringsAsFactors=FALSE),
               grouper(vec))
  lapply(dfs, function(df) df$x)
}

average_df <- function(df, col_names=NULL) {
  # averages df by similarity of columnnames
  if (is.null(col_names)) col_names <- colnames(df)
  column_indexes <- split_by_group(col_names, indexes=TRUE)
  res.df <- Reduce(cbind, lapply(column_indexes, function(idxs) rowMeans(df[, idxs]))) # rowSds
  colnames(res.df) <- names(column_indexes)
  res.df
} # works as expected!


trimr <- function(s, end) {
  gsub(paste0(end, "+$"), "", s)
}

extract_group_name <- function(pathroot, split_ending=c("_", "-")) {
  res <- basename(pathroot)
  gsub(paste0("(", paste(split_ending, collapse="|"), ")+$"), "", res)
}

extract_group_names <- function(pathroots, split_ending=c("_", "-")) {
  sapply(pathroots, extract_group_name, split_ending)
}



#' Group list of paths by change of common substrings automatically.
#'
#' Group list of paths.
#'
#' @param paths the file paths of the symbol counts of each sample
#' @param split_ending by which substring the path string should be split and grouped.
#' @return list of grouped file paths named by file name part
#' @export
generate_groups_list <- function(paths, split_ending=c("_", "-")) {
  res <- split_by_group(paths)
  names(res) <- extract_group_names(names(res), split_ending=split_ending)
  res
} # works!!

generate_one_to_one_dict <- function(names_) {
  # generate 1:1 dict - in our case list
  res <- as.list(names_)
  names(res) <- names_
  res
}

#' Group list of paths by two layers.
#'
#' Group 2-layer list of paths.
#'
#' @param vec the vector of paths.
#' @return list of grouped file paths - 2-layer grouping
#' @export
build_two_level_list <- function(vec) {
  js <- generate_groups_list(vec)
  js_ <- generate_groups_list(names(js))
  res <- lapply(names(js_), function(key) as.vector(sapply(js_[[key]], function(x) js[[x]])))
  names(res) <- names(js_)
  res
}

# the easiest method however:
# manually split filepaths vector by
# js <- split(vec, your_manual_conditions_vector) # it also gives the correct conditions names!


####################################
# filter by unique list
####################################

filter_previous_elements <- function(l) {
  names_ <- names(l)
  pool <- c()
  res <- list()
  for (i in seq(l)) {
    el <- l[[i]]
    res[[i]] <- setdiff(el, pool)
    pool <- unique(c(pool, el))
  }
  names(res) <- names_
  res
}

#' Group list by vector containing all beginning strings.
#'
#' Group list of paths by filename beginnings. It sorts by longest forms first and 
#' removes all elements which were already chosen previously in the list.
#'
#' @param vec the vector of paths.
#' @param kinds the vector of beginning strings which becomes name of the result list
#' @return list of grouped file paths
#' @export
match_beginning <- function(vec, kinds) {
  kinds <- rev(sort(kinds))
  vec_ <- basename(vec)
  res <- lapply(kinds, function(s) {
    vec[grepl(paste0("^", s), vec_)]
  })
  names(res) <- kinds
  # filter out those which were in previous
  filter_previous_elements(res)
}










