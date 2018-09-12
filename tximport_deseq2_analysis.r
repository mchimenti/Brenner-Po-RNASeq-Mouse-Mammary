## Analysis of Po/Brenner Mammary Gland Tissue from Mouse
## Date: 8.30.2018
## Author: Michael Chimenti
## Organism: mm10 / mouse
## Aligners: hisat2 / salmon
## Design: Diet + Treatment 
## Reps: 3

##########
## Imports
##########

#source("https://bioconductor.org/biocLite.R")
#biocLite("DEGreport")

#negative binomial GLM and related
library('DESeq2')
library('calibrate')
library('tximport')
library('readr')
#annotation
library('biomaRt')
library("AnnotationDbi")
library("org.Hs.eg.db")
#Exploratory analysis
library('tidyverse')
library('pcaExplorer')
#pathway and gene clusters
library('DEGreport')
#library(pathview)
#library(gage)
#library(gageData)
#library(ggplot2)

setwd("~/iihg/RNA_seq/brenner_lab/project_po_sept2018/") 

###########
##Function Defs
###########

get_annotation <- function(dds, biomart_dataset, idtype){
  if(is.null(biomart_dataset))
    stop("Select a species to generate the corresponding annotation.
         To obtain a list, type mart = useMart('ensembl'), followed by listDatasets(mart).")
  
  mart <- useMart(biomart="ENSEMBL_MART_ENSEMBL",
                  host="www.ensembl.org",
                  dataset=biomart_dataset)
  anns <- getBM(attributes = c(idtype, "external_gene_name", "description"),
                filters = idtype,
                values = rownames(dds),
                mart = mart)
  
  # keep and match with the ones that are actually there
  anns2 <- anns[match(rownames(dds), anns[, 1]), ]
  rownames(anns2) <- rownames(dds)
  # rename the columns rsp. add row names to be consistent with other function
  colnames(anns2) <- c("gene_id","gene_name","description")
  
  return(anns2)
}

## Volcano Plot function 
volcanoplot <- function (res, lfcthresh=2, sigthresh=0.05, main="Volcano Plot", legendpos="topright", labelsig=TRUE, textcx=1, ...) {
  with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main=main, ...))
  with(subset(res, padj<sigthresh ), points(log2FoldChange, -log10(pvalue), pch=20, col="red", ...))
  with(subset(res, abs(log2FoldChange)>lfcthresh), points(log2FoldChange, -log10(pvalue), pch=20, col="orange", ...))
  with(subset(res, padj<sigthresh & abs(log2FoldChange)>lfcthresh), points(log2FoldChange, -log10(pvalue), pch=20, col="green", ...))
  if (labelsig) {
    require(calibrate)
    with(subset(res, padj<sigthresh & abs(log2FoldChange)>lfcthresh), textxy(log2FoldChange, -log10(pvalue), labs=ext_gene, cex=textcx, offset=0.3, ...))
  }
  legend(legendpos, xjust=1, yjust=1, legend=c(paste("FDR<",sigthresh,sep=""), paste("|LogFC|>",lfcthresh,sep=""), "both"), pch=20, col=c("red","orange","green"))
}

#######################################
## tximport > DESeq2 
#######################################
samples <- read.table("samples.csv", sep=',', header=TRUE)
rownames(samples) <- samples$sample
samples$group <- paste0(samples$diet, samples$treatment)

files <- file.path(getwd(), samples$sample, 'salmon', 'quant.sf')
names(files) <- samples$sample

tx2gene <- read_csv(file.path(getwd(), "tx2gene.csv"), col_names = FALSE)

#tx2gene$X1 <- tx2gene$X1 %>%
#  strsplit(split = '.', fixed = TRUE) %>%
#  sapply( "[", 1)  ## obtuse sapply statement needed b/c of annoying way strsplit returns list of lists

txi <- tximport(files, type="salmon", tx2gene=tx2gene)
ddsTxi <- DESeqDataSetFromTximport(txi,
                                   colData = samples,
                                   design = ~ group)

ddsTxi <- ddsTxi[ rowSums(counts(ddsTxi)) > 5, ]
ddsTxi <- DESeq(ddsTxi)

anno <- get_annotation(ddsTxi, 'mmusculus_gene_ensembl','ensembl_gene_id')
anno <- na.omit(anno)

rldTxi <- rlog(ddsTxi, blind=FALSE)
pcaExplorer(dds=ddsTxi,annotation=anno,rlt=rldTxi)

## look at dispersion estimates 
plotDispEsts(ddsTxi)
plotMA(object = ddsTxi, alpha = 0.05)
plotPCA(object = rldTxi, intgroup = 'group')

# drop sample 5 outlier

ddsTxi <- ddsTxi[ , ddsTxi$sample != 's5']
ddsTxi$sample <- droplevels(ddsTxi$sample)
ddsTxi <- DESeq(ddsTxi)

##DE testing 
library("IHW")

## normal chow treated vs. non-treated
res_nc_treat <- results(ddsTxi, contrast = c("group","ncy","ncn"), filterFun = ihw)
res_nc_treat <- na.omit(res_nc_treat)  #drop NA rows
res_nc_treat_sig <- res_nc_treat[res_nc_treat$padj < 0.1 & res_nc_treat$baseMean > 5.0,]
res_nc_treat_ord <- res_nc_treat_sig[order(res_nc_treat_sig$padj),]
res_nc_treat_ord$ext_gene <- anno[row.names(res_nc_treat_ord), "gene_name"]

png("test.png", 1200, 1500, pointsize=20, res=100)
volcanoplot(res_nc_treat_ord, main = "DE genes, Normal Chow Treated vs. Non", lfcthresh=1.0, sigthresh=0.05, textcx=.8, xlim=c(-4, 4), ylim = c(4,10))
dev.off()

degPlot(dds = ddsTxi, res = res_nc_treat_ord, n = 9, xs = 'diet', group = 'treatment')

#high-fat treated vs. non-treated
res_hfd_treat <- results(ddsTxi, contrast = c("group","hfdy", "hfdn"), filterFun = ihw)
res_hfd_treat <- na.omit(res_hfd_treat)  #drop NA rows
res_hfd_treat_sig <- res_hfd_treat[res_hfd_treat$padj < 0.1 & res_hfd_treat$baseMean > 5.0,]
res_hfd_treat_ord <- res_hfd_treat_sig[order(res_hfd_treat_sig$padj),]
res_hfd_treat_ord$ext_gene <- anno[row.names(res_hfd_treat_ord), "gene_name"]

degPlot(dds = ddsTxi, res = res_hfd_treat_ord, n = 9, xs = 'diet', group = 'treatment')

png("test.png", 1200, 1500, pointsize=20, res=100)
volcanoplot(res_hfd_treat_ord, main = "DE Genes, HFD treated vs. Non-treated", lfcthresh=1.0, sigthresh=0.05, textcx=.8, xlim=c(-4, 4), ylim = c(4,10))
dev.off()


## high fat diet vs normal chow, non-treated
res_hfd_nc <- results(ddsTxi, contrast = c("group","hfdn", "ncn"), filterFun = ihw)
res_hfd_nc <- na.omit(res_hfd_nc)  #drop NA rows
res_hfd_nc_sig <- res_hfd_nc[res_hfd_nc$padj < 0.1 & res_hfd_nc$baseMean > 5.0,]
res_hfd_nc_ord <- res_hfd_nc_sig[order(res_hfd_nc_sig$padj),]
res_hfd_nc_ord$ext_gene <- anno[row.names(res_hfd_nc_ord), "gene_name"]

degPlot(dds = ddsTxi, res = res_hfd_nc_ord, n = 9, xs = 'diet', group = 'treatment')

png("vol_highfat_vs_chow_NO_NR.png", 1200, 1500, pointsize=20, res=100)
volcanoplot(res_hfd_nc_ord, main = "", lfcthresh=1.0, sigthresh=0.05, textcx=.7, xlim=c(-2.5, 6), ylim = c(3,10))
dev.off()

my_cols <- c("baseMean","log2FoldChange","padj", "ext_gene")
write.csv(x = res_nc_treat_ord[,my_cols], file = "de_genes_normal_chow_treated_v_non_padj_10percent.csv")
write.csv(x = res_hfd_treat_ord[,my_cols], file = "de_genes_highfat_treated_v_non_padj_10percent.csv")
write.csv(x = res_hfd_nc_ord[,my_cols], file = "de_genes_nontreated_highfat_vs_chow_padj_10percent.csv")

## high fat diet vs normal chow, treated
res_hfd_nc_treat <- results(ddsTxi, contrast = c("group", "hfdy", "ncy"), filterFun = ihw)
res_hfd_nc_treat <- na.omit(res_hfd_nc_treat)
res_hfd_nc_treat_sig <- res_hfd_nc_treat[res_hfd_nc_treat$padj < 0.1 & res_hfd_nc_treat$baseMean > 5.0,]
res_hfd_nc_treat_ord <- res_hfd_nc_treat_sig[order(res_hfd_nc_treat_sig$padj),]
res_hfd_nc_treat_ord$ext_gene <- anno[row.names(res_hfd_nc_treat_ord), "gene_name"]

write.csv(x = res_hfd_nc_treat_ord[,my_cols], file = "de_genes_treated_highfat_vs_chow_padj_10percent.csv")

png("vol_highfat_vs_chow_withNR.png", 1200, 1500, pointsize=20, res=100)
volcanoplot(res_hfd_nc_treat_ord, main = "", lfcthresh=1.0, sigthresh=0.05, textcx=.5, xlim=c(-4, 4), ylim = c(4,8))
dev.off()
