## credit to Wei ZHAO (zhaow17@uci.edu), UC Irvine 
## adapted from ArchR package
library(Seurat);library(dplyr);library(data.table);library(Matrix)
library(stats);library(reshape2)
library(Signac)
library(Seurat)
library(JASPAR2020)
library(TFBSTools)
library(BSgenome.Mmusculus.UCSC.mm10)
library(EnsDb.Mmusculus.v75)
library(patchwork)
library(Rcpp)
library(GenomicRanges)
library(ggplot2)
library(Rcpp)
sourceCpp("~/Documents/GitHub/CCI2TF2CCI/src/row_correlation.cpp") # download from https://github.com/GreenleafLab/ArchR/blob/master/src/Correlation.cpp
sourceCpp("~/Documents/GitHub/CCI2TF2CCI/src/knn.cpp") # download from https://github.com/GreenleafLab/ArchR/blob/master/src/KNN_Utils.cpp

setwd('/Users/weizhao/Documents/scMultiome/reseq/')
int <- readRDS(file="mettl_integrated_annotate.rds") 

# define functions
myArchR_peak2geneLinkage <- function(seu,peak.assay='shared_ATAC',peak.slot='counts',rna.assay='RNA',rna.slot='data',
                                     corCutOff = 0.45,
                                     FDRCutOff = 1e-04){
  #### peak-gene relation
  # preparing genome info
  library(EnsDb.Mmusculus.v75)
  keys <- keys(EnsDb.Mmusculus.v75)
  anno.result <- AnnotationDbi::select(EnsDb.Mmusculus.v75, keys=keys, columns=c("GENEID","SYMBOL",'SEQNAME',"TXSEQSTART","TXCDSSEQSTART"),keytype="GENEID")
  genome.info <- anno.result[,c(3,4,4,2)]
  names(genome.info) <- c("Chrom","Starts","Ends","genes"); genome.info$Ends <- genome.info$Starts +1
  genome.info <- genome.info[!duplicated(genome.info$genes),] # filter out different transcript
  genome.info$Chrom <- paste('chr',genome.info$Chrom, sep='')
  
  #Get Peak Set
  DefaultAssay(seu) <- peak.assay #'shared_ATAC'
  main.chroms <- standardChromosomes(BSgenome.Mmusculus.UCSC.mm10)
  peak.matrix <- slot(seu[[seu@active.assay]],name = peak.slot) # counts slot
  idx.keep <- rowSums(x = peak.matrix) > 0
  peak.matrix <- peak.matrix[idx.keep, , drop = FALSE]
  # motif.matrix <- motif.matrix[idx.keep, , drop = FALSE]
  peak.ranges <- StringToGRanges(rownames(peak.matrix), sep = c("-", "-"))
  idx.keep2 <- which(seqnames(peak.ranges) %in% standardChromosomes(peak.ranges, species="Mus_musculus"))
  peak.matrix <- peak.matrix[idx.keep2,]
  peak.ranges <- peak.ranges[idx.keep2,]
  mcols(peak.ranges)$peak_id <- rownames(peak.matrix)[idx.keep2]
  peakSet <- peak.ranges
  
  #Gene Info
  RNA_assay_name <- rna.assay #'RNA'
  genome.info <- genome.info[ genome.info$Chrom  %in% extractSeqlevels('Mus_musculus',style = 'UCSC'),]
  rna.matrix <-slot(seu[[RNA_assay_name]], name= rna.slot)
  rna.matrix <- rna.matrix[rowSums(x = slot(seu[[RNA_assay_name]], name= 'counts')) > 100, , drop = FALSE]
  geneSet <- genome.info[which(genome.info$genes %in% rownames(rna.matrix)),]
  rna.matrix <- rna.matrix[geneSet$genes,]
  geneStart <- GRanges(geneSet$Chrom, IRanges(geneSet$Starts, width = 1), name = geneSet$genes, idx = 1:dim(geneSet)[1])
  
  #Check Chromosomes
  chri <- gtools::mixedsort(unique(paste0(seqnames(peakSet))))
  chrj <- gtools::mixedsort(unique(paste0(seqnames(geneStart))))
  chrij <- intersect(chri, chrj)
  
  #Get Reduced Dims
  rD <- seu@reductions$integrated_lsi@cell.embeddings
  #Subsample
  idx <- sample(seq_len(nrow(rD)), 500, replace = !nrow(rD) >= 500)
  #KNN Matrix
  data = rD; query = rD[idx,]; k = min(100,round(1/10*dim(rD)[1]))
  library(nabor)
  knnObj <- nabor::knn(data = data, query = query, k = k)$nn.idx
  overlapCutoff = 0.8
  
  
  keepKnn <- determineOverlapCpp(knnObj, floor(overlapCutoff * k))
  knnObj <- knnObj[keepKnn==0,]
  #Convert To Names List
  knnObj <- lapply(seq_len(nrow(knnObj)), function(x){
    rownames(rD)[knnObj[x, ]]
  }) %>% SimpleList
  
  #Group Matrix RNA
  rna_mask <- sapply(seq_len(length(knnObj)), function(x) colnames(rna.matrix) %in%
                       knnObj[[x]])
  rna_mask <- Matrix::Matrix(rna_mask)
  rna_new <- as.matrix(rna.matrix %*% rna_mask)
  
  #Group Matrix ATAC
  atac_mask <- sapply(seq_len(length(knnObj)), function(x) colnames(peak.matrix) %in%
                        knnObj[[x]])
  atac_mask <- Matrix::Matrix(atac_mask)
  atac_new <- as.matrix(peak.matrix %*% atac_mask)
  
  rna_new <- t(t(rna_new) / colSums(rna_new)) * 1e6
  atac_new <- t(t(atac_new) / colSums(atac_new)) * 1e6
  
  rna_new <- log(rna_new + 1)
  atac_new <- log(atac_new + 1)
  
  #Overlaps
  o <- DataFrame(findOverlaps( resize(geneStart, 2 * 250000 + 1, "center"), peakSet, ignore.strand = TRUE))
  
  #Get Distance from Fixed point A B
  o$distance <- distance(geneStart[o[,1]] , peakSet[o[,2]] )
  colnames(o) <- c("B", "A", "distance")
  
  # Computing Correlations
  o$Correlation <- rowCorCpp(as.integer(o$A), as.integer(o$B), atac_new, rna_new)
  o$VarAssayA <- ArchR:::.getQuantiles(matrixStats::rowVars(atac_new))[o$A]
  o$VarAssayB <- ArchR:::.getQuantiles(matrixStats::rowVars(rna_new))[o$B]
  o$TStat <- (o$Correlation / sqrt((pmax(1-o$Correlation^2, 0.00000000000000001, na.rm = TRUE))/(ncol(atac_new)-2))) #T-statistic P-value
  o$Pval <- 2*pt(-abs(o$TStat), ncol(atac_new) - 2)
  o$FDR <- p.adjust(o$Pval, method = "fdr")
  out <- o[, c("A", "B", "Correlation", "FDR", "VarAssayA", "VarAssayB")]
  colnames(out) <- c("idxATAC", "idxRNA", "Correlation", "FDR", "VarQATAC", "VarQRNA")
  out$gene <- rownames(rna_new)[out$idxRNA]
  out$peak <- rownames(atac_new)[out$idxATAC]
  out <- out[,c(1,8,2,7,3:6)]
  out <- out[,c(2,4,5,6)]
  # corCutOff = 0.45; FDRCutOff = 1e-04
  #out <- out[which(out$Correlation>=corCutOff & out$FDR <= FDRCutOff & out$VarQATAC >=varCutOffATAC & out$VarQRNA >= varCutOffRNA),]
  out <- out[which(out$Correlation>=corCutOff & out$FDR <= FDRCutOff),]
  return(out)}
lm_eqn2 <- function(df){
  m <- lm(logFC_gene ~ logFC_peak, df);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3)))
  return(as.character(as.expression(eq)))
}

# for loop for cell types
celltypes <- names(table(int@active.ident)); 
dir.create('ArchR_peak2geneLinkage_DEG')
deg_list1 <- read.csv(file='volcano_plot/DEGs_bimod_KOCtrl_all_genes_subtypes.csv')
deg_list <- split(deg_list1, f = deg_list1$celltype)     
deg_list[[8]] <- read.csv(file='volcano_plot/DEGs_bimod_KOCtrl_all_genes_cholinergic.csv')
deg_list[[9]] <- read.csv(file='volcano_plot/DEGs_bimod_KOCtrl_all_genes_all_cells.csv')
names(deg_list)[8:9] <- celltypes[8:9]
deg_list <- deg_list[celltypes]

p_corr_list <- list()
p_corr_list0 <- list()

for(j in 1:length(celltypes)){
seu <- subset(int, subset=annotated_cluster %in% celltypes[j])#,"Visceral MNs","Cholinergic INs"))
out <- myArchR_peak2geneLinkage(seu)
## Differential peaks
da_peaks <- presto:::wilcoxauc.Seurat(seu, group_by = "genotype", assay = 'data', seurat_assay = "shared_ATAC", verbose= TRUE)
da_peaks <- da_peaks[da_peaks$group=='KO',]
dplyr::filter(da_peaks, padj < 0.01, logFC > 0.1, pct_in > 5) %>% arrange(group, -auc) -> da_peaks_tmp
## Differential  genes
de_genes<- presto:::wilcoxauc.Seurat(seu, group_by = "genotype", assay = 'data', seurat_assay = "RNA", verbose= TRUE)
de_genes <- de_genes[de_genes$group=='KO',]

out$logFC_gene <- de_genes$logFC[match(out$gene,de_genes$feature)]
out$logFC_peak <- da_peaks$logFC[match(out$peak,da_peaks$feature)]
out <- as.data.frame(out)

deg_j <- deg_list[[celltypes[j]]]$gene[deg_list[[celltypes[j]]]$sig == TRUE]
dap_j <- da_peaks_tmp$feature

out$DEG <- ifelse(out$gene %in% deg_j,TRUE,FALSE)
out$DAP <- ifelse(out$peak %in% dap_j,TRUE,FALSE)

write.csv(out,file=paste('./ArchR_peak2geneLinkage_DEG/','ArchR_peak2geneLinkage_',celltypes[j],'.csv',sep=''))

rm(seu);gc()
# scatter plot to show correlation
out %>% dplyr::filter(DEG==TRUE) %>% 
  group_by(gene) %>%
  slice_max(n=1, order_by = Correlation) -> hhd
cor_spearman <- cor.test(hhd$logFC_gene,hhd$logFC_peak,method = 'spearman')
library(cowplot)
library(scattermore)
library(ggrastr)

p_corr <-  ggplot(data = hhd, aes(x =logFC_peak, y = logFC_gene)) + theme_cowplot() +
  rasterise(geom_point(), dpi = 300) +
  geom_smooth(method = "lm", se = FALSE) +
  labs(title = celltypes[j], x = "logFC_peak", y = "logFC_gene") +
  theme(text = element_text(size = 12), axis.text = element_text(size = 12), legend.text = element_text(size=12))
# p_corr <- ggplot(data=hhd,aes(x =logFC_peak, y = logFC_gene)) +  theme_cowplot() +
#   geom_scattermore(aes(logFC_peak, logFC_gene),
#                    pointsize = 3,
#                    alpha = 1,
#                    pixels = c(512, 512),
#                    interpolate = TRUE) + geom_smooth(method = "lm", se = FALSE) +
#   labs(title = celltypes[j], x = "logFC_peak", y = "logFC_gene") + 
#   theme(text = element_text(size = 18), axis.text.x = element_text(size = 18), legend.text = element_text(size=18))
#   
p_corr1 <- p_corr +  geom_text(x = ggplot_build(p_corr)$layout$panel_scales_x[[1]]$range$range[1], 
                               y = ggplot_build(p_corr)$layout$panel_scales_y[[1]]$range$range[2]-0.1, hjust=0,
                               label = lm_eqn2(hhd), parse = TRUE) + 
  geom_text(x = ggplot_build(p_corr)$layout$panel_scales_x[[1]]$range$range[1],
            y = ggplot_build(p_corr)$layout$panel_scales_y[[1]]$range$range[2], hjust=0,
            label = paste('Spearman rho = ',round(cor_spearman$estimate,digits = 3),', p-value = ', signif(cor_spearman$p.value,digits = 3), sep=''), parse = FALSE) 
p_corr_list[[j]] <- p_corr1

p_corr <- p_corr + theme(axis.text.x=element_blank(),
                         axis.text.y=element_blank(),
                         axis.title.x=element_blank(),
                         axis.title.y=element_blank(),legend.position="none",
                         panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
                         panel.grid.minor=element_blank(),plot.background=element_blank())
p_corr_list0[[j]] <- p_corr
# ggsave(paste('./ArchR_peak2geneLinkage_DEG/scatter_gene_peak_correlation_',celltypes[j],'.tiff',sep=''),
#        plot =   p_corr1,
#        width=10,height=10,units='cm')
# ggsave(paste('./ArchR_peak2geneLinkage_DEG/scatter_gene_peak_correlation_',celltypes[j],'.pdf',sep=''),
#        plot =   p_corr1,
#        width=10,height=10,units='cm')

}
for(j in 1:length(p_corr_list)){
  ggsave(paste('./ArchR_peak2geneLinkage_DEG/scatter_gene_peak_correlation_',celltypes[j],'.pdf',sep=''),
         #plot =  rasterize(p_corr_list[[j]], layers='Point', dpi=50),
         plot =  p_corr_list[[j]],
         width=10,height=10,units='cm')
  ggsave(paste('./ArchR_peak2geneLinkage_DEG/scatter_gene_peak_correlation_',celltypes[j],'_notext.pdf',sep=''),
         #plot =  rasterize(p_corr_list[[j]], layers='Point', dpi=50),
         plot =  p_corr_list0[[j]],
         width=10,height=10,units='cm')
}
# ggsave(paste('./ArchR_peak2geneLinkage_DEG/scatter_gene_peak_correlation_combined','.tiff',sep=''),#dpi = 300,
#          plot =   wrap_plots(p_corr_list,ncol = 3),
#          width=30,height=20,units='cm')
ggsave(paste('./ArchR_peak2geneLinkage_DEG/scatter_gene_peak_correlation_combined_text','.pdf',sep=''),#dpi = 300,
       plot =   wrap_plots(p_corr_list,ncol = 3),
       width=30,height=20,units='cm')
ggsave(paste('./ArchR_peak2geneLinkage_DEG/scatter_gene_peak_correlation_combined_notext','.pdf',sep=''),#dpi = 300,
       plot =   wrap_plots(p_corr_list0,ncol = 3),
       width=30,height=20,units='cm')
