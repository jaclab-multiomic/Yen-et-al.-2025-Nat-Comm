library(Seurat)
library(Signac)
library(ggplot2)
library(dplyr)
library(VennDiagram)
library(RColorBrewer)
library(clusterProfiler) # clusterProfiler_4.6.2
library(org.Mm.eg.db) # org.Mm.eg.db_3.16.0
library(cowplot)
library(VennDiagram)
library(RColorBrewer)


# color setting
l <- c("#AC83B3","#A9ACCC","#A1CFEF", "#94D8CA", "#C5EAB1", "#FEF398","#F9A8A6","#B1BDF2")


# load preprocess data
mtl <- readRDS("./mettl_integrated_slim.rds")


## Fig. 7C
# DefaultAssay(mtl) <- "integrated" ## Integrated RNA by MNN-CCA
# mtl <- RunUMAP(mtl, dims = 1:60, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')
# mtl <- RunUMAP(mtl, reduction = 'integrated_lsi', dims = 2:10, reduction.name = "umap.atac", reduction.key = "atacUMAP_")

p1 <- DimPlot(mtl, reduction = "umap.rna",  cols = l, pt.size = 0.5, shuffle = T, label = T, repel = TRUE, label.size = 4) + ggtitle("RNA")
p2 <- DimPlot(mtl, reduction = "umap.atac", cols = l, pt.size = 0.5, shuffle = T, label = T,  repel = TRUE, label.size = 4) + ggtitle("ATAC")
p3 <- DimPlot(mtl, reduction = "wnn.umap",  cols = l, pt.size = 0.5, shuffle = T, label = T, repel = TRUE, label.size = 4) + ggtitle("Integrated ATAC +RNA (WNN)")
p1 + p2 + p3 & NoLegend() & theme(plot.title = element_text(hjust = 0.5))

#ggsave("celltype_assays.pdf", width = 12, height = 6, units = "in")

## Fig. 7D
# donut plot
df <- as.data.frame(prop.table(table(Idents(mtl), mtl$genotype), margin = 2))
df$percent <- round(df$Freq *100, 2)
names(df) <- c("Cell_type", "Genotype", "Freq", "Percent")
p <- ggplot(df, aes(fill=Cell_type, y=Percent, x=Genotype)) + 
  geom_bar(position="stack", stat="identity") + theme_classic()
pp <- ggplot(df, aes(x = 3,
                     y = Percent,
                     fill = Cell_type)) +
  geom_col(width = 1.5,
           color = 'white') +
  facet_grid(.~Genotype) +
  coord_polar(theta = "y") +
  xlim(c(0.2, 3.8)) +
  scale_fill_manual(values = l) +
  theme_void()+
  theme(
    strip.text.x = element_text(size = 14),
    legend.title = element_text(size = 15),
    legend.text = element_text(size = 14)
  )
pp
pp1 <- pp +
  geom_text(aes(label = ifelse(Percent>3, paste0(round(Percent,2),'%'), '')),
            position = position_stack(vjust = 0.5),
            size = 4)
pp1
# ggsave("celltype_proportion.pdf", width = 7.5, height = 7.5, units = "in")


## Fig. 7E, 7F GO analysis
# Fig. S9A, D
# DEG
table(mtl@active.ident)
DefaultAssay(mtl) <- "RNA"
mtl.markers <- FindAllMarkers(mtl, only.pos = TRUE, min.pct = 0.5, logfc.threshold = 0.5)
mtl.markers %>%
  group_by(cluster) %>%
  slice_max(n = 1e5, order_by = avg_log2FC) -> mtl.markers
write.csv(mtl.markers,file='celltype_markers.csv')

mtl.markers %>%
  group_by(cluster) %>%
  slice_max(n = 20, order_by = avg_log2FC) -> top10
print(top10,n=20*length(unique(mtl@active.ident)))

obj <- mtl; obj <- ScaleData(obj,features=row.names(obj))

# for each subtype
ALL <- NULL
p_list <- list()
dir.create('DEG_heatmap')
celltypes <- names(table(mtl@active.ident))
for(j in 1:length(celltypes)){
  obj.j <- subset(obj,subset= (annotated_cluster == celltypes[j]))
  Idents(obj.j) <- obj.j$genotype
  obj.j@active.ident <- factor(x = obj.j@active.ident, levels =c('WT','KO'))
  DE.j <-FindMarkers(obj.j, ident.1="KO",ident.2="WT",test.use="bimod",
                     assay="RNA") # latent.vars=c("Age","Sex"),
  DE.j$celltype <- celltypes[j]
  DE.j$gene<-rownames(DE.j)
  ALL <- rbind(ALL,DE.j)
  
  # plot heatmap
  library(ggplot2)
  DE.j %>%
    slice_max(n = 30, order_by = avg_log2FC) -> top10
  p_list[[j]] <- DoHeatmap(obj.j, features = top10$gene) + ggtitle(celltypes[j]) #+NoLegend()
  ggsave(paste(paste('DEG_heatmap/DEG_heatmap_',celltypes[j],sep=''),'.pdf',sep='_'),
         plot =   p_list[[j]],
         width=25,height=50*length(top10$gene)/90+10,units='cm')
}

ALL$cat<-ifelse(ALL$avg_log2FC >0, "up","down")
significant<-ALL[which(ALL$p_val_adj<0.05),]
write.csv(significant, "DEGs_bimod_KOCtrl_subtypes.csv")


subc_genes <- split(significant$gene, list(significant$cat, significant$celltype))

# GO
# the processed csv will be uploaded due to version differences

# m6a list from nanopore analysis
m6a <- read.csv("./epinano_m6anet_union.gene 231020 new update.csv")
m6alist <- unique(m6a$symbol)

subc_genes_m6a <- lapply(subc_genes, function(x) {
  intersect(x, m6alist)}) # intersect degs with m6a
chINs_down <- intersect(subc_genes$`down.Cholinergic INs`, m6alist)
skmns_down <- intersect(subc_genes$`down.Skeletal MNs`, m6alist)
vis_down <- intersect(subc_genes$`down.Visceral MNs`, m6alist)

chINs_up <- intersect(subc_genes$`up.Cholinergic INs`, m6alist)
skmns_up <- intersect(subc_genes$`up.Skeletal MNs`, m6alist)
vis_up <- intersect(subc_genes$`up.Visceral MNs`, m6alist)

# also run for MF, CC
ontology <- "BP"
outTitle <- paste0("clusterProfiler_GO-", ontology, "_simplify")
outTitle

# upregulated DEGs, m6a-tagged
IDs <- list(Upregulated_cINs = chINs_up, Upregulated_sMNs = skmns_up, Upregulated_vMNs = vis_up)
gl <- lapply(IDs, function(x) {
  bitr(x,fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Mm.eg.db")})
IDs <- unlist(gl)

ego <- enrichGO(gene      = IDs,
                OrgDb         = org.Mm.eg.db,
                ont           = ontology,
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.05)
ego <- setReadable(ego, OrgDb = org.Mm.eg.db)
upSimGO = clusterProfiler::simplify(ego, cutoff = 0.7, by = "p.adjust", select_fun = min, measure = "Wang", 
                                    semData = NULL)
# Only keep entries with p.adjust < 0.05
upSimGO_padj <- filter(upSimGO@result, p.adjust < 0.05)
upSimGO_padj <- upSimGO_padj[order(upSimGO_padj$Count, decreasing = TRUE),]
upSimGO_top10 <- upSimGO_padj[1:10, ]
ggplot(data=upSimGO_top10, aes(x=Count, y=reorder(Description, Count))) +
  geom_bar(stat="identity", fill="#d1cfe2", width = 0.5) +
  labs(y = " ", x = "Number of genes") +
  ggtitle(paste0("GO-", ontology," upregulated and m6a-tagged genes")) +
  ylab(" ") + xlab("Count") + 
  theme_cowplot() +
  theme_minimal(
    base_size = 30) + 
  theme_cowplot()
ggsave(paste0(outTitle, "up_m6a_chneurons.pdf"), width = 10, height = 5, units = "in")


# Add a new score (qscore) that is -log10(p.adjust)
upSimGO_padj_plot <- mutate(upSimGO_padj, logpval = -log(pvalue, base = 10))
upSimGO_padj_plot <-upSimGO_padj_plot[order(upSimGO_padj_plot$logpval, decreasing = T),]
upSimGO_padj_plot_top10 <- upSimGO_padj_plot[1:10, ]
ggplot(data=upSimGO_padj_plot_top10, aes(x=logpval, y=reorder(Description, logpval))) +
  geom_bar(stat="identity", fill="#d1cfe2", width = 0.5) +
  labs(y = " ", x = "-log10(PValue)") +
  ggtitle(paste0("GO-", ontology," upregulated and m6a-tagged genes")) +
  theme_cowplot() +
  theme_minimal(
    base_size = 30) + 
  theme_cowplot()
ggsave(paste0(outTitle, "up_m6a_chneurons_logpval.pdf"), width = 10, height = 5, units = "in")

# save top10 terms by count/-logpval
write.csv(upSimGO_padj_plot_top10, file = paste0(outTitle, "logpval_top10_up_m6a_chneurons.csv"), quote = F, 
          row.names = F)
write.csv(upSimGO_top10, file = paste0(outTitle, "count_top10_up_m6a_chneurons.csv"), quote = F, 
          row.names = F)

# downregulated genes
IDs_downs <- list(Downregulated_cINs = chINs_down, Downregulated_sMNs = skmns_down, Downregulated_vMNs = vis_down)
gl_down <- lapply(IDs_downs, function(x) {
  bitr(x,fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Mm.eg.db")})
IDs_downs <- unlist(gl_down)
ego_down <- enrichGO(gene      = IDs_downs,
                     OrgDb         = org.Mm.eg.db,
                     ont           = ontology,
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 0.05,
                     qvalueCutoff  = 0.05)
ego_down <- setReadable(ego_down, OrgDb = org.Mm.eg.db)
downSimGO = clusterProfiler::simplify(ego_down, cutoff = 0.7, by = "p.adjust", select_fun = min, measure = "Wang", 
                                      semData = NULL)

# Only keep entries with p.adjust < 0.05
downSimGO_padj <- filter(downSimGO@result, p.adjust < 0.05)
downSimGO_padj <- downSimGO_padj[order(downSimGO_padj$Count, decreasing = TRUE),]
downSimGO_top10 <- downSimGO_padj[1:10, ]
ggplot(data=downSimGO_top10, aes(x=Count, y=reorder(Description, Count))) +
  geom_bar(stat="identity", fill="#d1cfe2", width = 0.5) +
  labs(y = " ", x = "Number of genes") +
  ggtitle(paste0("GO-", ontology," downregulated and m6a-tagged genes")) +
  ylab(" ") + xlab("Count") + 
  theme_cowplot() +
  theme_minimal(
    base_size = 30) + 
  theme_cowplot()
ggsave(paste0(outTitle, "down_m6a_chneurons.pdf"), width = 10, height = 5, units = "in")


# Add a new score (qscore) that is -log10(p.adjust)
downSimGO_padj_plot <- mutate(downSimGO_padj, logpval = -log(pvalue, base = 10))
downSimGO_padj_plot <-downSimGO_padj_plot[order(downSimGO_padj_plot$logpval, decreasing = T),]
downSimGO_padj_plot_top10 <- downSimGO_padj_plot[1:10, ]
ggplot(data=downSimGO_padj_plot_top10, aes(x=logpval, y=reorder(Description, logpval))) +
  geom_bar(stat="identity", fill="#d1cfe2", width = 0.5) +
  labs(y = " ", x = "-log10(PValue)") +
  ggtitle(paste0("GO-", ontology," downregulated and m6a-tagged genes")) +
  theme_cowplot() +
  theme_minimal(
    base_size = 30) + 
  theme_cowplot()
ggsave(paste0(outTitle, "down_m6a_chneurons_logpval.pdf"), width = 10, height = 5, units = "in")

# save top10 terms by count/-logpval
write.csv(downSimGO_padj_plot, file = paste0(outTitle, "logpval_top10_down_m6a_chneurons.csv"), quote = F, 
          row.names = F)
write.csv(downSimGO_padj, file = paste0(outTitle, "count_top10_down_m6a_chneurons.csv"), quote = F, 
          row.names = F)


# KEGG 
kegg_down <- enrichKEGG(gene = IDs_downs, organism = "mmu", pvalueCutoff = 0.05, 
                        pAdjustMethod = "BH", qvalueCutoff  = 0.05)
# kegg_up_padj <- filter(kegg_up@result, p.adjust < 0.05)
kegg_down <- kegg_down@result
kegg_down$Descriptions <- gsub(" - Mus musculus \\(house mouse)$", "", as.character(kegg_down$Description))
kegg_down_padj <- filter(kegg_down, p.adjust < 0.05)
kegg_down_padj <- kegg_down_padj[order(kegg_down_padj$Count, decreasing = TRUE),]
# kegg_down_top10 <- kegg_down_padj[1:10, ]
ggplot(data=kegg_down_padj, aes(x=Count, y=reorder(Descriptions, Count))) +
  geom_bar(stat="identity", fill="#d1cfe2", width = 0.5) +
  labs(y = " ", x = "Number of genes") +
  ggtitle(paste0("KEGG - downregulated and m6a-tagged genes")) +
  ylab(" ") + xlab("Count") + 
  theme_cowplot() +
  theme_minimal(
    base_size = 30) + 
  theme_cowplot()
ggsave(paste0("KEGG_down_m6a_chneurons.pdf"), width = 10, height = 5, units = "in")

# Add a new score (qscore) that is -log10(p.adjust)
kegg_down_padj_plot <- mutate(kegg_down_padj, logpval = -log(pvalue, base = 10))
kegg_down_padj_plot <-kegg_down_padj_plot[order(kegg_down_padj_plot$logpval, decreasing = T),]
# kegg_down_padj_plot_top10 <- kegg_down_padj_plot[1:10, ]
ggplot(data=kegg_down_padj_plot, aes(x=logpval, y=reorder(Descriptions, logpval))) +
  geom_bar(stat="identity", fill="#d1cfe2", width = 0.5) +
  labs(y = " ", x = "-log10(PValue)") +
  ggtitle(paste0("KEGG - downregulated and m6a-tagged genes")) +
  theme_cowplot() +
  theme_minimal(
    base_size = 30) + 
  theme_cowplot()
ggsave(paste0("KEGG_down_m6a_chneurons_logpval.pdf"), width = 10, height = 5, units = "in")
write.csv(kegg_down_padj_plot, file = "KEGG-m6a-tagged_genes.csv", quote = F, 
          row.names = F)

kegg_down_padj %>% dplyr::select(ID, geneID) -> kegg2
kegg2 <- strsplit(kegg2$geneID, "/") 
names(kegg2) <- kegg_down_padj$ID
kegg4 <- lapply(
  kegg2, function(x) {
    bitr_conv <- bitr(x, fromType = "ENTREZID", toType = "SYMBOL", org.Mm.eg.db)
    # Symbols are stored in column "SYMBOL", so we extract that for intersection
    genes <- bitr_conv$SYMBOL
    return(genes)
  }
)
longest <- max(
  sapply(kegg4, function(x) {length(x)})
)
kegg4 <- lapply(kegg4, function(x) {
  pad <- rep("", longest - length(x))
  return(c(x, pad))
})
kegg4_df <- do.call(cbind, kegg4)
write.csv(kegg4_df, "kegg_down_padj_cutoff_symbol.csv")


## Fig. 7F
# m6a list from Nanopore analysis
m6a <- read.csv("./epinano_m6anet_union.gene 231020 new update.csv")
m6alist <- unique(m6a$symbol)


## Fig. 7H
DefaultAssay(mtl) <- "RNA"
Idents(mtl) <- "celltype_genotype"

mtl@active.ident <- 
  factor(mtl@active.ident, levels=c("WT_Skeletal MNs", "KO_Skeletal MNs", 
                                    "WT_Visceral MNs", "KO_Visceral MNs", 
                                    "WT_Cholinergic INs", "KO_Cholinergic INs",
                                    "WT_Excitatory neurons", "KO_Excitatory neurons",
                                    "WT_Inhibitory Neurons", "KO_Inhibitory Neurons",
                                    "WT_Astrocytes", "KO_Astrocytes", 
                                    "WT_Oligodendrocytes", "KO_Oligodendrocytes"))


genes <- c("Fus", "Bscl2", "Erbb4", "Pikfyve", "Dctn1")
DotPlot(mtl, genes, idents = c("WT_Skeletal MNs", "KO_Skeletal MNs", 
                               "WT_Visceral MNs", "KO_Visceral MNs", 
                               "WT_Cholinergic INs", "KO_Cholinergic INs")) + 
  coord_flip() + RotatedAxis() + 
  scale_colour_gradient2(low = "red", mid = "grey", high = "blue")
#ggsave("m6atagged_multiome_genes.pdf")


## Fig. S7 
# A. QC metrics
Idents(mtl) <- "orig.ident"
mtl@active.ident <- factor(mtl@active.ident, levels=c("KO_r1", "KO_r3",  "KO_r4", "WT_r3", "WT_r4", "WT_r5"))
VlnPlot(mtl, features = c("nCount_RNA", "nFeature_RNA", "percent.mt"), pt.size = 0, cols = l, ncol=1)
#ggsave("qcmetrics.pdf", width = 7, height = 14, units = "in")
DefaultAssay(mtl) <- "shared_ATAC"
VlnPlot(mtl, features = c("atac_fragments", "TSS.enrichment", "nucleosome_signal", "blacklist_fraction"), pt.size = 0, cols = l, ncol=1)
#ggsave("qcmetrics_atacII.pdf", width = 7, height = 14, units = "in")


# B and C. UMAP group.by samples/genotype and batch
color <- c("#A083D5", "#F4E064")
p1 <- DimPlot(mtl, group.by = "genotype", shuffle = T, pt.size = 0.5, cols = color, reduction = "wnn.umap")
p2 <- DimPlot(mtl, shuffle = T, pt.size = 0.5, cols = l, reduction = "wnn.umap")
p1 + p2
#ggsave("replicates.pdf", width = 14, height = 7, units = "in")


# D. Pca of cell cluster replicates
# subset cholinergic INs, skeletal MNs and visceral MNs
ch <- readRDS("./cholinergicneurons_Mettl14.rds")
DefaultAssay(ch) <- "RNA"

ch_avg2 <- AverageExpression(object = ch, 
                             group.by = c('annotated_cluster', 'orig.ident'), 
                             return.seurat = TRUE)
ch_avg2 <- FindVariableFeatures(ch_avg2)
ch_avg2 <- ScaleData(ch_avg2, assay = "RNA")
ch_avg2 <- RunPCA(ch_avg2, npcs = 17)

# PC1 and PC2
r1 <- DimPlot(ch_avg2, pt.size = 1.5, reduction = "pca",
        cols = c("Cholinergic INs_KO-r1" = "#AC83B3" , 
                 "Cholinergic INs_KO-r3" = "#AC83B3", 
                 "Cholinergic INs_KO-r4" = "#AC83B3", 
                 "Cholinergic INs_WT-r3" = "#A9ACCC", 
                 "Cholinergic INs_WT-r4" = "#A9ACCC", 
                 "Cholinergic INs_WT-r5" = "#A9ACCC", 
                 "Skeletal MNs_KO-r1" = "#A1CFEF", 
                 "Skeletal MNs_KO-r3" = "#A1CFEF", 
                 "Skeletal MNs_KO-r4" = "#A1CFEF", 
                 "Skeletal MNs_WT-r3" = "#94D8CA", 
                 "Skeletal MNs_WT-r4" = "#94D8CA", 
                 "Skeletal MNs_WT-r5" = "#94D8CA",
                 "Visceral MNs_KO-r1" = "#C5EAB1", 
                 "Visceral MNs_KO-r3"= "#C5EAB1", 
                 "Visceral MNs_KO-r4"= "#C5EAB1", 
                 "Visceral MNs_WT-r3" = "#FEF398",
                 "Visceral MNs_WT-r4"= "#FEF398", 
                 "Visceral MNs_WT-r5"= "#FEF398"))

# PC2 and PC3
r2 <- DimPlot(ch_avg2, pt.size = 1.5, reduction = "pca", dims = c(2:3),
        cols = c("Cholinergic INs_KO-r1" = "#AC83B3" , 
                 "Cholinergic INs_KO-r3" = "#AC83B3", 
                 "Cholinergic INs_KO-r4" = "#AC83B3", 
                 "Cholinergic INs_WT-r3" = "#A9ACCC", 
                 "Cholinergic INs_WT-r4" = "#A9ACCC", 
                 "Cholinergic INs_WT-r5" = "#A9ACCC", 
                 "Skeletal MNs_KO-r1" = "#A1CFEF", 
                 "Skeletal MNs_KO-r3" = "#A1CFEF", 
                 "Skeletal MNs_KO-r4" = "#A1CFEF", 
                 "Skeletal MNs_WT-r3" = "#94D8CA", 
                 "Skeletal MNs_WT-r4" = "#94D8CA", 
                 "Skeletal MNs_WT-r5" = "#94D8CA",
                 "Visceral MNs_KO-r1" = "#C5EAB1", 
                 "Visceral MNs_KO-r3"= "#C5EAB1", 
                 "Visceral MNs_KO-r4"= "#C5EAB1", 
                 "Visceral MNs_WT-r3" = "#FEF398",
                 "Visceral MNs_WT-r4"= "#FEF398", 
                 "Visceral MNs_WT-r5"= "#FEF398"))

r1 + r2       
#ggsave("pca_clusters_replicates.pdf", width = 14, height = 7, units = "in")           


# E. Violinplot of known markers
DefaultAssay(mtl) <- "RNA"
Idents(mtl) <- "annotated_cluster"
mtl@active.ident <- factor(mtl@active.ident, 
                           levels=c("Skeletal MNs", "Visceral MNs",  "Cholinergic INs",
                                    "Excitatory neurons", "Inhibitory Neurons", "Oligodendrocytes", "Astrocytes"))
genes3 <- c("Slc5a7","Chat", "Tns1", "Bcl6", "Nos1", "Pax2", 
            "Slc17a6", "Gad1", "Gad2", "Slc6a5", 
            "Mbp", "Mobp", "Gfap", "Aqp4")
VlnPlot(mtl, features = genes3, pt.size = 0, stack = T, flip = T)
#ggsave("vlnplot.pdf", width = 7.5, height = 7.5, units = "in")



## Fig. S8
# A. UMAP
# skeletal motor neurons are subset based on previous label transfer annotation. 
# 14 cells were dropped based on their low label transfer prediction score from Blum et al annotation
# some low prediction scores cell were annotated as "undetermined" later
# note: Seurat version used were 2.3.4, label transfer prediction may differ with later version
# [1] "KO_r1_CGGTTCCGTTGGGTTA-1" "KO_r1_TGTGGCTCATTCAGCA-1" "KO_r4_AATTGGACATCGCTCC-1"
# [4] "KO_r4_GGTAAGGGTATCTGGA-1" "KO_r4_TGAGCTTAGCCTGAGC-1" "WT_r3_AAGACAAGTTGTAACG-1"
# [7] "WT_r3_TGAAGCAAGCGGATAA-1" "WT_r4_GCTAGCCAGGCTATGT-1" "WT_r4_GGCCATCAGGAGGACT-1"
# [10] "WT_r5_AGTGTGGCATAGGCGA-1" "WT_r5_CGAACAAAGGCTGTCA-1" "WT_r5_CTATGACAGAGAGGCT-1"
# [13] "WT_r5_GACGTAAAGGCGCTAC-1" "WT_r5_TGTTGTAAGGCACAGG-1"

sk <- readRDS("./skeletalmns.rds")
DefaultAssay(sk) <- "RNA"
sk <- NormalizeData(sk)
sk <- FindVariableFeatures(sk)

# load ref object (Blum et al)
load("./skeletalmns_Blum.RData")
refobj <- NormalizeData(skeletalmns)
refobj <- FindVariableFeatures(refobj)

# cell label transfer
ann_anchors <- FindTransferAnchors(
  reference = refobj, query = sk, normalization.method = "LogNormalize"
)
pred <- TransferData(ann_anchors, refdata = refobj$subcluster2)
sk <- AddMetaData(sk, pred)
head(pred)

sk$pred.scores <- "false"
sk$pred.scores[which(sk$prediction.score.max >= 0.85)] <- "true" # tried 0.8, 0.9
sk$pred_id_modify <- sk$predicted.id
sk$pred_id_modify[which(sk$pred.scores == "false")] <- "undetermined"
Idents(sk) <- "pred_id_modify"
sk@active.ident <- 
  factor(sk@active.ident, levels=c("alpha", "gamma", "gamma*", "undetermined"
  ))
DimPlot(sk, cols = l)
# ggsave("umap_predictid_all.pdf", width = 7, height = 7, units = "in")

# B. Violin plot of prediction scores
VlnPlot(sk, features = c("prediction.score.alpha", "prediction.score.gamma", 
                         "prediction.score.gamma.", "prediction.score.max"), 
        ncol = 2, cols = l)
#ggsave("vlnplot_predictid_all.pdf", width = 5, height = 5, units = "in")


# C. Violin plot of reported markers
genes <- c("Rbfox3", "Spp1", "Vipr2", "Npas1", "Htr1d", "Creb5", "Pard3b","Stxbp6", "Plch1")
VlnPlot(sk, genes, stack = T, flip = T) + NoLegend()
# ggsave("vlnplot_markers_subset.pdf", width = 5, height = 7, units = "in")

# D. donut plot
df <- as.data.frame(prop.table(table(Idents(sk), sk$genotype), margin = 2))
df$percent <- round(df$Freq *100, 2)
names(df) <- c("Cell_type", "Genotype", "Freq", "Percent")
p <- ggplot(df, aes(fill=Cell_type, y=Percent, x=Genotype)) + 
  geom_bar(position="stack", stat="identity") + theme_classic()
pp <- ggplot(df, aes(x = 3,
                     y = Percent,
                     fill = Cell_type)) +
  geom_col(width = 1.5,
           color = 'white') +
  facet_grid(.~Genotype) +
  coord_polar(theta = "y") +
  xlim(c(0.2, 3.8)) +
  scale_fill_manual(values = l) +
  theme_void()+
  theme(
    strip.text.x = element_text(size = 14),
    legend.title = element_text(size = 15),
    legend.text = element_text(size = 14)
  )
pp
pp1 <- pp +
  geom_text(aes(label = ifelse(Percent>3, paste0(round(Percent,2),'%'), '')),
            position = position_stack(vjust = 0.5),
            size = 4)
pp1
#ggsave("celltype_proportion.pdf", width = 7.5, height = 7.5, units = "in")


# E and F. subtypes markers
genes1 <- c("Chodl", "Sv2a", "Kcnq5")
FeaturePlot(sk, genes1, cells = WhichCells(sk, idents = "alpha"))
# ggsave("umap_alphasubcluster.pdf", width = 7, height = 5, units = "in")



# B, C, E. selected DE genes
# stackvlnplot with reference to https://rdrr.io/github/sqjin/CellChat/
source("stackvlnplot_edited.R")
ch <- readRDS("./cholinergicneurons_Mettl14.rds")
markers <- read.csv("./DEGs_bimod_KOCtrl_subtypes.csv")
markers <- markers[markers$celltype %in% c("Skeletal MNs", "Cholinergic INs", "Visceral MNs"), ]


DefaultAssay(ch) <- "RNA"
Idents(ch) <- "annotated_cluster"
ch@active.ident <- 
  factor(ch@active.ident, levels=c("Skeletal MNs", "Visceral MNs", "Cholinergic INs"))

genes3 <- c("Malat1", "Srsf2", "Fus", "Srek1", "Hnrnph1", "Tra2a")
StackedVlnPlot(ch, features = genes3, split.by = "genotype", pt.size = 0.000001,
               deg.tbl = markers, alpha = 0.5, title = "RNA splicing related genes")
# ggsave(filename = "RNA-splicing-related-genes.pdf", width = 8, height = 10)

genes4 <- c("Nefm", "Nefl", "Kif5c", "Tubb3", "Fus", "Nup93")
StackedVlnPlot(ch, features = genes4, split.by = "genotype", pt.size = 0.000001,
               deg.tbl = markers, alpha = 0.5, title = "KEGG ALS pathway")
# ggsave(filename = "KEGG-ALS-pathway-genes.pdf", width = 8, height = 10)

genes5 <- c("Chd2", "Chd6", "Bcl7c", "Ncoa6", "Ube2b")
StackedVlnPlot(ch, features = genes5, split.by = "genotype", pt.size = 0.000001,
               deg.tbl = markers, alpha = 0.5, title = "ATP-dependent chromatin remodeler")
# ggsave(filename = "ATP-dependent-chromatin-remodeler-genes.pdf", width = 8, height = 10)




