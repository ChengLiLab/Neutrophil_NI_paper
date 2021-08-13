# Samples: wt.ctl.bm1, wt.ctl.bm2, wt.ctl.pb2 and wt.ctl.sp2
# Aim: select neu by basic analysis
# 
# 2020-03-20, Qiang Shi

rm(list = ls())
setwd("/lustre/user/liclab/shiq/single_cell/luolab/neu_profile/3_wt_ctl/1_ref")

library(Seurat, lib.loc = "/lustre/user/liclab/shiq/software/R_library/")
library(dplyr)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
library(RColorBrewer)
library(ggsci)

########################################
# (1) Setup the Seurat Object
#####################################
# 1) merge single sample passing QC
wt.ctl.bm1 <- readRDS("../0_QC/wt_ctl_bm1/wt_ctl_bm1_qc.rds")
wt.ctl.bm2 <- readRDS("../0_QC/wt_ctl_bm2/wt_ctl_bm2_qc.rds")
wt.ctl.pb2 <- readRDS("../0_QC/wt_ctl_pb2/wt_ctl_pb2_qc.rds")
wt.ctl.sp2 <- readRDS("../0_QC/wt_ctl_sp2/wt_ctl_sp2_qc.rds")

tmp1 <- MergeSeurat(object1 = wt.ctl.bm1, 
                    object2 = wt.ctl.bm2, 
                    add.cell.id1 = 'wt.ctl.bm1', 
                    add.cell.id2 = 'wt.ctl.bm2',
                    do.normalize = FALSE)

tmp2 <- MergeSeurat(object1 = tmp1,
                    object2 = wt.ctl.pb2, 
                    add.cell.id2 = 'wt.ctl.pb2',
                    do.normalize = FALSE)

ref <- MergeSeurat(object1 = tmp2, 
                   object2 = wt.ctl.sp2, 
                   add.cell.id2 = 'wt.ctl.sp2',
                   do.normalize = FALSE,
                   project      = 'ref')

ref@raw.data@Dim
# 27998 19582

ref@data@Dim
# 27998 19582

head(ref@meta.data)

table(ref@meta.data$orig.ident)
# wt_ctl_bm1 wt_ctl_bm2 wt_ctl_pb2 wt_ctl_sp2
# 3063       5240       5655       5624

# 2) add tissue info
ref@meta.data$tissue_marker <- ref@meta.data$orig.ident
ref@meta.data$tissue        <- ref@meta.data$orig.ident

old.id  <- c('wt_ctl_bm1', 'wt_ctl_bm2', 'wt_ctl_pb2', 'wt_ctl_sp2')
new.id1 <- c('WT_BM_C-kit+Gr1', 'WT_BM_Gr1', 'WT_PB_Gr1', 'WT_SP_Gr1')
new.id2 <- c('WT_BM_C-kit+Gr1', 'WT_BM', 'WT_PB', 'WT_SP')

ref@meta.data$tissue_marker <- plyr::mapvalues(x    = ref@meta.data$tissue_marker, 
                                               from = old.id, 
                                               to   = new.id1)

ref@meta.data$tissue <- plyr::mapvalues(x    = ref@meta.data$tissue, 
                                        from = old.id, 
                                        to   = new.id2)

table(ref@meta.data$tissue_marker)
table(ref@meta.data$tissue)

# renew factor order
ref@meta.data$tissue_marker <- factor(ref@meta.data$tissue_marker)
ref@meta.data$tissue        <- factor(x      = ref@meta.data$tissue, 
                                      levels = new.id2)

# 3) add batch info
ref@meta.data$batch <- 'Batch2'

ref@meta.data$batch[ref@meta.data$orig.ident=='wt_ctl_bm1'] <- 'Batch1'

ref@meta.data$batch <- factor(ref@meta.data$batch)

table(ref@meta.data$batch)

head(ref@meta.data)

# 4) save object
ref <- SetAllIdent(ref, id = 'tissue')
saveRDS(ref, file = "./ref_raw.rds")
# ref <- readRDS("./ref_raw.rds")

rm(wt.ctl.bm1, wt.ctl.bm2, wt.ctl.pb2, wt.ctl.sp2, tmp1, tmp2)
rm(old.id, new.id1, new.id2)

#####


########################################
# (2) QC
#####################################
pdf("./ref_qc.pdf")
VlnPlot(object         = ref, 
        features.plot  = c('nGene', 'nUMI', 'percent.mt', 'nUMI.per.Gene'),
        nCol           = 2,
        size.x.use     = 10,
        point.size.use = 0,
        cols.use       = brewer.pal(12, 'Paired')[c(1,2,4,8)],
        group.by       = 'tissue_marker',
        x.lab.rot      = TRUE)
dev.off()

aggregate(nGene ~ tissue_marker, data = ref@meta.data, FUN = "median")
#   tissue_marker nGene
# WT_BM_C-kit+Gr1  2421
#       WT_BM_Gr1  1375
#       WT_PB_Gr1   772
#       WT_SP_Gr1  1541

median(ref@meta.data$nGene)
# 1241

# GenePlot 
GenePlot(object = ref, gene1 = 'nUMI', gene2 = 'percent.mt')
GenePlot(object = ref, gene1 = 'nUMI', gene2 = 'nGene')

#####


########################################
# (3) Normalization, Feature, Scale
#####################################
# normalization
ref <- NormalizeData(ref)

# variable genes
ref <- FindVariableGenes(object        = ref, 
                         x.low.cutoff  = 0.0125, 
                         x.high.cutoff = 4, 
                         y.cutoff      = 1)

length(x = ref@var.genes)
# 2553

write.table(ref@var.genes, "./ref_vargene.txt", row.names = FALSE)

# scale
head(ref@meta.data)
# no.regress has the best clusters
ref <- ScaleData(object = ref, #genes.use = ref@var.genes,
                 display.progress = FALSE)

######


########################################
# (4) Dimension reduction
#####################################
# 1) pca
ref <- RunPCA(object = ref, pcs.compute = 50, do.print = FALSE)

# Examine and visualize PCA results
PCAPlot(object = ref, dim.1 = 1, dim.2 = 2)

PCHeatmap(object        = ref, 
          pc.use        = 1, 
          cells.use     = 500, 
          do.balanced   = TRUE, 
          label.columns = FALSE)

PCHeatmap(object        = ref, 
          pc.use        = 1:12, 
          cells.use     = 500, 
          do.balanced   = TRUE, 
          label.columns = FALSE)

# 2) select PCs
PCElbowPlot(object = ref, num.pc = 50)
PCElbowPlot(object = ref, num.pc = 30)

pcs.use <- 1:20

# 3) umap
ref <- RunUMAP(object = ref, dims.use = pcs.use, seed.use = 42)

#####


########################################
# (5) Cluster
#####################################
ref <- FindClusters(object       = ref, 
                    dims.use     = pcs.use, 
                    resolution   = seq(0.2,0.4,0.1), 
                    print.output = FALSE, 
                    save.SNN     = TRUE)

# investigating how many clusters each resolution produces 
sapply(X   = grep('^res',colnames(ref@meta.data),value = TRUE),
       FUN = function(x) length(unique(ref@meta.data[,x])))

# set resolutions
ref <- SetAllIdent(ref, id = 'res.0.3')

# save the object
saveRDS(ref, file = "./ref_cluster.rds")
# ref <- readRDS("./ref_cluster.rds")

#####


########################################
# (6) canonical markers
#####################################
# set resolutions
ref <- SetAllIdent(ref, id = 'res.0.3')

# pca visualization
DimPlot(object        = ref, 
        reduction.use = 'pca', 
        do.label      = TRUE, 
        coord.fixed   = TRUE, 
        group.by      = 'tissue_marker')

DimPlot(object        = ref, 
        reduction.use = 'pca', 
        do.label      = TRUE, 
        coord.fixed   = TRUE, 
        group.by      = 'ident')

# umap visualization
DimPlot(object        = ref, 
        reduction.use = 'umap', 
        do.label      = TRUE, 
        coord.fixed   = TRUE, 
        group.by      = 'tissue_marker')

DimPlot(object        = ref, 
        reduction.use = 'umap', 
        do.label      = TRUE, 
        coord.fixed   = TRUE, 
        group.by      = 'ident')

# QC index
FeaturePlot(object = ref, nCol = 2, reduction.use = 'umap', coord.fixed = TRUE,
            features.plot = c('nGene', 'nUMI', 'percent.mt', 'nUMI.per.Gene'))

# T 1
FeaturePlot(object = ref, nCol = 2, reduction.use = 'umap', coord.fixed = TRUE,
            features.plot = c('Cd3g', 'Cd3e', 'Bcl11b', 'Cd8a'))

# B 5,11
FeaturePlot(object = ref, nCol = 3, reduction.use = 'umap', coord.fixed = TRUE,
            features.plot = c('Cd79a', 'Cd79b', 'Fcmr', 'Ms4a1', 'Cd19', 'Cd22'))

# Mono 3,8
FeaturePlot(object = ref, nCol = 2, reduction.use = 'umap', coord.fixed = TRUE,
            features.plot = c('Fn1', 'Ccr2', 'Cx3cr1', 'S100a4'))

# DC 9
FeaturePlot(object = ref, nCol = 3, reduction.use = 'umap', coord.fixed = TRUE,
            features.plot = c('Siglech', 'Cd209d', 'Klk1', 'Ccr9', 'Cd209a'))

# platelet 7,4,2,0
FeaturePlot(object = ref, nCol = 2, reduction.use = 'umap', coord.fixed = TRUE,
            features.plot = c('Cd9','Pf4'))

# red cell
FeaturePlot(object = ref, nCol = 3, reduction.use = 'umap', coord.fixed = TRUE,
            features.plot = c('Hbb-bt', 'Hbb-bs', 'Hba-a1', 'Alas2', 'Snca'))

# hsc and progenitor 10,6
FeaturePlot(object = ref, nCol = 5, reduction.use = 'umap', coord.fixed = TRUE,
            features.plot = c('Mpl','Cpox','Gypa','Klf1','Ly6a','Gata2','Kit',
                              'Cd34','Arap2','Mpo','Elane'))

# Neu 7,4,2,0
FeaturePlot(object = ref, nCol = 3, reduction.use = 'umap', coord.fixed = TRUE,
            features.plot = c('Ly6g', 'S100a8', 'S100a9', 'Cd177', 'Mmp8'))

# granules 6,7,4,2,0
FeaturePlot(object = ref, nCol = 3, reduction.use = 'umap', coord.fixed = TRUE,
            features.plot=c('Elane', 'Ltf', 'Mmp8', 
                            'Prtn3', 'Camp', 'Mmp9'))

FeaturePlot(object = ref, nCol = 3, reduction.use = 'umap', coord.fixed = TRUE,
            features.plot=c('Ccr2', 'C1qc', 'Cx3cr1', 
                            'Ly6g', 'S100a8', 'Cxcr2'))

#####


########################################
# (7) rename cluster and visualization
#####################################
# 1) rename cluster
ref <- SetAllIdent(ref, id = 'res.0.3')

old.id <- c(0,2,4,7,6,10,9,3,8,5,11,1,12)
new.id <- c(rep('Neu',4),
            'MP', 
            'HSPC',
            'DC', 
            rep('Mono',2),
            rep('B',2), 
            'T', 
            'Cont')

# rename
ref@ident <- plyr::mapvalues(x    = ref@ident, 
                             from = old.id, 
                             to   = new.id)

# 2) save this to meta.data
ref@meta.data$cell_type <- ref@ident

head(ref@meta.data)
table(ref@meta.data$cell_type)

ref@meta.data$cell_type <- factor(x     = ref@meta.data$cell_type, 
                                  levels = unique(new.id))

# 3) visualization umap
DimPlot(object        = ref, 
        reduction.use = 'umap', 
        do.label      = FALSE, 
        coord.fixed   = TRUE,
        do.return     = TRUE,
        group.by      = 'tissue_marker',
        cols.use      =  brewer.pal(12, 'Paired')[c(1,2,4,8)],
        pt.size       = 0.1) -> p1

DimPlot(object        = ref, 
        reduction.use = 'umap', 
        do.label      = TRUE, 
        coord.fixed   = TRUE,
        do.return     = TRUE,
        group.by      = 'cell_type',
        cols.use      = brewer.pal(8, 'Dark2')[1:8],
        pt.size       = 0.1) -> p2

library(ggpubr)
ggarrange(p1,p2, ncol = 2, align  = 'hv')

ggsave(filename = "./ref_umap_tissue.pdf", plot = p1, width = 7, height = 7)
ggsave(filename = "./ref_umap_celltype.pdf", plot = p2, width = 7, height = 7)

ggsave(filename = "./ref_umap_tissue.eps", plot = p1, width = 7, height = 7)
ggsave(filename = "./ref_umap_celltype.eps", plot = p2, width = 7, height = 7)

# save the object
head(ref@meta.data)
ref@meta.data$res.0.2 <- NULL
ref@meta.data$res.0.3 <- NULL
ref@meta.data$res.0.4 <- NULL
ref@scale.data <- NULL

saveRDS(ref, file = "./ref_rename.rds")
# ref <- readRDS("./ref_rename.rds")

# save metadata
write.table(ref@meta.data, "./ref_meta.txt", quote = FALSE)

#####


########################################
# (8) DEG and visualization
#####################################
# 1) filter contaminations
ref      <- SetAllIdent(ref, id = 'cell_type')
ref.pure <- SubsetData(object = ref, ident.remove = 'Cont')

table(ref.pure@ident)
ref.pure@ident <- factor(x      = ref.pure@ident, 
                         levels = c('Neu', 'MP', 'HSPC', 'DC', 'Mono', 'B', 'T'))

# 2) ave exp
cell.type.ave.exp <- AverageExpression(ref.pure)
write.csv(cell.type.ave.exp, "./cell_type_ave_exp.csv", quote = F)
# cell.type.ave.exp <- read.csv("./cell_type_ave_exp.csv", row.names = 1)

# 3) find markers
cell.type.markers <- FindAllMarkers(object          = ref.pure, 
                                    logfc.threshold = log(2),
                                    test.use        = 't',
                                    only.pos        = TRUE,
                                    do.print        = FALSE)

cell.type.markers <- cell.type.markers[cell.type.markers$p_val_adj < 0.05, ]
table(cell.type.markers$cluster)

cell.type.markers.top30 <- cell.type.markers %>% group_by(cluster) %>% top_n(30, avg_logFC)
table(cell.type.markers.top30$cluster)

write.csv(cell.type.markers, "./cell_type_markers_log2.csv", quote = F)
# cell.type.markers <- read.csv("./cell_type_markers.csv", row.names = 1, stringsAsFactors = FALSE)
write.csv(cell.type.markers.top30, "./cell_type_markers_top30.csv", quote = F)

# 4) preapare gene
cell.type.markers %>% group_by(cluster) %>% top_n(5, avg_logFC) -> topgene

gene.plot <- unique(topgene$gene)
gene.plot <- c(gene.plot[1:13],'Atpif1','Hsp90aa1',gene.plot[16:35])

# 5) heatmap using ComplexHeatmap
library(ComplexHeatmap)

order.cell <- rownames(ref.pure@meta.data[order(ref.pure@meta.data$cell_type),])

# annotation column
ann.col = data.frame(Cell_type = ref.pure@meta.data[order.cell, 'cell_type'])
ann.col = HeatmapAnnotation(df  = ann.col,
                            col = list(Cell_type = c(Neu=brewer.pal(8, 'Dark2')[1], 
                                                     MP=brewer.pal(8, 'Dark2')[2], 
                                                     HSPC=brewer.pal(8, 'Dark2')[3], 
                                                     DC=brewer.pal(8, 'Dark2')[4], 
                                                     Mono=brewer.pal(8, 'Dark2')[5],
                                                     B=brewer.pal(8, 'Dark2')[6],
                                                     T=brewer.pal(8, 'Dark2')[7])))

# plot
library(circlize)

pdf("./ref_heatmap_ComplexHeatmap.pdf")
Heatmap(matrix            = t(scale(t(ref.pure@data[gene.plot, order.cell]))),
        col = colorRamp2(seq(-2,2,length.out = 7),rev(brewer.pal(n = 7, name = "RdBu"))),
        name              = 'Exp',
        width             = unit(10,'cm'),
        height            = unit(12,'cm'),
        cluster_rows      = FALSE,
        cluster_columns   = FALSE,
        show_row_names    = TRUE,
        show_column_names = FALSE,
        top_annotation    = ann.col)
dev.off()

#####


########################################
# (9) extract neutrophils
#####################################
# 1) including bm_ckit, set up reference
ref <- SetAllIdent(ref, id = 'cell_type')

ref.neu <- SubsetData(object     = ref, 
                      ident.use  = c('MP', 'Neu'),
                      do.clean   = TRUE,
                      subset.raw = TRUE)

ref.neu@raw.data@Dim
# 27998 12285

table(ref.neu@meta.data$tissue_marker)
table(ref.neu@meta.data$cell_type)

ref.neu@meta.data$cell_type <- factor(ref.neu@meta.data$cell_type)

# 2) excluding bm_ckit
ref.neu    <- SetAllIdent(ref.neu, id = 'tissue')
wt.ctl.neu <- SubsetData(object       = ref.neu, 
                         ident.remove = 'WT_BM_C-kit+Gr1',
                         do.clean     = TRUE,
                         subset.raw   = TRUE)

wt.ctl.neu@raw.data@Dim
# 27998 10127

table(wt.ctl.neu@meta.data$tissue)
table(wt.ctl.neu@meta.data$cell_type)

wt.ctl.neu@meta.data$tissue <- factor(wt.ctl.neu@meta.data$tissue)

# save the object
saveRDS(ref.neu, file = "./ref_neu_raw.rds")
saveRDS(wt.ctl.neu, file = "./wt_ctl_neu_raw.rds")

#####

