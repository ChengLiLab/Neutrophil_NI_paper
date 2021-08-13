# Sample: neu selected from wt.ctl.bm1, wt.ctl.bm2, wt.ctl.pb2 and wt.ctl.sp2
# Aim: establish reference by basic analysis
#
# 2020-03-20, Qiang Shi

rm(list = ls())
setwd("/lustre/user/liclab/shiq/single_cell/luolab/neu_profile/3_wt_ctl/2_ref_neu")

library(Seurat, lib.loc = "/lustre/user/liclab/shiq/software/R_library/")
library(dplyr)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
library(RColorBrewer)
library(ggsci)

ann.color <- list(tissue = c('WT_BM_C-kit+Gr1' = brewer.pal(12, "Paired")[1],
                             WT_BM_Gr1 = brewer.pal(12, "Paired")[2],
                             WT_PB_Gr1 = brewer.pal(12, "Paired")[4],
                             WT_SP_Gr1 = brewer.pal(12, "Paired")[8]),
                  cluster = c(G0 = pal_npg()(10)[1], 
                              G1 = pal_npg()(10)[2], 
                              G2 = pal_npg()(10)[3], 
                              G3 = pal_npg()(10)[4], 
                              G4 = pal_npg()(10)[5],
                              G5a = pal_npg()(10)[6], 
                              G5b = pal_npg()(10)[7], 
                              G5c = pal_npg()(10)[8],
                              GM = pal_npg()(10)[9],
                              Cont = pal_npg()(10)[10]))

########################################
# (1) Setup the Seurat Object
#####################################
ref.neu <- readRDS(file = "../1_ref/ref_neu_raw.rds")

ref.neu@raw.data@Dim
# 27998 12285

ref.neu@data@Dim
# 27998 12285

head(ref.neu@meta.data)

table(ref.neu@meta.data$tissue_marker)
# WT_BM_C-kit+Gr1       WT_BM_Gr1       WT_PB_Gr1       WT_SP_Gr1 
# 2158                  3591            4955            1581

#####


########################################
# (2) QC
#####################################
pdf("./ref_neu_qc.pdf")
VlnPlot(object         = ref.neu, 
        features.plot  = c('nGene', 'nUMI', 'percent.mt', 'nUMI.per.Gene'),
        nCol           = 2,
        size.x.use     = 10,
        point.size.use = 0,
        cols.use       = brewer.pal(12, "Paired")[c(1,2,4,8)],
        group.by       = 'tissue_marker',
        x.lab.rot      = TRUE)
dev.off()

# GenePlot 
GenePlot(object = ref.neu, gene1 = 'nUMI', gene2 = 'percent.mt')
GenePlot(object = ref.neu, gene1 = 'nUMI', gene2 = 'nGene')

#####


########################################
# (3) Normalization, Feature, Scale
#####################################
# normalization
ref.neu <- NormalizeData(ref.neu)

# variable genes
ref.neu <- FindVariableGenes(object        = ref.neu, 
                             x.low.cutoff  = 0.0125, 
                             x.high.cutoff = 4, 
                             y.cutoff      = 1.2)

length(x = ref.neu@var.genes)
# 2392

write.table(ref.neu@var.genes, "./ref_neu_vargene.txt", row.names = FALSE)

# scale
head(ref.neu@meta.data)
# regress c("tissue", "percent.mt") has the best clusters, each tissue corresponds to each dataset
ref.neu <- ScaleData(object           = ref.neu, 
                     vars.to.regress  = c('tissue','percent.mt'), #c('batch',"nUMI",),
                     display.progress = F, 
                     do.par           = TRUE,
                     num.cores        = 24)

######


########################################
# (4) Dimention reduction
#####################################
# 1) pca
ref.neu <- RunPCA(object = ref.neu, pcs.compute = 50, do.print = FALSE)

# Examine and visualize PCA results
PCAPlot(object = ref.neu, dim.1 = 1, dim.2 = 2)

PCHeatmap(object        = ref.neu, 
          pc.use        = 1, 
          cells.use     = 500, 
          do.balanced   = TRUE, 
          label.columns = FALSE)

PCHeatmap(object        = ref.neu, 
          pc.use        = 1:12, 
          cells.use     = 500, 
          do.balanced   = TRUE, 
          label.columns = FALSE)

# 2) select PCs
PCElbowPlot(object = ref.neu, num.pc = 50)
PCElbowPlot(object = ref.neu, num.pc = 30)

pcs.use <- 1:15

# 3) umap
ref.neu <- RunUMAP(object = ref.neu, dims.use = pcs.use, seed.use = 42)

#####


########################################
# (5) Cluster
#####################################
ref.neu <- FindClusters(object       = ref.neu, 
                        dims.use     = pcs.use, 
                        resolution   = seq(0.5,0.7,0.1),
                        print.output = FALSE,
                        save.SNN     = TRUE)

# investigating how many clusters each resolution produces 
sapply(X   = grep('^res', colnames(ref.neu@meta.data), value = TRUE),
       FUN = function(x) length(unique(ref.neu@meta.data[,x])))

# set resolutions
ref.neu <- SetAllIdent(ref.neu, id = 'res.0.6')

# save the object
saveRDS(ref.neu, file = "./ref_neu_cluster.rds")
# ref.neu <- readRDS("./ref_neu_cluster.rds")

#####


########################################
# (6) canonical markers
#####################################
# set resolutions
ref.neu <- SetAllIdent(ref.neu, id = 'res.0.6')

# pca visualization
DimPlot(object        = ref.neu, 
        reduction.use = 'pca', 
        do.label      = TRUE, 
        coord.fixed   = TRUE, 
        group.by      = 'tissue')

DimPlot(object        = ref.neu, 
        reduction.use = 'pca', 
        do.label      = TRUE, 
        coord.fixed   = TRUE, 
        group.by      = 'ident')

# umap visualization
DimPlot(object        = ref.neu, 
        reduction.use = 'umap', 
        do.label      = TRUE, 
        coord.fixed   = TRUE, 
        group.by      = 'tissue')

DimPlot(object        = ref.neu, 
        reduction.use = 'umap', 
        do.label      = TRUE, 
        coord.fixed   = TRUE, 
        group.by      = 'ident')

# QC index
FeaturePlot(object = ref.neu, nCol = 2, reduction.use = 'umap', coord.fixed = TRUE,
            features.plot = c('nGene', 'nUMI', 'percent.mt', 'nUMI.per.Gene'))

# T
FeaturePlot(object = ref.neu, nCol = 2, reduction.use = 'umap', coord.fixed = TRUE,
            features.plot = c('Cd3g', 'Cd3e', 'Bcl11b', 'Cd8a'))

# B
FeaturePlot(object = ref.neu, nCol = 3, reduction.use = 'umap', coord.fixed = TRUE,
            features.plot = c('Cd79a', 'Cd79b', 'Fcmr', 'Ms4a1', 'Cd19', 'Cd22'))

# Mono
FeaturePlot(object = ref.neu, nCol = 2, reduction.use = 'umap', coord.fixed = TRUE,
            features.plot = c('Fn1', 'Ccr2', 'Cx3cr1', 'S100a4'))

# DC
FeaturePlot(object = ref.neu, nCol = 3, reduction.use = 'umap', coord.fixed = TRUE,
            features.plot = c('Siglech', 'Cd209d', 'Klk1', 'Ccr9', 'Cd209a'))

# platelet 5,1,8,2,0,3,4,9
FeaturePlot(object = ref.neu, nCol = 2, reduction.use = 'umap', coord.fixed = TRUE,
            features.plot = c('Cd9','Pf4'))

# red cell
FeaturePlot(object = ref.neu, nCol = 3, reduction.use = 'umap', coord.fixed = TRUE,
            features.plot = c('Hbb-bt', 'Hbb-bs', 'Hba-a1', 'Alas2', 'Snca'))

# hsc and progenitor 6,7
FeaturePlot(object = ref.neu, nCol = 5, reduction.use = 'umap', coord.fixed = TRUE,
            features.plot = c('Mpl','Cpox','Gypa','Klf1','Ly6a','Gata2','Kit',
                              'Cd34','Arap2','Mpo','Elane'))

# Neu 5,1,8,2,0,3,4,9
FeaturePlot(object = ref.neu, nCol = 3, reduction.use = 'umap', coord.fixed = TRUE,
            features.plot = c('Ly6g', 'S100a8', 'S100a9', 'Cd177', 'Mmp8'))

# granules 7,5,1,8,2,0,3,4,9
FeaturePlot(object = ref.neu, nCol = 3, reduction.use = 'umap', coord.fixed = TRUE,
            features.plot=c('Elane', 'Ltf', 'Mmp8', 
                            'Prtn3', 'Camp', 'Mmp9'))

# find markers before rename
neu.cluster.deg <- FindAllMarkers(object          = ref.neu, 
                                  logfc.threshold = log(2),
                                  test.use        = 't',
                                  only.pos        = FALSE,
                                  do.print        = FALSE)

write.csv(neu.cluster.deg, "./neu_cluster_deg.csv", quote = F)

#####


########################################
# (7) rename cluster and visualization
#####################################
# 1) rename
ref.neu <- SetAllIdent(ref.neu, id = 'res.0.6')

old.id  <- c(6,7,5,1,2,0,4,3,8,9)

new.id <- c('G0','G1','G2','G3','G4','G5a','G5b','G5c','GM','Cont')

# rename
ref.neu@ident <- plyr::mapvalues(x    = ref.neu@ident, 
                                 from = old.id, 
                                 to   = new.id)

# stash this into meta.data
ref.neu@meta.data$cluster <- ref.neu@ident

table(ref.neu@meta.data$cluster)

ref.neu@meta.data$cluster <- factor(x      = ref.neu@meta.data$cluster, 
                                    levels = new.id)

# 2) visualization
DimPlot(object        = ref.neu, 
        reduction.use = 'umap', 
        do.label      = FALSE, 
        coord.fixed   = TRUE,
        group.by      = 'tissue_marker',
        cols.use      = ann.color$tissue,
        pt.size       = 0.5) -> p1

DimPlot(object        = ref.neu, 
        reduction.use = 'umap', 
        do.label      = TRUE, 
        coord.fixed   = TRUE,
        group.by      = 'cluster',
        cols.use      = ann.color$cluster,
        pt.size       = 0.5) -> p2

DimPlot(object        = ref.neu, 
        reduction.use = 'umap', 
        do.label      = TRUE,
        coord.fixed   = TRUE,
        group.by      = 'cell_type',
        cols.use      = brewer.pal(8, 'Dark2')[c(1:2)],
        pt.size       = 0.5) -> p3

ggsave("./ref_neu_umap_tissue.pdf", plot = p1, width = 7, height = 7)
ggsave("./ref_neu_umap_cluster.pdf", plot = p2, width = 7, height = 7)
ggsave("./ref_neu_umap_celltype.pdf", plot = p3, width = 7, height = 7)

ggsave("./ref_neu_umap_tissue.eps", plot = p1, width = 7, height = 7)
ggsave("./ref_neu_umap_cluster.eps", plot = p2, width = 7, height = 7)
ggsave("./ref_neu_umap_celltype.eps", plot = p3, width = 7, height = 7)

# save the object
head(ref.neu@meta.data)
ref.neu@meta.data$res.0.5 <- NULL
ref.neu@meta.data$res.0.6 <- NULL
ref.neu@meta.data$res.0.7 <- NULL
ref.neu@scale.data <- NULL

saveRDS(ref.neu, file = "./ref_neu_rename.rds")
# ref.neu <- readRDS("./ref_neu_rename.rds")

# save metadata
write.table(ref.neu@meta.data, "./ref_neu_meta.txt", quote = FALSE)

meta <- ref.neu@meta.data
meta <- cbind(meta, ref.neu@dr$umap@cell.embeddings)
write.table(meta, "./ref_neu_meta2.txt", quote = FALSE)

#####


########################################
# (8) find DEF and visualization
#####################################
# 1) filter contaminations
ref.neu      <- SetAllIdent(ref.neu, id = 'cluster')
ref.neu.pure <- SubsetData(object = ref.neu, ident.remove = 'Cont')

table(ref.neu.pure@ident)

# 2) ave exp
neu.cluster.ave.exp <- AverageExpression(ref.neu.pure)
write.csv(neu.cluster.ave.exp, "./neu_cluster_ave_exp.csv", quote = F)
# neu.cluster.ave.exp <- read.csv("./neu_cluster_ave_exp.csv", row.names = 1)

# 3) find markers after rename
neu.cluster.deg <- FindAllMarkers(object          = ref.neu.pure, 
                                  logfc.threshold = log(1.2),
                                  test.use        = 't',
                                  only.pos        = TRUE)

neu.cluster.deg <- subset(neu.cluster.deg, p_val_adj<=0.05)

table(neu.cluster.deg$cluster)
#   G0   G1   G2   G3   G4  G5a  G5b  G5c   GM 
# 2787 2222  826  183   93  478  367  702   86 

nrow(neu.cluster.deg) # 7744

write.csv(neu.cluster.deg, "./neu_cluster_deg_log1.2.csv", quote = F)
# neu.cluster.deg <- read.csv("./neu_cluster_deg1.2.csv", row.names = 1, stringsAsFactors = FALSE)

# 4) preapare gene
neu.cluster.markers <- subset(neu.cluster.deg, avg_logFC>log(1.5) & cluster!='GM')
table(neu.cluster.markers$cluster)
#   G0   G1   G2   G3   G4  G5a  G5b  G5c   GM 
# 1109  932  257   53   18  110  132  264    0 

neu.cluster.markers %>% group_by(cluster) %>% top_n(20, avg_logFC)-> topgene

gene.plot <- unique(topgene$gene)
gene.plot <- c('Cd34', gene.plot[1:20], 'Kit', gene.plot[21:139])

# 5) heatmap using ComplexHeatmap
library(ComplexHeatmap)

meta <- ref.neu.pure@meta.data
meta <- subset(meta, cluster!='GM')

order.cell <- rownames(meta[order(meta$cluster),])

# annotation column
ann.col = data.frame(Cluster = meta[order.cell, 'cluster'])
ann.col = HeatmapAnnotation(df  = ann.col,
                            col = list(Cluster = c(G0=pal_npg()(10)[1], 
                                                   G1=pal_npg()(10)[2], 
                                                   G2=pal_npg()(10)[3], 
                                                   G3=pal_npg()(10)[4], 
                                                   G4=pal_npg()(10)[5],
                                                   G5a=pal_npg()(10)[6], 
                                                   G5b=pal_npg()(10)[7], 
                                                   G5c=pal_npg()(10)[8])))

# annotation row
gene.mark <- c('Cd34','Rpl12','Npm1','Kit','Elane','Mpo','Ctsg','Chil3','Camp','Fcnb','Ifitm6',
               'Ltf',
               'Ngp','Gm5483','Ccl6','Mmp8','Isg15','Rsad2','Ifit1','Il1b','Gm2a','Cxcr4','Wfdc17')

ann.row = rowAnnotation(foo = anno_mark(at = match(gene.mark, gene.plot), labels = gene.mark))

# plot
library(circlize)

pdf("./ref_neu_heatmap.pdf")
Heatmap(matrix = t(scale(t(ref.neu.pure@data[gene.plot, order.cell]))),
        col    = colorRamp2(seq(-2,2,length.out = 7),rev(brewer.pal(n = 7, name = "RdBu"))),
        name              = 'Exp',
        width             = unit(10,'cm'),
        height            = unit(12,'cm'),
        cluster_rows      = FALSE,
        cluster_columns   = FALSE,
        show_row_names    = FALSE,
        show_column_names = FALSE,
        top_annotation    = ann.col,
        right_annotation  = ann.row)
dev.off()

#####

