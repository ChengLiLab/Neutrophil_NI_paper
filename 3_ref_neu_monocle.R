# monocle for ref.neu
#
# 2020-03-02, Qiang Shi

rm(list = ls())
setwd("/lustre/user/liclab/shiq/single_cell/luolab/neu_profile/wt_ctl/ref_neu/analysis/monocle/")

library(monocle)
library(Seurat, lib.loc = "/lustre/user/liclab/shiq/software/R_library/")
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
library(RColorBrewer)
library(ggsci)

#######################################
# (0) loading data
####################################
# 1) load data
ref.neu <- readRDS( file = "../../ref_neu_rename.rds")

ref.neu@data@Dim
# 27998 12285

# 2) filter cont
head(ref.neu@meta.data)
ref.neu <- SetAllIdent(ref.neu, id = 'cluster')

table(ref.neu@meta.data$cluster)
#  G0   G1   G2   G3   G4  G5a  G5b  G5c   GM Cont 
# 509  436  699 2496 1837 3545  823 1647  163  130 

ref.neu <- SubsetData(object       = ref.neu, 
                      ident.remove = c('GM','Cont'),
                      do.clean     = FALSE,
                      subset.raw   = TRUE)

ref.neu@data@Dim
# 27998 11992

ref.neu@meta.data$cluster <- factor(ref.neu@meta.data$cluster)

#####


#######################################
# (1) set cds
####################################
library(reshape2)

pd <- new('AnnotatedDataFrame', data = ref.neu@meta.data)
fd <- rownames(ref.neu@data)
names(fd) <- fd
fd <- new('AnnotatedDataFrame', data = data.frame(gene_short_name = fd))

cds <- newCellDataSet(cellData = ref.neu@raw.data,
                       phenoData   = pd,
                       featureData = fd,
                       lowerDetectionLimit = 0.5,
                       expressionFamily = negbinomial.size())

head(pData(cds))
head(fData(cds))

cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)

#####


#######################################
# (2) construct trajectory
####################################
# 1) gene selection

# a) monocle
# disp_table <- dispersionTable(cds)
# head(disp_table)
# 
# unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.05 & dispersion_empirical >= 1.5 * dispersion_fit)
# cds <- setOrderingFilter(cds, unsup_clustering_genes$gene_id)

# b) var.gene by seurat, has similar results with monocle
ref.neu.var.gene <- read.table("../../ref_neu_vargene.txt")
cds <- setOrderingFilter(cds, ref.neu.var.gene$V1)

# c) combination of top 100 degs of each cluster
# neu.cluster.markers <- read.csv(file      = "../../neu_cluster_deg.csv", 
#                                 row.names = 1, stringsAsFactors = F)
# 
# neu.cluster.markers <- neu.cluster.markers[neu.cluster.markers$avg_logFC>0,]
# 
# table(neu.cluster.markers$cluster)
# 
# library(dplyr)
# neu.cluster.markers <- neu.cluster.markers %>% group_by(cluster) %>% top_n(100, avg_logFC)
# cds <- setOrderingFilter(cds, unique(neu.cluster.markers$gene))

plot_ordering_genes(cds)

# 2) dimensionality reduction
cds <- reduceDimension(cds)

# 3) order cells along the trajectory
set.seed(7)
cds <- orderCells(cds, reverse = TRUE)

head(pData(cds))
head(fData(cds))

# 4) plot_cell_trajectory
colScale1 <- scale_colour_manual(values = brewer.pal(12, 'Paired')[c(1,2,4,8)])
colScale2 <- scale_colour_manual(values = pal_npg()(10)[1:8])

p1 <- plot_cell_trajectory(cds, color_by = 'tissue_marker', cell_size = 1, show_branch_points = FALSE) + colScale1
p2 <- plot_cell_trajectory(cds, color_by = 'cluster', cell_size = 1, show_branch_points = FALSE) + colScale2
p3 <- plot_cell_trajectory(cds, color_by = 'Pseudotime', cell_size = 1, show_branch_points = FALSE)

ggsave("./monocle_cell_trajectory_tissue.pdf", plot = p1, width = 7, height = 7)
ggsave("./monocle_cell_trajectory_cluster.pdf", plot = p2, width = 7, height = 7)
ggsave("./monocle_cell_trajectory_pseudotime.pdf", plot = p3, width = 7, height = 7)

ggsave("./monocle_cell_trajectory_tissue.eps", plot = p1, width = 7, height = 7)
ggsave("./monocle_cell_trajectory_cluster.eps", plot = p2, width = 7, height = 7)
ggsave("./monocle_cell_trajectory_pseudotime.eps", plot = p3, width = 7, height = 7)

# save rds
saveRDS(cds, "./monocle_cds.rds")
# cds <- readRDS("./monocle_cds.rds")

#####

