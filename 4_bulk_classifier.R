# deconvolution from jiayu
#
# 2020-01-15, Qiang Shi

rm(list = ls())
setwd("/lustre/user/liclab/shiq/single_cell/luolab/neu_profile/4_bulk/")

library(ComplexHeatmap)
library(ggplot2)
library(ggsci)
library(RColorBrewer)
library(edgeR)
library(Seurat)
library(dplyr)

######################################################################
## (1) preprocess data
###################################################################

files <- dir(path="./bulk_data/", pattern="*\\.tab$")
files <- paste("./bulk_data/", files, sep = "")
bulk <- readDGE(files, columns=c(1,3), header=FALSE)

# Read in 10X single-cell feature annotations
feature <- read.table("./features.tsv.gz", row.names = 1)
gene.intersect <- intersect(rownames(bulk), rownames(feature))
bulk <- bulk[gene.intersect,]
rownames(bulk) <- feature[rownames(bulk),]$V2

dim(bulk)
# 30309    15

rm(files, feature)

#####


######################################################################
## (2) bulk process
###################################################################

FQnorm <- function(counts){
  # Quantile normalization function
  rk <- apply(counts,2,rank,ties.method='min')
  counts.sort <- apply(counts,2,sort)
  refdist <- apply(counts.sort,1,median)
  norm <- apply(rk,2,function(r){ refdist[r] })
  rownames(norm) <- rownames(counts)
  return(norm)
}

# Quantile normalization
bulk <- FQnorm(bulk$counts)

# Combine bulk replicates
bulk <- cbind(rowMeans(bulk[,1:3]),
                    rowMeans(bulk[,4:6]),
                    rowMeans(bulk[,7:9]),
                    rowMeans(bulk[,10:12]),
                    rowMeans(bulk[,13:15]))

colnames(bulk)
colnames(bulk) <- c("Neu", "MB", "MC", "MM", "PM")
bulk <- bulk[,c(2,5,3,4,1)]

sc.meta = read.table("../wt_ctl/ref_neu/ref_neu_meta.txt")
str(sc.meta)

sc.meta = sc.meta[sc.meta$cluster %in% c("G0", "G1", "G2", "G3","G4"),]

sc.meta$cluster <- factor(sc.meta$cluster)

#####


######################################################################
## (3) Compare bulk and single-cell expression profiles
###################################################################

# Load single-cell defined signatures in bulk
deg = read.csv("../wt_ctl/ref_neu/neu_cluster_deg_log1.2.csv", row.names = 1)
library(dplyr)
top20 <- deg %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)
genes.used <- as.character(top20$gene)[top20$gene %in% rownames(bulk)[rowSums(bulk)>0]]

# Load merged single-cell expression data
ave.exp = read.csv("../wt_ctl/ref_neu/neu_cluster_ave_exp.csv", row.names = 1)
ave.exp = ave.exp[,c("G0","G1","G2","G3","G4")]
ave.exp <- ave.exp[genes.used,]
ave.exp.scaled = t(scale(t(ave.exp)))

#####


######################################################################
## (4) Bulk deconvolution by nnls regression,
###################################################################
library(nnls)
# Scale each gene to have the same mean and variance
bulk.scaled = t(scale(t(bulk[genes.used,])))
deconv.result = c()
# Apply nnls regression
for (i in 1:(dim(bulk.scaled)[2])){
  out = nnls(ave.exp.scaled, bulk.scaled[,i])
  deconv.result = rbind(deconv.result, out$x)
}

# deconv.result
cmatrix = rbind(c(1.047437,0.6470013,0,0.4004393,0.3953729),
                c(0.3011903,0.8216298,0.1412882,0,0.4173733),
                c(0.1861662,0.3115887,0.7696528,0.1423949,0),
                c(0.1994302,0,0.8283148,0.491128,0.3154364),
                c(0.04599633,0,0.04096413,0.7462578,0.6520374))

cmatrix <- deconv.result
rownames(cmatrix) = c("MB", "PM", "MC", "MM", "Neu")
colnames(cmatrix) = c("G0", "G1", "G2", "G3", "G4")

# Scale the composition matrix by group sizes
# scale.prop = matrix(table(sc.meta$cluster)/dim(sc.meta)[1])[,1]
# cmatrix = t(cmatrix)/scale.prop
# c = t(cmatrix)/colSums(cmatrix)
# c <- t(c)

c <- t(cmatrix)
c = c/colSums(c)

# Visualization
hm1 = Heatmap(c, col = brewer.pal(n=9, "OrRd"), name = "Proportion", 
              cluster_rows = F, cluster_columns = F,
              width = unit(2.65, "cm"), height = unit(6.5, "cm"),
              row_names_side = "left",
              heatmap_legend_param = list(title_gp = gpar(fontsize=7),
                                          labels_gp = gpar(fontsize=7),
                                          grid_height = unit(0.2, "cm"),
                                          legend_direction = "horizontal"),
              column_names_rot = 90,
              row_names_gp = gpar(fontsize=7),
              column_names_gp = gpar(fontsize=7))
draw(hm1, heatmap_legend_side="top")
decorate_heatmap_body("Proportion", {
  grid.rect(gp = gpar(fill = "transparent", col = "black", lwd = 0.5))
})

#####

