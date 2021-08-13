# aging scores
# 
# 2020-03-25, Qiang Shi

rm(list = ls())
setwd("/lustre/user/liclab/shiq/single_cell/luolab/neu_profile/3_wt_ctl/4_analysis/5_score")
options(stringsAsFactors = FALSE)

library(Seurat, lib.loc = "/lustre/user/liclab/shiq/software/R_library/")
library(dplyr)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
library(pheatmap)
library(ComplexHeatmap)
library(RColorBrewer)
library(ggsci)
library(ggpubr)

ann.color <- list(cluster = c(G0=pal_npg()(10)[1], 
                              G1=pal_npg()(10)[2], 
                              G2=pal_npg()(10)[3], 
                              G3=pal_npg()(10)[4],
                              G4=pal_npg()(10)[5], 
                              G5a=pal_npg()(10)[6], 
                              G5b=pal_npg()(10)[7], 
                              G5c=pal_npg()(10)[8]),
                  tissue = c(WT_BM=brewer.pal(12, 'Paired')[2],
                             WT_PB=brewer.pal(12, 'Paired')[4],
                             WT_SP=brewer.pal(12, 'Paired')[8]))

#######################################
# (0) loading g5 rds data
####################################
ref.neu <- readRDS( file = "../../2_ref_neu/ref_neu_rename.rds")
ref.neu@data@Dim
# 27998 12285
head(ref.neu@meta.data)

# filter cont
table(ref.neu@meta.data$cluster)
#  G0   G1   G2   G3   G4  G5a  G5b  G5c   GM Cont 
# 509  436  699 2496 1837 3545  823 1647  163  130 

ref.neu <- SetAllIdent(ref.neu, id = 'cluster')
g5 <- SubsetData(object = ref.neu, ident.use = c('G5a','G5b','G5c'))

g5@data@Dim
# 27998 6015

g5@meta.data$cluster <- factor(g5@meta.data$cluster)
table(g5@meta.data$cluster)
#  G5a  G5b  G5c 
# 3545  823 1647 

#####


#######################################
# (1-1) cluster aging score
####################################
gene <- c("Sell","Itgam","Itga4","Cxcr4","Cxcr2",
          "Cd47","Cd24a","Tlr4","Icam1","Itgax")
weight <- c(-1,1,1,1,-1,-1,1,1,1,1)/10

z_matrix = t(scale(t(g5@data[gene,])))

aging.score = t(z_matrix[gene,]) %*% weight 
g5@meta.data$aging.score = aging.score[,1]

meta <- g5@meta.data

# Visualization
ggplot(meta, aes(x=cluster,y=aging.score, fill=cluster)) + 
  geom_violin() + geom_boxplot(width=0.2,outlier.shape = NA) + 
  labs(title = 'All tissue', y = "Aging Score") +
  stat_compare_means(comparisons = list(c("G5a", "G5b"), c("G5a", "G5c"),
                                        c("G5b", "G5c")), 
                     label = "p.signif", #label.y = c(1.7,1.9,2.1),
                     method = "t.test")+
  scale_fill_manual(values = ann.color$cluster)+
  theme_minimal() + 
  theme(
    axis.text = element_text(size=7, color = "black"),
    plot.title = element_text(hjust = 0.5),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size=7,color = "Black"),
    panel.grid=element_blank(),
    axis.line = element_line(),
    legend.position = "none") -> ps1

#####


#######################################
# (1-2) cluster aging proportion
####################################
# Define aged subgroups by Gaussian mixture model
require(mixtools)
out = normalmixEM(aging.score, k=2)
plot(out,2) # check the order
aging.group <- apply(out$posterior, 1, which.max)
meta$aging.group <- aging.group

# Check the aged group's composition
table(meta$aging.group, meta$cluster)
tapply(meta$aging.score, meta$cluster, summary)

# proportion
prop <- t(table(meta$aging.group, meta$cluster)) / 
  colSums(table(meta$aging.group, meta$cluster))
prop

df <- data.frame(cluster = c("G5a","G5b","G5c"), prop = prop[,2])

pp = ggplot(df, aes(x=cluster, y=prop, fill=cluster)) +
  geom_bar(width = 0.75, stat = "identity", color="black") +
  labs(title = 'All tissue', y = "Proportion of aged cells")+
  scale_fill_manual(values=ann.color$cluster) +
  theme_minimal() +
  theme(
    axis.text = element_text(size=7, color = "black"),
    plot.title = element_text(hjust = 0.5),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size=7,color = "Black"),
    panel.grid=element_blank(),
    axis.line = element_line(),
    legend.position = "none") -> pp1

ggarrange(ps1,pp1,nrow = 1, align = "hv")

#####


#######################################
# (2-1) tissue aging score
####################################
g5 <- SetAllIdent(g5, id = 'tissue')
g5 <- SubsetData(object = g5, ident.remove = 'WT_BM_C-kit+Gr1')
g5@meta.data$tissue <- factor(g5@meta.data$tissue)

gene <- c("Sell","Itgam","Itga4","Cxcr4","Cxcr2",
          "Cd47","Cd24a","Tlr4","Icam1","Itgax")
weight <- c(-1,1,1,1,-1,-1,1,1,1,1)/10

z_matrix = t(scale(t(g5@data[gene,])))

aging.score = t(z_matrix[gene,]) %*% weight 
g5@meta.data$aging.score = aging.score[,1]

meta <- g5@meta.data

# Visualization
ggplot(meta, aes(x=tissue,y=aging.score, fill=tissue)) + 
  geom_violin() + geom_boxplot(width=0.2,outlier.shape = NA) + 
  labs(title = 'G5', y = "Aging Score") +
  stat_compare_means(comparisons = list(c("WT_BM", "WT_PB"), 
                                        c("WT_BM", "WT_SP"),
                                        c("WT_PB", "WT_SP")), 
                     label = "p.signif", #label.y = c(1.7,1.9,2.1),
                     method = "t.test")+
  scale_fill_manual(values = ann.color$tissue)+
  theme_minimal() + 
  theme(
    axis.text = element_text(size=7, color = "black"),
    plot.title = element_text(hjust = 0.5),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size=7,color = "Black"),
    panel.grid=element_blank(),
    axis.line = element_line(),
    legend.position = "none") -> ps2

#####


#######################################
# (2-2) tissue aging proportion
####################################
# Define aged subgroups by Gaussian mixture model
require(mixtools)
out = normalmixEM(aging.score, k=2)
plot(out,2) # check the order
aging.group <- apply(out$posterior, 1, which.max)
meta$aging.group <- aging.group

# Check the aged group's composition
table(meta$aging.group, meta$tissue)
tapply(meta$aging.score, meta$tissue, summary)

# proportion
prop <- t(table(meta$aging.group, meta$tissue)) / 
  colSums(table(meta$aging.group, meta$tissue))
prop

df <- data.frame(tissue = c("WT_BM","WT_PB","WT_SP"), prop = prop[,2])

ggplot(df, aes(x=tissue, y=prop, fill=tissue)) +
  geom_bar(width = 0.75, stat = "identity", color="black") +
  labs(title = 'G5', y = "Proportion of aged cells")+
  scale_fill_manual(values=ann.color$tissue) +
  theme_minimal() +
  theme(
    axis.text = element_text(size=7, color = "black"),
    plot.title = element_text(hjust = 0.5),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size=7,color = "Black"),
    panel.grid=element_blank(),
    axis.line = element_line(),
    legend.position = "none") -> pp2

ggarrange(ps2,pp2,nrow = 1, align = "hv")

#####

