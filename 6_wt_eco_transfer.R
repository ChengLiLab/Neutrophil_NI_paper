# define the cluster name of wt.eco.neu by ref, Seurat v3
#
# 2020-04-03, Qiang Shi

rm(list = ls())
setwd("/lustre/user/liclab/shiq/single_cell/luolab/neu_profile/7_wt_eco/3_transfer")

library(Seurat)

#######################################
# (1) construct ref
####################################
# matrix
ref <- readRDS("../../3_wt_ctl/1_ref/ref_neu_raw.rds")

ref.data <- ref@raw.data
dim(ref.data)
# 27998 12285

# meta
meta <- read.table("../../3_wt_ctl/2_ref_neu/ref_neu_meta.txt", header = T)
head(meta)

sum(colnames(ref.data) != rownames(meta))

# ref object
ref <- CreateSeuratObject(counts    = ref.data, 
                          meta.data = meta, 
                          project   = 'ref')

# excluding Cont
ref <- subset(ref, cluster!='Cont')

ref@assays$RNA@counts@Dim
# 27998 12155

ref <- NormalizeData(ref)

head(ref@meta.data)

table(ref@meta.data$cluster)
ref@meta.data$cluster <- factor(ref@meta.data$cluster)

rm(ref.data, meta)

#####


#######################################
# (2) construct wt.eco
####################################
# matrix
wt.eco <- readRDS("../1_wt_eco/wt_eco_neu_raw.rds")

# wt.eco object
wt.eco <- CreateSeuratObject(counts    = wt.eco@raw.data, 
                             meta.data = wt.eco@meta.data, 
                             project   = 'wt.eco')
dim(wt.eco)
# 27998 11640

wt.eco <- NormalizeData(wt.eco)

head(wt.eco@meta.data)

table(wt.eco$tissue)

#####


#######################################
# (3) prepare features, marker gene
####################################
deg <- read.csv(file      = "../../3_wt_ctl/2_ref_neu/neu_cluster_deg_log1.2.csv", 
                row.names = 1, stringsAsFactors = F)

table(deg$cluster)

library(dplyr)
deg <- deg %>% group_by(cluster) %>% top_n(100, avg_logFC)

#####


#######################################
# (4) transfer
####################################
# anchor
transfer.anchors <- FindTransferAnchors(reference = ref, 
                                        query     = wt.eco, 
                                        features  = unique(c(deg$gene)),
                                        dims      = 1:15)

# predict
prediction <- TransferData(anchorset = transfer.anchors, 
                           refdata   = ref$cluster,
                           dims      = 1:15)

head(wt.eco@meta.data)
sum(rownames(wt.eco@meta.data)!=rownames(prediction))

wt.eco$prediction.score.max <- prediction$prediction.score.max
wt.eco$predicted.id <- prediction$predicted.id

# unassigned
plot(density(wt.eco$prediction.score.max))

sum(wt.eco$prediction.score.max<0.5)
# 1266
# 1112

wt.eco$predicted.cluster <- wt.eco$predicted.id

wt.eco@meta.data[wt.eco$prediction.score.max<0.5, 'predicted.cluster'] <- 'Unassigned'

table(wt.eco@meta.data[wt.eco$prediction.score.max<0.5, 'tissue'])
table(wt.eco$tissue, wt.eco$predicted.cluster)

table(ref$cluster)
table(wt.eco$predicted.cluster)

#####


#######################################
# (5) save
####################################
head(wt.eco@meta.data)

# save prediction
write.table(prediction, "./wt_eco_predictions_marker.txt", quote = FALSE)

# save metadata
write.table(wt.eco@meta.data, "./wt_eco_neu_meta_transfer_new.txt", quote = FALSE)

#####

