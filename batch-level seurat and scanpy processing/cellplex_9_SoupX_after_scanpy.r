library(SoupX)
library(Seurat)

tod <- Seurat::Read10X('/share/crsp/lab/dalawson/share/1_Lab_sequencing/1_Tatyana/Cellranger/both_blood_sorted_epithelial_very_small_tumor/outs/multi/count/raw_feature_bc_matrix')

toc_1 <- Seurat::Read10X('/share/crsp/lab/dalawson/share/1_Lab_sequencing/1_Tatyana/Cellranger/both_blood_sorted_epithelial_very_small_tumor/outs/per_sample_outs/blood_MUT_1/count/sample_filtered_feature_bc_matrix')

toc_2 <- Seurat::Read10X('/share/crsp/lab/dalawson/share/1_Lab_sequencing/1_Tatyana/Cellranger/both_blood_sorted_epithelial_very_small_tumor/outs/per_sample_outs/blood_MUT_2/count/sample_filtered_feature_bc_matrix')

toc_3 <- Seurat::Read10X('/share/crsp/lab/dalawson/share/1_Lab_sequencing/1_Tatyana/Cellranger/both_blood_sorted_epithelial_very_small_tumor/outs/per_sample_outs/blood_MUT_3/count/sample_filtered_feature_bc_matrix')

toc_4 <- Seurat::Read10X('/share/crsp/lab/dalawson/share/1_Lab_sequencing/1_Tatyana/Cellranger/both_blood_sorted_epithelial_very_small_tumor/outs/per_sample_outs/blood_WT_1/count/sample_filtered_feature_bc_matrix')

toc_5 <- Seurat::Read10X('/share/crsp/lab/dalawson/share/1_Lab_sequencing/1_Tatyana/Cellranger/both_blood_sorted_epithelial_very_small_tumor/outs/per_sample_outs/blood_WT_2/count/sample_filtered_feature_bc_matrix')

toc_6 <- Seurat::Read10X('/share/crsp/lab/dalawson/share/1_Lab_sequencing/1_Tatyana/Cellranger/both_blood_sorted_epithelial_very_small_tumor/outs/per_sample_outs/blood_WT_3/count/sample_filtered_feature_bc_matrix')

toc_7 <- Seurat::Read10X('/share/crsp/lab/dalawson/share/1_Lab_sequencing/1_Tatyana/Cellranger/both_blood_sorted_epithelial_very_small_tumor/outs/per_sample_outs/else_epithelial_MUT/count/sample_filtered_feature_bc_matrix')

toc_8 <- Seurat::Read10X('/share/crsp/lab/dalawson/share/1_Lab_sequencing/1_Tatyana/Cellranger/both_blood_sorted_epithelial_very_small_tumor/outs/per_sample_outs/else_epithelial_WT/count/sample_filtered_feature_bc_matrix')

toc_9 <- Seurat::Read10X('/share/crsp/lab/dalawson/share/1_Lab_sequencing/1_Tatyana/Cellranger/both_blood_sorted_epithelial_very_small_tumor/outs/per_sample_outs/intermediate_MUT/count/sample_filtered_feature_bc_matrix')

toc_10 <- Seurat::Read10X('/share/crsp/lab/dalawson/share/1_Lab_sequencing/1_Tatyana/Cellranger/both_blood_sorted_epithelial_very_small_tumor/outs/per_sample_outs/intermediate_WT/count/sample_filtered_feature_bc_matrix')

toc_11 <- Seurat::Read10X('/share/crsp/lab/dalawson/share/1_Lab_sequencing/1_Tatyana/Cellranger/both_blood_sorted_epithelial_very_small_tumor/outs/per_sample_outs/very_small_tumor_966/count/sample_filtered_feature_bc_matrix')


toc_1 <- toc_1[[1]]
toc_2 <- toc_2[[1]]
toc_3 <- toc_3[[1]]
toc_4 <- toc_4[[1]]
toc_5 <- toc_5[[1]]
toc_6 <- toc_6[[1]]
toc_7 <- toc_7[[1]]
toc_8 <- toc_8[[1]]
toc_9 <- toc_9[[1]]
toc_10 <- toc_10[[1]]
toc_11 <- toc_11[[1]]

tod <- tod[[1]]

toc <- cbind(toc_1,toc_2,toc_3,toc_4,toc_5,toc_6,toc_7,toc_8,toc_9,toc_10,toc_11)
toc <- toc[ order(match(rownames(toc), rownames(tod))), ]

cluster_labels <- read.csv("cellplex9.bigsur.celltype.clusters.csv", header = T)
colnames(cluster_labels) <- c("Barcode","Cluster")

toc <- toc[,which(colnames(toc) %in% cluster_labels$Barcode)]

sc <- SoupChannel(tod,toc)

sc <- setClusters(sc,setNames(cluster_labels$Cluster, cluster_labels$Barcode))


sc <- autoEstCont(sc) #rho = 0.02
out <- adjustCounts(sc)

obj <- CreateSeuratObject(counts = out)

saveRDS(obj, "cellplex9_soupX.rds")
