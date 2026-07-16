library(SoupX)
library(Seurat)

tod <- Seurat::Read10X('/share/crsp/lab/dalawson/share/1_Lab_sequencing/1_Tatyana/Cellranger/nMFP_plus_pMFP_cellplex3/outs/multi/count/raw_feature_bc_matrix')

toc_1 <- Seurat::Read10X('/share/crsp/lab/dalawson/share/1_Lab_sequencing/1_Tatyana/Cellranger/nMFP_plus_pMFP_cellplex3/outs/per_sample_outs/nMFP_1/count/sample_filtered_feature_bc_matrix')

toc_2 <- Seurat::Read10X('/share/crsp/lab/dalawson/share/1_Lab_sequencing/1_Tatyana/Cellranger/nMFP_plus_pMFP_cellplex3/outs/per_sample_outs/nMFP_2/count/sample_filtered_feature_bc_matrix')

toc_3 <- Seurat::Read10X('/share/crsp/lab/dalawson/share/1_Lab_sequencing/1_Tatyana/Cellranger/nMFP_plus_pMFP_cellplex3/outs/per_sample_outs/nMFP_3/count/sample_filtered_feature_bc_matrix')

toc_4 <- Seurat::Read10X('/share/crsp/lab/dalawson/share/1_Lab_sequencing/1_Tatyana/Cellranger/nMFP_plus_pMFP_cellplex3/outs/per_sample_outs/nMFP_4/count/sample_filtered_feature_bc_matrix')

toc_5 <- Seurat::Read10X('/share/crsp/lab/dalawson/share/1_Lab_sequencing/1_Tatyana/Cellranger/nMFP_plus_pMFP_cellplex3/outs/per_sample_outs/nMFP_5/count/sample_filtered_feature_bc_matrix')

toc_6 <- Seurat::Read10X('/share/crsp/lab/dalawson/share/1_Lab_sequencing/1_Tatyana/Cellranger/nMFP_plus_pMFP_cellplex3/outs/per_sample_outs/pMFP_1/count/sample_filtered_feature_bc_matrix')

toc_7 <- Seurat::Read10X('/share/crsp/lab/dalawson/share/1_Lab_sequencing/1_Tatyana/Cellranger/nMFP_plus_pMFP_cellplex3/outs/per_sample_outs/pMFP_2/count/sample_filtered_feature_bc_matrix')

toc_8 <- Seurat::Read10X('/share/crsp/lab/dalawson/share/1_Lab_sequencing/1_Tatyana/Cellranger/nMFP_plus_pMFP_cellplex3/outs/per_sample_outs/pMFP_3/count/sample_filtered_feature_bc_matrix')

toc_9 <- Seurat::Read10X('/share/crsp/lab/dalawson/share/1_Lab_sequencing/1_Tatyana/Cellranger/nMFP_plus_pMFP_cellplex3/outs/per_sample_outs/pMFP_4/count/sample_filtered_feature_bc_matrix')

toc_10 <- Seurat::Read10X('/share/crsp/lab/dalawson/share/1_Lab_sequencing/1_Tatyana/Cellranger/nMFP_plus_pMFP_cellplex3/outs/per_sample_outs/pMFP_5/count/sample_filtered_feature_bc_matrix')


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

tod <- tod[[1]]

toc <- cbind(toc_1,toc_2,toc_3,toc_4,toc_5,toc_6,toc_7,toc_8,toc_9,toc_10)
toc <- toc[ order(match(rownames(toc), rownames(tod))), ]

cluster_labels <- read.csv("cellplex3.bigsur.celltype.clusters.csv", header = T)
colnames(cluster_labels) <- c("Barcode","Cluster")

toc <- toc[,which(colnames(toc) %in% cluster_labels$Barcode)]

sc <- SoupChannel(tod,toc)

sc <- setClusters(sc,setNames(cluster_labels$Cluster, cluster_labels$Barcode))


sc <- autoEstCont(sc) #rho = 0.05
out <- adjustCounts(sc)

obj <- CreateSeuratObject(counts = out)

saveRDS(obj, "cellplex3_soupX.rds")
