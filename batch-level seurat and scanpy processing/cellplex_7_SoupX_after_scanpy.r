library(SoupX)
library(Seurat)

tod <- Seurat::Read10X('/share/crsp/lab/dalawson/share/1_Lab_sequencing/1_Tatyana/Cellranger/large_and_smallest_tumors_cellplex7/outs/multi/count/raw_feature_bc_matrix')

toc_1 <- Seurat::Read10X('/share/crsp/lab/dalawson/share/1_Lab_sequencing/1_Tatyana/Cellranger/large_and_smallest_tumors_cellplex7/outs/per_sample_outs/possible_smallest_tumor/count/sample_filtered_feature_bc_matrix')

toc_2 <- Seurat::Read10X('/share/crsp/lab/dalawson/share/1_Lab_sequencing/1_Tatyana/Cellranger/large_and_smallest_tumors_cellplex7/outs/per_sample_outs/tumor_081720/count/sample_filtered_feature_bc_matrix')

toc_3 <- Seurat::Read10X('/share/crsp/lab/dalawson/share/1_Lab_sequencing/1_Tatyana/Cellranger/large_and_smallest_tumors_cellplex7/outs/per_sample_outs/tumor_120519/count/sample_filtered_feature_bc_matrix')


toc_1 <- toc_1[[1]]
toc_2 <- toc_2[[1]]
toc_3 <- toc_3[[1]]

tod <- tod[[1]]

toc <- cbind(toc_1,toc_2,toc_3)
toc <- toc[ order(match(rownames(toc), rownames(tod))), ]

cluster_labels <- read.csv("cellplex7.bigsur.celltype.clusters.csv", header = T)
colnames(cluster_labels) <- c("Barcode","Cluster")

toc <- toc[,which(colnames(toc) %in% cluster_labels$Barcode)]

sc <- SoupChannel(tod,toc)

sc <- setClusters(sc,setNames(cluster_labels$Cluster, cluster_labels$Barcode))


sc <- autoEstCont(sc) #rho = 0.02
out <- adjustCounts(sc)

obj <- CreateSeuratObject(counts = out)

saveRDS(obj, "cellplex7_soupX.rds")
