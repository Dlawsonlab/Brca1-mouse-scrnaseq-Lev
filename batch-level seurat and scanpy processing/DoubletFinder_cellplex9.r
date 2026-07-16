library(Seurat)
library(DoubletFinder)

#R/4.2.2
#do in Rstudio 

setwd("/share/crsp/lab/dalawson/tzhuravl/CellPlex_Tumors/Scanpy_Seurat/cellplex_9")

X <- readRDS("/share/crsp/lab/dalawson/tzhuravl/CellPlex_Tumors/Scanpy_Seurat/cellplex_9/cellplex9_soupX.rds")
sweep.res <- paramSweep_v3(X, PCs = 1:10, sct = TRUE) 
sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
bcmvn <- find.pK(sweep.stats)

nExp_poi <- round(0.05*nrow(X@meta.data))
annotations <- X@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)

X <- doubletFinder_v3(X, PCs = 1:10, pN = 0.25, pK = 0.01, nExp = nExp_poi, reuse.pANN = FALSE, sct = TRUE)
singlet_or_doublet = colnames(X@meta.data)[grepl("DF.classification", colnames(X@meta.data))]

DimPlot(X, group.by = singlet_or_doublet) + NoAxes()
VlnPlot(X, features = "nFeature_RNA", group.by = singlet_or_doublet, pt.size = 0.1)

obj <- readRDS("/share/crsp/lab/dalawson/tzhuravl/CellPlex_Tumors/Scanpy_Seurat/cellplex_9/cellplex9_soupX.rds")

obj[["singlet_or_doublet"]] <- X@meta.data$DF.classification
DimPlot(obj, group.by = "singlet_or_doublet") + NoAxes()
VlnPlot(obj, features = "nFeature_RNA", group.by = "singlet_or_doublet", pt.size = 0.1)

saveRDS(obj, "/share/crsp/lab/dalawson/tzhuravl/CellPlex_Tumors/Scanpy_Seurat/cellplex_9/cellplex9_soupX_doubletfinder.rds")
