#!/usr/bin/env Rscript
# Seurat preprocessing pipeline for hESC_scRNASeq (iSCORE-PD irradiation)
# Methods from: bioRxiv 10.1101/2024.02.12.579917v2 (PMC10888955)
# "Seurat v4 according to default parameters for normalization...
#  Droplets with more than 15% mitochondrial reads detected were excluded"

library(Seurat)

data_dir <- "/mnt/c/Users/jesse/Downloads/CRN_METADATA_TABLES/test/matrix_GEX"
out_dir  <- "/mnt/c/Users/jesse/Downloads/CRN_METADATA_TABLES/test"

# 1. Load 10x data and create Seurat object
cat("Loading 10x data from:", data_dir, "\n")
counts <- Read10X(data.dir = data_dir)
obj <- CreateSeuratObject(counts = counts, project = "hESC_scRNASeq")
cat("After CreateSeuratObject:", ncol(obj), "cells,", nrow(obj), "genes\n")

# 2. Calculate percent mitochondrial
obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")

# QC summary before filtering
cat("\nPre-filter QC summary:\n")
print(summary(obj$nFeature_RNA))
cat("percent.mt summary:\n")
print(summary(obj$percent.mt))

# 3. Save QC violin plot
png(file.path(out_dir, "QC_violin.png"), width = 1200, height = 400)
VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

# 4. QC filtering (per iSCORE-PD hESC methods: only percent.mt > 15% excluded)
obj <- subset(obj, subset = percent.mt < 15)
cat("After QC filtering:", ncol(obj), "cells,", nrow(obj), "genes\n")

# 5. Normalize
obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)

# 6. Find variable features
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000)

# 7. Scale data (all genes)
all.genes <- rownames(obj)
obj <- ScaleData(obj, features = all.genes)

# 8. PCA
obj <- RunPCA(obj, features = VariableFeatures(object = obj))

# Save elbow plot
png(file.path(out_dir, "ElbowPlot.png"), width = 600, height = 400)
ElbowPlot(obj)
dev.off()

# 9. Clustering
obj <- FindNeighbors(obj, dims = 1:10)
obj <- FindClusters(obj, resolution = 0.5)

# 10. UMAP
obj <- RunUMAP(obj, dims = 1:10)

# Save UMAP plot
png(file.path(out_dir, "UMAP.png"), width = 700, height = 600)
DimPlot(obj, reduction = "umap")
dev.off()

# Final summary
cat("\n=== Final Seurat Object ===\n")
print(obj)
cat("\nCluster distribution:\n")
print(table(Idents(obj)))

# 11. Save RDS
rds_path <- file.path(out_dir, "irradiated_hESC_object.rds")
saveRDS(obj, file = rds_path)
cat("\nSaved RDS to:", rds_path, "\n")
