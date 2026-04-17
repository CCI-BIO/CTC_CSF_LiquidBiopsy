# 3 Seurat object merge and integration

# load packages
library(Seurat) # Seurat object related
library(SeuratWrappers) # Seurat object related
library(SeuratObject) # Seurat object related
library(Azimuth) # cell type annotation - immune cells only
library(tibble) # data table related
library(data.table) # data table related
library(readr) # data table related
library(scCustomize) # Seurat Assay version conversion

set.seed(99)

#options(Seurat.object.assay.version = "v5")

setwd("~/workspace/Nature_Comm_Use_ZCCID")
# (if sobj not loaded)
# load the sobj from 2_Cancer_Prediction output

# convert Seurat assay from V4 to V5
for (i in names(sobj)) {
  sobj[[i]] <- scCustomize::Convert_Assay(seurat_object = sobj[[i]], convert_to = "V5")
}

# merge Seurat objects
sobj_merge <- merge(x=sobj[[1]], y=sobj[c(2:33)])

# save data for backup
saveRDS(sobj_merge, file = 'out/3_Merge/merge_1.rds')

### To fix issue between Azimuth 0.5.0 and Seurat 5.3.0 (only needs to be done once, if there's an error showing ---
### "Reference assay is SCT, but query assay is RNA. Mixing SCT and non-SCT in FindTransferAnchors is not supported.")
### Uncomment the following line and execute
# devtools::install_github("satijalab/seurat", "fix/v.5.3.1", force=TRUE)

# run Azimuth
sobj_merge[["RNA"]] <- JoinLayers(sobj_merge[["RNA"]])
sobj_merge <- RunAzimuth(sobj_merge, reference = "pbmcref")
sobj_merge[["RNA"]] <- split(sobj_merge[["RNA"]], f = sobj_merge$SampleID)

# Seurat standard pre-processing workflow
sobj_merge <- NormalizeData(sobj_merge)
sobj_merge <- FindVariableFeatures(sobj_merge)
sobj_merge <- ScaleData(sobj_merge)
sobj_merge <- RunPCA(sobj_merge)
sobj_merge <- FindNeighbors(sobj_merge, dims = 1:30, reduction = "pca")
sobj_merge <- FindClusters(sobj_merge, resolution = 2, 
                           cluster.name = "unintegrated_clusters")

# Get unintegrated UMAP
sobj_merge <- RunUMAP(sobj_merge, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")

# integrated with Harmony
sobj_merge <- IntegrateLayers(
              object = sobj_merge, method = HarmonyIntegration,
              orig.reduction = "pca", new.reduction = "harmony",
              verbose = T
              )

sobj_merge <- FindNeighbors(sobj_merge, reduction = "harmony", dims = 1:30)
sobj_merge <- FindClusters(sobj_merge, resolution = 2, cluster.name = "harmony_clusters")
sobj_merge <- RunUMAP(sobj_merge, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony")

# join layers
sobj_merge <- JoinLayers(sobj_merge)

# save data
saveRDS(sobj_merge, file ='out/3_Merge/merge_2.rds')
















