# CNS CTC CSF Paper
# CCIA-Liquid Biopsy

# 1 QC

setwd('~/workspace/Nature_Comm_Use_ZCCID')

# load packages
library(Seurat)
library(reticulate)
library(anndata)
library(sceasy)
library(SeuratData)
library(SeuratWrappers)
library(Azimuth)
library(ggplot2)
library(patchwork)
library(scCustomize)
library(dplyr)
library(SingleCellExperiment)
library(SingleR)
library(celldex)
library(scuttle)
library(scDblFinder)

# load python packages
sc <- import("scanpy", convert = FALSE)
ct <- import("celltypist", convert = FALSE)

set.seed(99)

# --------------------------
# load data
# --------------------------

# get file path
file_path <- "data/seurat"
sample_list <- list.files(file_path)
sample_name <- gsub(pattern = "_Seurat.rds",
                    replacement = "", x = sample_list)

sample_name <- gsub(pattern = "-", replacement = "_", sample_name)

names(sample_list) <- sample_name

# annotate sample with disease type
disease_long <- c("zccs905_CSF383372"="Pineoblastoma", 
                  "zccs905_CSF387478"="Pineoblastoma", 
                  "zccs905_CSF396065"="Pineoblastoma",
                  "zccs905_CSF403779"="Pineoblastoma",
                  "zccs905_CSF428379"="Pineoblastoma",
                  "zccs905_CSF434052"="Pineoblastoma",
                  "zccs905_CSF434893"="Pineoblastoma",
                  
                  "zccs906_CSF384921"="Posterior fossa A ependymoma",
                  "zccs906_CSF425809"="Posterior fossa A ependymoma",
                  "zccs906_CSF432740"="Posterior fossa A ependymoma",
                  
                  "zccs907_CSF383908"="Infant-type hemispheric glioma", 
                  
                  "zccs920_CSF387517"="Pilocytic astrocytoma",
                  
                  "zccs929_CSF390503"="Medulloblastoma", 
                  
                  "zccs934_CSF394836"="Pilocytic astrocytoma",
                  
                  "zccs948_CSF398038"="Posterior fossa A ependymoma",
                  
                  "zccs1344_CSF395512"= "Diffuse leptomeningeal glioneuronal tumour",
                  
                  "zccs527_CSF337322"= "Diffuse midline glioma",
                  
                  "zccs1520_CSF343889"= "Medulloblastoma",
                  
                  "zccs1576_CSF341894"= "Medulloblastoma",
                  "zccs1576_CSF360038"="Medulloblastoma", 
                  
                  "zccs1649_CSF353421"= "Diffuse midline glioma",
                  "zccs1649_CSF361231"="Diffuse midline glioma", 
                  
                  "zccs1681_CSF363102"= "Atypical teratoid rhabdoid tumour",
                  "zccs1681_CSF384438"="Atypical teratoid rhabdoid tumour",
                  "zccs1681_CSF387541"="Atypical teratoid rhabdoid tumour",
                  "zccs1681_CSF391810"="Atypical teratoid rhabdoid tumour",
                  "zccs1681_CSF395366"="Atypical teratoid rhabdoid tumour",
                  "zccs1681_CSF423677"="Atypical teratoid rhabdoid tumour",
                  "zccs1681_CSF431954"="Atypical teratoid rhabdoid tumour",
                  
                  "zccs1695_CSF367534"="Diffuse midline glioma",
                  
                  "zccs3414_CSF376420"= "Medulloblastoma",
                  "zccs3414_CSF413089"= "Medulloblastoma",
                  
                  "zccs1713_CSF378332"="Glioneuronal tumour")

disease_short <- c("zccs905_CSF383372"="PB", 
                  "zccs905_CSF387478"="PB", 
                  "zccs905_CSF396065"="PB",
                  "zccs905_CSF403779"="PB",
                  "zccs905_CSF428379"="PB",
                  "zccs905_CSF434052"="PB",
                  "zccs905_CSF434893"="PB",
                  
                  "zccs906_CSF384921"="PFA",
                  "zccs906_CSF425809"="PFA",
                  "zccs906_CSF432740"="PFA",
                  
                  "zccs907_CSF383908"="IHG", 
                  
                  "zccs920_CSF387517"="PA",
                  
                  "zccs929_CSF390503"="MB", 
                  
                  "zccs934_CSF394836"="PA",
                  
                  "zccs948_CSF398038"="PFA",
                  
                  "zccs1344_CSF395512"= "DLGNT",
                  
                  "zccs527_CSF337322"= "DMG",
                  
                  "zccs1520_CSF343889"= "MB",
                  
                  "zccs1576_CSF341894"= "MB",
                  "zccs1576_CSF360038"="MB", 
                  
                  "zccs1649_CSF353421"= "DMG",
                  "zccs1649_CSF361231"="DMG", 
                  
                  "zccs1681_CSF363102"="ATRT",
                  "zccs1681_CSF384438"="ATRT",
                  "zccs1681_CSF387541"="ATRT",
                  "zccs1681_CSF391810"="ATRT",
                  "zccs1681_CSF395366"="ATRT",
                  "zccs1681_CSF423677"="ATRT",
                  "zccs1681_CSF431954"="ATRT",
                  
                  "zccs1695_CSF367534"="DMG",
                  
                  "zccs3414_CSF376420"= "MB",
                  "zccs3414_CSF413089"= "MB",
                  
                  "zccs1713_CSF378332"="GNT")

saveRDS(disease_long, file = "out/1_QC/disease_long_list.rds")
saveRDS(disease_short, file = "out/1_QC/disease_short_list.rds")


# load the wetlab data
wetlab_meta <- read.csv("misc/Nature_Comm_CNS_sample_table_20251205_ZCCID.csv")
wetlab_meta <- wetlab_meta %>% mutate(uni_name = paste(zcc_id, bioid, sep = "_"))

# Load and process each sobj
sobj <- c()

# - SingleR
ref_bed <- celldex::BlueprintEncodeData()
ref_hpcad <- celldex::HumanPrimaryCellAtlasData()

# - celltypist
model_immue = ct$models$Model$load(model = 'Immune_All_Low.pkl')
model_brain = ct$models$Model$load(model = 'Developing_Human_Brain.pkl')

# Processing
for (i in names(sample_list)) {
  print("Processing sample:")
  print(i)
  
  sobj[[i]] <- readRDS(paste(file_path, sample_list[i], sep = "/"))
  
  sobj[[i]]$QC <- "pass"
  sobj[[i]]$QC[sobj[[i]]$nFeature_RNA < 200] <- "fail"
  
  sobj[[i]]<- subset(sobj[[i]], subset = QC == "pass", invert = F)
  
  Project(sobj[[i]]) <- i
  sobj[[i]]$CellID <- colnames(sobj[[i]])
  sobj[[i]]$SampleID <- i
  sobj[[i]]$PatientID <- strsplit(i, split = "_")[[1]][1]
  sobj[[i]]$Disease_long <- disease_long[i] |> as.character()
  sobj[[i]]$Disease_short <- disease_short[i] |> as.character()
  
  sobj[[i]]$bioid <- wetlab_meta$bioid[wetlab_meta$uni_name == i]
  sobj[[i]]$public_subject_id <- wetlab_meta$public_subject_id[wetlab_meta$uni_name == i]
  sobj[[i]]$zcc_id <- wetlab_meta$zcc_id[wetlab_meta$uni_name == i]
  sobj[[i]]$lm_enrollmnet <- wetlab_meta$lm_enrollmnet[wetlab_meta$uni_name == i]
  sobj[[i]]$lm_mortality <- wetlab_meta$lm_mortality[wetlab_meta$uni_name == i]
  sobj[[i]]$lm_biomat_id <- wetlab_meta$lm_biomat_id[wetlab_meta$uni_name == i]
  sobj[[i]]$lm_biomat_collection <- wetlab_meta$lm_biomat_collection[wetlab_meta$uni_name == i]
  sobj[[i]]$V.mL. <- wetlab_meta$V.mL.[wetlab_meta$uni_name == i]
  sobj[[i]]$RBC <- wetlab_meta$RBC[wetlab_meta$uni_name == i]
  sobj[[i]]$Enrichment <- wetlab_meta$Enrichment[wetlab_meta$uni_name == i]
  sobj[[i]]$Preservation <- wetlab_meta$Preservation[wetlab_meta$uni_name == i]
  sobj[[i]]$Cytology <- wetlab_meta$Cytology[wetlab_meta$uni_name == i]
  
  # create adata for calculation in the python envrionment
  adata <- sceasy::convertFormat(sobj[[i]], from="seurat", to="anndata", main_layer="counts", drop_single_values=FALSE)
  
  # mt genes
  sobj[[i]] <- Seurat::PercentageFeatureSet(sobj[[i]], 
                                            pattern = "^MT[-|.]", 
                                            col.name = "percent.mt") 
  
  # ribosomal genes
  sobj[[i]] <- Seurat::PercentageFeatureSet(sobj[[i]], 
                                            pattern = "^RP[SL]",
                                            col.name = "percent.ribo")
  
  # hemoglobin genes (but not HBP)
  sobj[[i]] <- Seurat::PercentageFeatureSet(sobj[[i]],
                                            pattern = "^HB[^(P)]",
                                            col.name = "percent.globin")
  
  # define doublet
  sobj.sce <- as.SingleCellExperiment(sobj[[i]])
  sobj.sce <- scDblFinder(sobj.sce)
  sobj[[i]]$scDblFinder_class <- sobj.sce$scDblFinder.class
  sobj[[i]]$scDblFinder_score <- sobj.sce$scDblFinder.score
  sobj[[i]]$scDblFinder_weighted <- sobj.sce$scDblFinder.weighted
  sobj[[i]]$scDblFinder_cxds_score <- sobj.sce$scDblFinder.cxds_score
  
  print("scDblFinder done!")
  
  # singleR prediction
  sobj.sce <- scuttle::logNormCounts(sobj.sce)
  
  pred.bed <- SingleR(test=sobj.sce, ref=ref_bed, 
                      labels=ref_bed$label.main, de.method="wilcox")
  
  pred.hpcad <- SingleR(test=sobj.sce, ref=ref_hpcad, 
                        labels=ref_hpcad$label.main, de.method="wilcox")
  
  sobj[[i]]@misc[["singleR_bed"]] <- pred.bed
  sobj[[i]]@misc[["singleR_hpcad"]] <- pred.hpcad
  
  pred.bed.labels <- pred.bed$labels
  pred.hpcad.labels <- pred.hpcad$labels
  
  names(pred.bed.labels) <- rownames(pred.bed$labels)
  names(pred.hpcad.labels) <- rownames(pred.hpcad$labels)
  
  sobj[[i]] <- AddMetaData(sobj[[i]], 
                           metadata = pred.bed.labels,
                           col.name = "singleR_bed")
  
  sobj[[i]] <- AddMetaData(sobj[[i]], 
                           metadata = pred.hpcad.labels,
                           col.name = "singleR_hpcad")
  
  print("SingleR done!")
  
  # Celltypist annotation
  adata$layers["counts"] = adata$X$copy()
  sc$pp$normalize_per_cell(adata, counts_per_cell_after=10**4)
  sc$pp$log1p(adata)
  adata$X = adata$X$toarray()
  
  predictions = ct$annotate(adata, model = 'Immune_All_Low.pkl', majority_voting = TRUE)
  predictions_adata = predictions$to_adata()
  adata$obs["ct_immune"] = predictions_adata$obs$loc[adata$obs$index, "majority_voting"]
  adata$obs["ct_immune_score"] = predictions_adata$obs$loc[adata$obs$index, "conf_score"]
  
  predictions = ct$annotate(adata, model = 'Developing_Human_Brain.pkl', majority_voting = TRUE)
  predictions_adata = predictions$to_adata()
  adata$obs["ct_brain"] = predictions_adata$obs$loc[adata$obs$index, "majority_voting"]
  adata$obs["ct_brain_score"] = predictions_adata$obs$loc[adata$obs$index, "conf_score"]
  
  tmp_meta <- py_to_r(adata$obs)
  tmp_meta <- tmp_meta[, c("ct_immune", "ct_immune_score", "ct_brain", "ct_brain_score")]
  identical(rownames(tmp_meta), colnames(sobj[[i]]))
  sobj[[i]] <- AddMetaData(sobj[[i]], tmp_meta)
  
  print("Celltypist done!")
  
  print("=======")
  print("Done!")
}

# save data
saveRDS(sobj, file = "out/1_QC/Prepocessed.rds")

# 


