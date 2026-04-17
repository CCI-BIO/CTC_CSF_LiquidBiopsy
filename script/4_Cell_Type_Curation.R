# Finalise cell type annotation

library(Seurat)
library(dittoSeq)
library(scCustomize)
library(SCpubr)
library(dplyr)
library(tidyseurat)
library(tidyHeatmap)
library(SeuratExtend)
library(decoupleR)
library(readxl)

set.seed(99)
setwd("~/workspace/Nature_Comm_Use_ZCCID")

# (if sobj not loaded)
# load the merge_2.rds from 3_Merge output folder
# sobj_merge <- readRDS("out/3_Merge/merge_2.rds")

# add QC annotation
sobj_merge$QC[sobj_merge$percent.mt >= 40] <- "fail"

# define RBC
sobj_merge$isRBC <- sobj_merge$percent.globin > 1

# define platelets
sobj_merge$isPlatelet <- ((sobj_merge$predicted.celltype.l2 == "Platelet") | 
                            (sobj_merge$singleR_hpcad == "Platelets") |
                            (sobj_merge@assays$RNA$counts["PPBP", ] > 0))

# update QC annotation
sobj_merge$QC[sobj_merge$isRBC == TRUE] <- "fail"
sobj_merge$QC[sobj_merge$isPlatelet == TRUE] <- "fail"
sobj_merge$QC[sobj_merge$scDblFinder_class == "doublet"] <- "fail"

# subset to keep Singlets
sobj_merge_s <- subset(sobj_merge, subset = QC == "pass")

# rename clusters
Idents(sobj_merge_s, drop = T) <- "harmony_clusters"
sobj_merge_s <- RenameIdents(sobj_merge_s, 
                                 "0" = "CD4+ T",
                                 "1" = "Neutrophil",
                                 "2" = "Neutrophil",
                                 "3" = "Neutrophil",
                                 "4" = "Monocyte",
                                 "5" = "Microglia/Macrophage",
                                 "6" = "CD4+ T",
                                 "7" = "Tumour",
                                 "8" = "CD8+ T",
                                 "9" = "Tumour",
                                 "10" = "DC",
                                 "11" = "Microglia P2RY12/TMEM119",
                                 "12" = "Tumour",
                                 "13" = "Neutrophil",
                                 "14" = "Neutrophil",
                                 "15" = "Tumour",
                                 "16" = "Tumour",
                                 "17" = "Tumour",
                                 "18" = "Tumour",
                                 "19" = "Microglia P2RY12/TMEM119",
                                 "20" = "Microglia P2RY12/TMEM119",
                                 "21" = "CD8+ T",
                                 "22" = "Tumour",
                                 "23" = "NK",
                                 "24" = "Tumour",
                                 "25" = "Tumour",
                                 "26" = "CD4+ T",
                                 "27" = "Tumour",
                                 "28" = "Tumour",
                                 "29" = "B",
                                 "30" = "Microglia P2RY12/TMEM119",
                                 "31" = "Neutrophil CD66b",
                                 "32" = "Tumour",
                                 "33" = "pDC",
                                 "34" = "Tumour",
                                 "35" = "Tumour",
                                 "36" = "Microglia P2RY12/TMEM119",
                                 "37" = "T",
                                 "38" = "Tumour",
                                 "39" = "Proliferation T/NK",
                                 "40" = "Unknown Immune",
                                 "41" = "Microglia P2RY12/TMEM119",
                                 "42" = "Neutrophil",
                                 "43" = "Proliferation Microglia/Macrophage",
                                 "44" = "Tumour",
                                 "45" = "Tumour",
                                 "46" = "CD4+ T",
                                 "47" = "Tumour",
                                 "48" = "DC",
                                 "49" = "Microglia P2RY12/TMEM119",
                                 "50" = "Tumour", 
                                 "51" = "DC",
                                 #"52" = "",
                                 "53" = "Monocyte",
                                 "54" = "Monocyte",
                                 #"55" = "",
                                 "56" = "DC",
                                 "57" = "Neutrophil"
            )

sobj_merge_s[["cell_type_rough"]] <- Idents(sobj_merge_s)


# refine tumour cell annotation
sobj_merge_s$cell_type_fine1 <- ifelse(sobj_merge_s$cell_type_rough == "Tumour" & sobj_merge_s$SCF_predict == "Tumour", "Tumour", "Non-malignant")

saveRDS(sobj_merge_s, file ='out/4_Cell_Type_Curation/curated_1.rds')

