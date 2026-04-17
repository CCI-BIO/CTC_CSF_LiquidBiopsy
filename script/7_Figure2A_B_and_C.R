# CNS CTC CSF Paper
# Figure 2A, 2B and 2C

library(Seurat)
library(dittoSeq)
library(cowplot)
library(ggplot2)
library(dplyr)
library(tidyr)
library(copykat)
library(tibble)
library(patchwork)
library(ComplexHeatmap)
library(presto)

set.seed(99)
setwd("~/workspace/Nature_Comm_Use_ZCCID")


### CNV can be calculated by copykat using the code below
# raw_counts <- as.matrix(zccs905@assays$RNA$counts)
# copykat.zccs905 <- copykat(rawmat=raw_counts,
#                          id.type="S",
#                          cell.line = "no",
#                          ngene.chr=5,
#                          win.size=25,
#                          KS.cut=0.1,
#                          sam.name="zccs905",
#                          distance="euclidean",
#                          norm.cell.names="",
#                          output.seg="FALSE",
#                          plot.genes="FALSE",
#                          genome="hg20",
#                          n.cores=20)
# 
# saveRDS(copykat.zccs905, file = 'misc/copykat_results/copykat_output.rds')




# load copykat results
copykat.zccs905 <- readRDS("misc/copykat_results/copykat_output.rds")

# check if cell IDs match those in the Seurat object
cell_id_seurat <- colnames(zccs905)

cell_id_seurat <- gsub("_19$", "_2", cell_id_seurat)
cell_id_seurat <- gsub("_20$", "_1", cell_id_seurat)
cell_id_seurat <- gsub("_21$", "_3", cell_id_seurat)
cell_id_seurat <- gsub("_22$", "_4", cell_id_seurat)
cell_id_seurat <- gsub("_23$", "_5", cell_id_seurat)
cell_id_seurat <- gsub("_24$", "_6", cell_id_seurat)
cell_id_seurat <- gsub("_25$", "_7", cell_id_seurat)

CellID_check <- identical(cell_id_seurat |> sort(), copykat.zccs905$prediction$cell.names |> sort())

if (CellID_check) {
  print("Cell IDs in the copykat results and in the Seurat object are identical.")
}else{
  print("Found mismatched cell IDs between copykat results and the Seurat object.")
}




# visualise data with a heatmap
pred.test <- data.frame(copykat.zccs905$prediction)
pred.test <- pred.test %>% dplyr::filter(copykat.pred != "not.defined")  ##remove undefined cells
CNA.test <- data.frame(copykat.zccs905$CNAmat)

# get the data matrix
mat <- t(CNA.test[,4:ncol(CNA.test)])
com.preN <- pred.test$copykat.pred

# group cells by days from diagnosis
enroll <- unique(zccs905$lm_enrollmnet)
temp_time_name <- unique(zccs905$lm_biomat_collection)
temp_time <- as.Date(temp_time_name, format = "%d/%m/%Y") - as.Date(enroll, format = "%d/%m/%Y")
names(temp_time) <- temp_time_name
temp_time <- temp_time[order(temp_time)]

Idents(zccs905) <- "lm_biomat_collection"
zccs905 <- RenameIdents(zccs905, temp_time)
zccs905[["Days"]] <- Idents(zccs905)

# annotate cells with CNV prediction
temp_meta <- FetchData(zccs905, vars = c("SampleID", "Days")) %>% rownames_to_column(var = "cell.names")

temp_meta$cell.names.batch <- temp_meta$cell.names

temp_meta$cell.names.batch <- gsub("_19$", "_2", temp_meta$cell.names.batch)
temp_meta$cell.names.batch <- gsub("_20$", "_1", temp_meta$cell.names.batch)
temp_meta$cell.names.batch <- gsub("_21$", "_3", temp_meta$cell.names.batch)
temp_meta$cell.names.batch <- gsub("_22$", "_4", temp_meta$cell.names.batch)
temp_meta$cell.names.batch <- gsub("_23$", "_5", temp_meta$cell.names.batch)
temp_meta$cell.names.batch <- gsub("_24$", "_6", temp_meta$cell.names.batch)
temp_meta$cell.names.batch <- gsub("_25$", "_7", temp_meta$cell.names.batch)

pred.test.meta <- left_join(pred.test, temp_meta, by = c("cell.names" = "cell.names.batch"))

# heatmap row annotation
days_id_leftAnno <- pred.test.meta$Days
col_days_id <- ggsci::pal_jama()(7)
names(col_days_id) <- unique(days_id_leftAnno)
names(days_id_leftAnno) <- paste0("X", pred.test.meta$cell.names)

# get clustering info
hc.umap <- cutree(copykat.zccs905$hclustering,2)
hc.cut <- cutree(copykat.zccs905$hclustering,2)
hcc <- readRDS("misc/copykat_results/copykat_hcc_result.rds")
roworder <- hcc$order

# Chr color assignment
chr_cols <- c("1" = "black", "2" = "grey", "3" = "black", "4" = "grey", "5" = "black", "6" = "grey", "7" = "black", 
              "8" = "grey", "9" = "black", "10" = "grey", "11" = "black", "12" = "grey", "13" = "black", 
              "14" = "grey", "15" = "black", "16" = "grey", "17" = "black", "18" = "grey", "19" = "black", 
              "20" = "grey", "21"= "black", "22" = "grey", "23" = "black", "chrX" = "grey", "chrY" = "black", "chrMT" = "grey")

# heatmap column annotation
topAnnotation <- HeatmapAnnotation(Chr = CNA.test$chrom, 
                                   show_legend = F, 
                                   show_annotation_name = T, 
                                   col = list(Chr = chr_cols))

# use the following leftAnnotation for displaying days apart among the samples
leftAnnotation <- rowAnnotation(Days = days_id_leftAnno, 
                                Prediction=com.preN,
                                col = list(Days = col_days_id, 
                                           Prediction = c("diploid" = "#5773CC", "aneuploid" = "#FFB900") ) )

col_fun = circlize::colorRamp2(c(min(mat), 0, max(mat)), c("#08088A", "white", "#B40431"))


ht <- ComplexHeatmap::Heatmap(mat,
               left_annotation = leftAnnotation,
               row_order = roworder, # roworder is by copykat clustering
               top_annotation=topAnnotation,
               column_split = CNA.test$chrom,
               column_gap =  unit(0, "mm"),
               column_title = c(1:23),
               column_title_rot = 45,
               col=col_fun,
               name = "Score",
               show_row_names=F,
               show_column_names=F,
               cluster_columns=F,
               clustering_distance_columns = "euclidean",
               clustering_method_columns = "ward.D2",
               use_raster = T
)


pdf("out/5_Figures/Figure_2C.pdf",width = 8, height = 5.5) 
draw(ht)
dev.off()





# ----- UMAP plots

# add CNV prediction to the Seurat object
temp <- FetchData(object = zccs905, vars = c("orig.ident", "SampleID")) |> tibble::rownames_to_column("cell.names")

copykat_predict <- copykat.zccs905$prediction

copykat_predict$cell.names.batch <- copykat_predict$cell.names

copykat_predict$cell.names.batch <- gsub("_2$", "_19", copykat_predict$cell.names.batch)
copykat_predict$cell.names.batch <- gsub("_1$", "_20", copykat_predict$cell.names.batch)
copykat_predict$cell.names.batch <- gsub("_3$", "_21", copykat_predict$cell.names.batch)
copykat_predict$cell.names.batch <- gsub("_4$", "_22", copykat_predict$cell.names.batch)
copykat_predict$cell.names.batch <- gsub("_5$", "_23", copykat_predict$cell.names.batch)
copykat_predict$cell.names.batch <- gsub("_6$", "_24", copykat_predict$cell.names.batch)
copykat_predict$cell.names.batch <- gsub("_7$", "_25", copykat_predict$cell.names.batch)

temp <- dplyr::left_join(temp, copykat_predict, by = c("cell.names" = "cell.names.batch"))
temp_pred <- temp %>% pull(copykat.pred)
names(temp_pred) <- temp %>% pull(cell.names)

zccs905 <-AddMetaData(zccs905, metadata = temp_pred, col.name = "CNV Prediction")


# annotation formatting
Idents(zccs905) <- "CNV Prediction"
zccs905 <- RenameIdents(zccs905, 
                      "aneuploid" = "Aneuploid",
                      "diploid" = "Diploid",
                      "not.defined" = "Not defined")
zccs905[["CNV_predict"]] <- Idents(zccs905)


Idents(zccs905) <- "SCF_predict"
zccs905 <- RenameIdents(zccs905, 
                      "Normal" = "Non-malignant")
zccs905[["SCF_predict2"]] <- Idents(zccs905)




# umap plot showing CNV prediction results
p_cnv_umap <-dittoSeq::dittoDimPlot(zccs905,
                                    var = c("CNV_predict"), main= NULL, legend.title = "CNV Prediction",
                                    reduction = "umap.harmony", theme = theme_cowplot() + theme(legend.position = "bottom", 
                                                                                                legend.justification = "centre", 
                                                                                                legend.text.position = "right",
                                                                                                legend.title.position = "top", 
                                                                                                legend.title = element_text(hjust = 0.5, size = 18, face = "bold"),
                                                                                                legend.text = element_text(size = 16),
                                                                                                panel.border = element_rect(fill = "transparent", color = "grey")) + NoAxes(),
                                    legend.show = T, min.color = "grey", max.color = "red",
                                    color.panel = c("Aneuploid" = "#FFB900", "Diploid" = "#5773CC", "Not defined" = "green")
)



# umap plot showing Cancer-Finder prediction results
p_scf_umap <-dittoSeq::dittoDimPlot(zccs905,
                                    var = c("SCF_predict2"), main= NULL, legend.title = "Cancer-Finder Prediction",
                                    reduction = "umap.harmony", theme = theme_cowplot() + theme(legend.position = "bottom", 
                                                                                                legend.justification = "centre", 
                                                                                                legend.text.position = "right",
                                                                                                legend.title.position = "top", 
                                                                                                legend.title = element_text(hjust = 0.5, size = 18, face = "bold"),
                                                                                                legend.text = element_text(size = 16),
                                                                                                panel.border = element_rect(fill = "transparent", color = "grey")) + NoAxes(),
                                    legend.show = T, min.color = "grey", max.color = "red",
                                    color.panel = c("Tumour" = "#FFB900", "Non-malignant" = "#5773CC")
)



# umap plot showing gene expressions
p_foxr2_umap <-dittoSeq::dittoDimPlot(zccs905,
                                      var = c("FOXR2"), main= NULL, legend.title = "FOXR2",
                                      reduction = "umap.harmony", theme = theme_cowplot() + theme(legend.position = "bottom", 
                                                                                                  legend.justification = "centre", 
                                                                                                  legend.title.position = "top", 
                                                                                                  legend.title = element_text(hjust = 0.5, size = 18, face = "bold"),
                                                                                                  legend.text = element_text(size = 16),
                                                                                                  panel.border = element_rect(fill = "transparent", color = "grey")) + NoAxes(),
                                      legend.show = T, min.color = "grey", max.color = "red"
)




p_myc_umap <-dittoSeq::dittoDimPlot(zccs905,
                                    var = c("MYC"), main= NULL, legend.title = "MYC",
                                    reduction = "umap.harmony", theme = theme_cowplot() + theme(legend.position = "bottom", 
                                                                                                legend.justification = "centre", 
                                                                                                legend.title.position = "top", 
                                                                                                legend.title = element_text(hjust = 0.5, size = 18, face = "bold"),
                                                                                                legend.text = element_text(size = 16),
                                                                                                panel.border = element_rect(fill = "transparent", color = "grey")) + NoAxes(),
                                    legend.show = T, min.color = "grey", max.color = "red"
)



p_mycn_umap <-dittoSeq::dittoDimPlot(zccs905,
                                    var = c("MYCN"), main= NULL, legend.title = "MYCN",
                                    reduction = "umap.harmony", theme = theme_cowplot() + theme(legend.position = "bottom", 
                                                                                                legend.justification = "centre", 
                                                                                                legend.title.position = "top", 
                                                                                                legend.title = element_text(hjust = 0.5, size = 18, face = "bold"),
                                                                                                legend.text = element_text(size = 16),
                                                                                                panel.border = element_rect(fill = "transparent", color = "grey")) + NoAxes(),
                                    legend.show = T, min.color = "grey", max.color = "red"
)




p_ptprc_umap <-dittoSeq::dittoDimPlot(zccs905,
                                    var = c("PTPRC"), main= NULL, legend.title = "PTPRC",
                                    reduction = "umap.harmony", theme = theme_cowplot() + theme(legend.position = "bottom", 
                                                                                                legend.justification = "centre", 
                                                                                                legend.title.position = "top", 
                                                                                                legend.title = element_text(hjust = 0.5, size = 18, face = "bold"),
                                                                                                legend.text = element_text(size = 16),
                                                                                                panel.border = element_rect(fill = "transparent", color = "grey")) + NoAxes(),
                                    legend.show = T, min.color = "grey", max.color = "red"
)



# arrange gene expression plots
(p_cnv_umap + p_scf_umap) / (p_ptprc_umap + p_foxr2_umap) / (p_myc_umap + p_mycn_umap)
p_all <- (p_cnv_umap + p_scf_umap) / (p_ptprc_umap + p_foxr2_umap) / (p_myc_umap + p_mycn_umap)



# save gene expression plots
ggsave(filename = "out/5_Figures/Figure_2A_and_B.svg", plot = p_all, device = "svg", width = 8, height = 14)
ggsave(filename = "out/5_Figures/Figure_2A_and_B.pdf", plot = p_all, device = "pdf", width = 8, height = 14)








