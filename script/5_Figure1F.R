# CNS CTC CSF Paper
# Figure 1F

library(Seurat)
library(dittoSeq)
library(scCustomize)
library(SCpubr)
library(dplyr)
library(tidyseurat)
library(tidyHeatmap)
library(tidyheatmaps)
library(SeuratExtend)
library(readxl)
library(patchwork)
library(ggbump)
library(cowplot)

set.seed(99)
setwd("~/workspace/Nature_Comm_Use_ZCCID")

# load Seurat object from 4_Cell_Type_Curation
# sobj_merge_s <- readRDS("out/4_Cell_Type_Curation/curated_1.rds")

# rename idents to get Tumour and non-Tumour annotation
Idents(sobj_merge_s) <- "cell_type_fine1"
sobj_merge_s <- RenameIdents(sobj_merge_s, "Non-malignant" = "Other")
sobj_merge_s[["cell_type_fine2"]] <- Idents(sobj_merge_s)

# make plot
y1 <- SCpubr::do_BarPlot(sobj_merge_s, 
                         group.by = "cell_type_fine2",
                         split.by = "Disease_short", 
                         xlab = "Disease", 
                         legend.title = "Group",
                         position = "fill",
                         add.n = T,
                         return_data = T,
                         flip = FALSE)

label_df <- y1$Data %>% dplyr::filter(cell_type_fine2 == "Tumour") %>%
            rowwise() %>%
            mutate(label = paste("Tumour = ", n, " (", round(freq*100,2),"%)", sep = ""))

ggplot(y1$Data, aes(x = Disease_short, y = n, fill = cell_type_fine2)) +
  geom_bar(stat = "identity", position = "fill") +
  labs(x = "Disease", y = "Frequency", fill = "Group") +
  scale_x_discrete(guide = guide_axis(angle = 0)) +  
  scale_fill_manual(values = c("Tumour" = "#FFB900", "Other" = "#5773CC"), name = "Cells") +
  labs(color = "Cells")+
  coord_flip() +
  cowplot::theme_cowplot() + 
  theme(legend.position = "top", legend.justification = "centre", plot.margin = margin(r=0.5, t = 0.1, b = 0.1, l = 0.1, unit = "in"),
        axis.title = element_text(face = "bold", size = 16), legend.title = element_text(face = "bold", size = 16), 
        legend.text = element_text(size = 16))

ggsave(filename = "out/5_Figures/Figure_1F.svg", device = "svg", width = 6, height = 4.5)
ggsave(filename = "out/5_Figures/Figure_1F.pdf", device = "pdf", width = 6, height = 4.5)
