# CNS CTC CSF Paper
# Figure 2D and Figure 2E

library(Seurat)
library(cowplot)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggpubr)
library(ggrepel)

set.seed(99)
setwd("~/workspace/Nature_Comm_Use_ZCCID")

# (if sobj not loaded)
# load the merge_2.rds from 3_Merge output folder
# sobj_merge <- readRDS("out/3_Merge/merge_2.rds")

# subset samples
zccs905 <- subset(sobj_merge_s, subset = zcc_id == "zccs905")

# Extra infomation of this patient
enroll <- unique(zccs905$lm_enrollmnet)

treatments <- read.csv("misc/zccs905_Treatments.csv", header = T) %>%
              select(-Treatment.method.specific) %>%
              mutate(enrollment = as.Date(enroll, format = "%d/%m/%Y"), 
                     Start.date = as.Date(treatment.start, format = "%d/%m/%Y"), 
                     End.date = as.Date(treatment.end, format = "%d/%m/%Y"),
                     Days.from.enroll = Start.date - enrollment,
                     Days.from.enroll.num = as.numeric(Days.from.enroll),
                     Duration = End.date - Start.date,
                     Duration.num = as.numeric(Duration)
                     )

tumour_eval <- read.csv("misc/zccs905_Tumour_Evaluation.csv", header = T) %>% 
               mutate(Tumour.evaluation.simplified = c("Diagnosis", "Diagnosis", "Progression", "Progression", "Progression"),
                      enrollment = as.Date(enroll, format = "%d/%m/%Y"),
                      Date.of.evaluation = as.Date(Date, format = "%d/%m/%Y"),
                      Days.from.enroll = Date.of.evaluation - enrollment,
                      Days.from.enroll.num = as.numeric(Days.from.enroll),
                      seg_end_point = c(Days.from.enroll.num[2:length(Days.from.enroll.num)], 736) # same date as the last recorded event
                      )
tumour_eval$Modality.of.evaluation[1] <- "MRI"

# create table for plotting
df_zccs905_meta <- zccs905@meta.data %>%
                   dplyr::select(SampleID, 
                                 zcc_id,
                                 bioid,
                                 lm_biomat_collection, 
                                 V.mL., 
                                 RBC, 
                                 Enrichment, 
                                 Preservation, 
                                 Cytology) %>%
                   dplyr::filter(!duplicated(lm_biomat_collection)) %>%
                   dplyr::arrange(bioid) %>%
                   dplyr::mutate(lm_biomat_collection = as.Date(lm_biomat_collection, format = "%d/%m/%Y"))




df_zccs905_join <- as.data.frame.matrix(table(zccs905$SampleID, zccs905$cell_type_fine1)) %>%
                   mutate(Sum_count = rowSums(.)) %>%
                   tibble::rownames_to_column(var = "Sample") %>%
                   tidyr::separate(Sample, c("pid", "bioid"), "_") %>%
                   select(-pid) %>%
                   mutate(enrollment = as.Date(enroll, format = "%d/%m/%Y")) %>% 
                   left_join(y = df_zccs905_meta, by = "bioid") %>%
                   mutate(Days.from.enroll = lm_biomat_collection - enrollment,
                          Days.from.enroll.num = as.numeric(Days.from.enroll)
                         )



df_zccs905_join_long <- tidyr::pivot_longer(data = df_zccs905_join, cols = c("Non-malignant", "Tumour"), names_to = "Cell group", values_to = "Count") %>%
                      mutate(Count_log10 = log10(Count))




# make plot -------- CTC count (line plot, clinical plots)
# setting up plot spacers
CTC_track_ymax <- max(df_zccs905_join_long %>% dplyr::filter(`Cell group` == "Tumour") %>% pull(Count_log10))
Cytology_track_ymiddle <- (CTC_track_ymax/2) * 3 * 0.95 # use the last number to adjust the track posistion
Cytology_track_y_width <- CTC_track_ymax * 0.3 # manuall setting
Tumour_eval_track_ymiddle <- (CTC_track_ymax/2) * 5 * 0.75 # use the last number to adjust the track posistion
Tumour_eval_track_y_width <- CTC_track_ymax * 0.45 # manual setting
treatment_track_ymiddle <- (CTC_track_ymax/2) * 8 * 0.75 # use the last number to adjust the track posistion
treatment_track_y_TreatmentSection_width <- CTC_track_ymax * 1.5



# track 1 -- CTC
hline_label <- c(1, 10, 100, 1000, 10000)
hline_loc_log10 <- log10(hline_label)

p1 <- ggplot(df_zccs905_join_long, aes(x = Days.from.enroll.num, y = Count_log10)) + 
     geom_point(size = 5, aes(color = `Cell group`)
                ) + 
     geom_line(linewidth = 1, alpha = 0.5, aes(color = `Cell group`), 
               linetype = "dashed") +
     geom_vline(xintercept = df_zccs905_join$Days.from.enroll.num, color = "grey50", size = 0.5, na.rm = T, linetype = 3, alpha = 0.7) +
     geom_hline(yintercept = hline_loc_log10, linetype = 3, alpha = 0.7, col = "grey50", size = 0.5) +  ## depends on max CTC of all time points
     geom_text(data = data.frame(xloc =  rep(Inf, length(hline_loc_log10)), 
                                 yloc = hline_loc_log10, 
                                 label = hline_label), 
               aes(x = xloc, y = yloc, label =label), check_overlap = T, hjust = 1, vjust = -.5) + 
     geom_label_repel(aes(x = Days.from.enroll.num, y = Count_log10, label = Count)) +
     theme_cowplot(line_size = .5)+ 
     #labs(fill = "Tumour Evaluation") +
     scale_y_log10()+
     scale_color_manual(values = c("Non-malignant" = "#5773CC", "Tumour" = "#FFB900"), name = "CSF Cell Count")





# track 2 -- add cytology data
df_zccs905_join$Cytology_rename <- NA
df_zccs905_join$Cytology_rename[df_zccs905_join$Cytology == "+"] <- "Positive"
df_zccs905_join$Cytology_rename[df_zccs905_join$Cytology == "-"] <- "Negative"
df_zccs905_join$Cytology_rename[df_zccs905_join$Cytology == ""] <- "Not Performed"

p2 <- p1 + geom_point(data = df_zccs905_join, aes(x = Days.from.enroll.num, 
                                                y = Cytology_track_ymiddle, 
                                                shape = Cytology_rename), size = 3.5, stroke = 1.5) +
           scale_shape_manual(values = c("Not Performed" = 4, "Negative" = 25, "Positive" = 21), name = "Cytology") 





# track 3 -- add tumour evaluation data
p3 <- p2 + geom_rect(data = tumour_eval, aes(xmin = Days.from.enroll.num, 
                                             xmax = seg_end_point,
                                             ymin = rep(Tumour_eval_track_ymiddle - Tumour_eval_track_y_width/2, dim(tumour_eval)[1]),
                                             ymax = rep(Tumour_eval_track_ymiddle + Tumour_eval_track_y_width/2, dim(tumour_eval)[1]),
                                             fill = Tumour.evaluation.simplified, colour = "grey99",
                                            ), inherit.aes = F, alpha = 0.7
                     ) +
           geom_label_repel(data = tumour_eval, 
                            aes(x = Days.from.enroll.num, 
                                y = Tumour_eval_track_ymiddle + Tumour_eval_track_y_width/2, 
                                label = Modality.of.evaluation), 
                      inherit.aes = F) +
          scale_fill_manual(values = c("Diagnosis" = "grey", "Progression" = "#FF9932", "Progression" = "#FF9932", "Progression" = "#FF9932", "Progression" = "#D43F3A"), 
                            breaks = c("Diagnosis", "Progression", "Progression", "Progression", "Progression"), name = "Tumour\nEvaluation")





# track 4 -- add treatment layer
### Surgery treatment is one-off, for plotting purpose, treatment duration made as 1
treatments$Duration.num[which(is.na(treatments$Duration.num))] <- 1
NumOfTreatment <- dim(treatments)[1]
treatmentSpacer <- treatment_track_y_TreatmentSection_width/NumOfTreatment
treatmentSingleTrackWidth <- treatmentSpacer*0.9
treatmentRectColor <- rep("#FFB2B2", NumOfTreatment)
names(treatmentRectColor) <- treatments$Treatment.type


if (NumOfTreatment%%2 == 1) {
  y_min <- (-(NumOfTreatment%/%2):(NumOfTreatment%/%2)) * treatmentSpacer - treatmentSingleTrackWidth/2 + treatment_track_ymiddle
  y_max <- (-(NumOfTreatment%/%2):(NumOfTreatment%/%2)) * treatmentSpacer + treatmentSingleTrackWidth/2 + treatment_track_ymiddle
  y_middle <- (-(NumOfTreatment%/%2):(NumOfTreatment%/%2)) * treatmentSpacer + treatment_track_ymiddle
}else{
  y_min <- (-((NumOfTreatment%/%2) - 1):(NumOfTreatment%/%2)) * treatmentSpacer - treatmentSingleTrackWidth/2 + treatment_track_ymiddle
  y_max <- (-((NumOfTreatment%/%2) - 1): (NumOfTreatment%/%2)) * treatmentSpacer + treatmentSingleTrackWidth/2 + treatment_track_ymiddle
  y_middle <- (-((NumOfTreatment%/%2) - 1): (NumOfTreatment%/%2)) * treatmentSpacer + treatment_track_ymiddle
}

TreatmentBoarderColor <- rep("grey", NumOfTreatment)
TreatmentBoarderColor[treatments$Treatment.type == "Surgery"] <- "#4B1F6F"

p4 <- p3 + annotate("rect", xmin = treatments$Days.from.enroll.num, 
                            xmax = treatments$Duration.num + treatments$Days.from.enroll.num, 
                            ymin = y_min, 
                            ymax = y_max, 
                            alpha = 0.7, 
                            fill = "#FFB2B2", colour = TreatmentBoarderColor) +
           annotate("text", x = treatments$Days.from.enroll.num, 
                            y = y_middle, 
                            label = treatments$Treatment.type, 
                            vjust = 0.5, hjust=0)




# plot settings
y_breaks <- c(CTC_track_ymax/2, Cytology_track_ymiddle, Tumour_eval_track_ymiddle, treatment_track_ymiddle)
y_label <- c("CSF cell count", "Cytology diagnosis", "Tumour evaluation", "Treatment")
p5 <- p4 + scale_y_continuous(breaks = y_breaks, labels = y_label) + 
     xlab("Days from Diagnosis") + 
     scale_x_continuous(position = "top", n.breaks = 10, limits = c(0, 850)) + 
     theme(axis.title.y = element_blank(), axis.line = element_blank(), 
           axis.title.x = element_text(face ="bold"), legend.title = element_text(face = "bold"),
           panel.border = element_rect(colour = "black", linewidth = 0.1) #, panel.grid.major.x = element_line(colour = "grey", linewidth = 0.5, linetype = "dotted")
           )

p5


ggsave(filename = "out/5_Figures/Figure_2D.svg", plot = p5, device = "svg", width = 11, height = 6.5)
ggsave(filename = "out/5_Figures/Figure_2D.pdf", plot = p5, device = "pdf", width = 11, height = 6.5)






# make plot (CTC UMAP plots)
temp_time_name <- unique(zccs905$lm_biomat_collection)
temp_time <- as.Date(temp_time_name, format = "%d/%m/%Y") - as.Date(enroll, format = "%d/%m/%Y")
names(temp_time) <- temp_time_name
temp_time <- temp_time[order(temp_time)]

Idents(zccs905) <- "lm_biomat_collection"
zccs905 <- RenameIdents(zccs905, temp_time)
zccs905[["Days"]] <- Idents(zccs905)

p6 <- dittoSeq::dittoDimPlot(zccs905, var = "cell_type_fine1", reduction.use = "umap.harmony", 
                                               legend.title = "Cell group",
                                               color.panel = c("Non-malignant" = "#5773CC", "Tumour" = "#FFB900"), split.by = "Days", split.show.all.others = T, 
                                               split.ncol = 7, theme = theme_cowplot(), do.contour = F, xlab = "Days From Diagnosis") & scale_y_continuous(breaks = 0) &
                        ggplot2::theme(strip.background =element_blank(), 
                        strip.text = element_text(size = 12),
                        axis.line = element_blank(), 
                        axis.ticks = element_blank(), 
                        axis.text = element_blank(), 
                        axis.title = element_text(size = 14, face = "bold"),
                        plot.title = element_blank(),
                        legend.position = "right", 
                        legend.title = element_text(size = 14, face = "bold"),
                        legend.text = element_text(size = 12),
                        legend.justification = "centre",
                        panel.border = element_rect(fill = "transparent", colour = "grey80", linetype = "solid"))


# add another layer to p6
dot_highlight <- as.data.frame(Embeddings(zccs905, "umap.harmony"))
dot_highlight$others <- "Other timepoints"

p6 <- p6 + geom_point(data = dot_highlight, 
                      aes(x=umapharmony_1, y = umapharmony_2, color = others),
                      inherit.aes = F) &
         scale_color_manual(values = c("Non-malignant" = "#5773CC", 
                                       "Tumour" = "#FFB900", 
                                       "Other timepoints" = "grey90"),
                            breaks = c("Non-malignant", "Tumour", "Other timepoints"),
                            name = "Cell group")

# reorder the new layer to the back
p6$layers <- c(p6[["layers"]][[3]], p6[["layers"]][[1]], p6[["layers"]][[2]])

p6[["facet"]][["params"]][["strip.position"]] <- "bottom"
p7 <- p6 + theme(axis.title.y = element_blank())





# save image
ggsave(filename = "out/5_Figures/Figure_2E.svg", plot = p7, device = "svg", width = 15, height = 2.5)
ggsave(filename = "out/5_Figures/Figure_2E.pdf", plot = p7, device = "pdf", width = 15, height = 2.5)

# combine plots
ggpubr::ggarrange(p5, p7, ncol = 1, heights = c(3.2,1), legend = "right")
ggsave(filename = "out/5_Figures/Figure_2D_and_E.svg", device = "svg", width = 15, height = 10)
ggsave(filename = "out/5_Figures/Figure_2D_and_E.pdf", device = "pdf", width = 15, height = 10)
