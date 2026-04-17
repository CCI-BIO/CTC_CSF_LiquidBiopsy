#CNS CTC CSF Paper
# save sample summary

library(Seurat)
library(tibble)
library(dplyr)

# get tumour cell number of each sample
stat_ctc <- as.data.frame.matrix(table(sobj_merge_s$SampleID, sobj_merge_s$cell_type_fine1)) %>% 
            as.data.frame() %>%
            tibble::rownames_to_column(var = "SampleID") %>%
            mutate(Count = `Non-malignant`+Tumour)


# get wetlab and clinical info
stat_timepoint_anno <- sobj_merge_s@meta.data %>% 
                       select(SampleID, PatientID, Disease_long, Disease_short,              
                              bioid,public_subject_id,zcc_id,                 
                              lm_enrollmnet,lm_mortality,lm_biomat_id,            
                              lm_biomat_collection,
                              V.mL.,RBC,Enrichment,               
                              Preservation,Cytology) %>%
                       dplyr::filter(!duplicated(SampleID))

# merge data frames
stat_out <- stat_ctc %>% left_join(stat_timepoint_anno, by = "SampleID")


# save table as csv
write.csv(x = stat_out, file = "out/6_Table/Summary_table.csv")
