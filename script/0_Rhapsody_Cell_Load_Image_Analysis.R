# BD Rhapsody Cell Load Image Analysis
# Author: Wenyan Li
# Contact: wli@ccia.org.au

# / load libraries
library(magick)
library(EBImage)
library(tidyverse)
library(ggplot2)
library(ggExtra)
library(patchwork)
library(ggsankey)

# /// ------------ User defined parameters ------------ 
#Set up working directory
cur_path_image <- getwd()

# These numbers mean [0 to 8um], [8um - #um] and [#um to Infinity]. 
# The cells will be divided into 3 groups based on sizes.
cell_size_cutpoint <- c(0, 8, 20, Inf)

# parameters to remove images with fluorescence noise only
max_image_intensity <- 0.3
mean_image_intensity <- 0.2
fluo_region_pct <- 0.5

# cluster images in how many communities
n_clusters <- 6

# set image crop size to 30px by 30px
crop_size <- c(30, 30) 

# set image zoom in to 3x (for GIF output only)
image_zoom_times <- 3

# /// ------------ Create output folder ------------ 
ifelse(!dir.exists(file.path(cur_path_image, "Image Output")), 
        dir.create(file.path(cur_path_image, "Image Output")), FALSE) #stop("Image Output folder is already exist!"))

# /// ------------ Image folder ------------ 
subfolder_name <- "Cell Load"

# /// ------------ Format Parameters ------------ 
crop_size_3_times_zoom <- crop_size*image_zoom_times
crop_size_3_times_zoom <- paste(crop_size_3_times_zoom, collapse = 'x')

crop_shift <- floor(crop_size[1]/2)
crop_size <- paste(crop_size, collapse = "x")

cell_group_labels <- sapply(cell_size_cutpoint, 
                            FUN = function(y) 
                              ifelse(y != cell_size_cutpoint[length(cell_size_cutpoint)],
                                     paste0(y, "-", cell_size_cutpoint[grep(y,cell_size_cutpoint)+1],
                                            " um"), y))[-length(cell_size_cutpoint)]

folder_name <- strsplit(cur_path_image, split = "/") %>% 
               unlist() %>% 
               .[length(.)]

# /// ------------ Gather image information ------------ 
input_info <- list(
                image.files.path = "./Cell Load",
                csv.files.path = "./Cell Load/IA_Result",
                BF.image.files = list.files("./Cell Load", pattern = "*.BF.png"),
                FG.image.files = list.files("./Cell Load", pattern = "*.FG.png"),
                csv.files = list.files("./Cell Load/IA_Result", pattern = "*.BF.csv"),
                file.length = length(list.files("./Cell Load/IA_Result", pattern = "*.BF.csv")),
                subfolder = "Cell Load",
                crop_size = crop_size,
                crop_shift = crop_shift,
                zoom = crop_size_3_times_zoom,
                cuts = cell_size_cutpoint,
                cut_labels = cell_group_labels,
                fl_max = max_image_intensity,
                fl_mean = mean_image_intensity,
                fl_ratio = fluo_region_pct,
                n_clusters = n_clusters,
                folder_name = folder_name
              )

# /// ------------ Build functions ------------ 

######## curate cell position tables
######## support function for func_build_image_crops
func_build_image_df <- function(file_path){
  
  x.csv <- read.csv(file_path,
                    skip = 49,
                    col.names = c(1:50),
                    na.strings = "",
                    header = F)

  # remove empty columns
  col.with.content <- x.csv %>% 
                      apply(., MARGIN = 2, FUN = function(x) sum(!is.na(x))) > 0
  x.csv <- x.csv[,col.with.content]
  
  # mark section rows
  well.content.row <- match("Well Col", x.csv[,1])
  bead.content.row <- match("Bead Diameter(um)", x.csv[,2])
  cell.content.row <- match("Cell Diameter(um)",x.csv[,2])

  # / well content data frame
  well.df <- x.csv[(well.content.row+1) : (bead.content.row-1), ] %>% 
             dplyr::filter(!is.na(X6))
  
  # / multiplet IDs
  if (dim(well.df)[2] > 6) {
    well.df.multi <- well.df %>% dplyr::filter(!is.na(X7))
    well.df.single <- well.df %>% dplyr::filter(is.na(X7))
  
    multi.ID <- well.df.multi[,6:dim(well.df.multi)[2]] %>% 
                unlist() %>% 
                {.[!is.na(.)]} %>% 
                data.frame(CellID = ., Multi = TRUE)
  
    multi.to.Plot <- data.frame(CellID = well.df.multi$X6, MultiPlot = TRUE)
    single.ID <- data.frame(CellID = well.df.single$X6, Single = TRUE)
    CellInWell <- data.frame(CellID = c(multi.ID$CellID, single.ID$CellID), CellInWell = TRUE)
  
    # / cell content data frame
    cell.df <- x.csv %>% dplyr::filter(stringr::str_starts(X1, "FG#")) %>% 
      select(CellID = X1, CellSize = X2, X_centre = X3, Y_centre = X4) %>%
      left_join(CellInWell, by = "CellID") %>% 
      left_join(single.ID, by = "CellID") %>% 
      left_join(multi.ID, by = "CellID") %>% 
      left_join(multi.to.Plot, by = "CellID")
    
    
    }else{
      well.df.multi <- well.df[NULL, ]
      well.df.single <- well.df
      
      single.ID <- well.df %>% select(CellID = X6) %>% mutate(CellInWell = TRUE, Single = TRUE)
      
      cell.df <- x.csv %>% dplyr::filter(stringr::str_starts(X1, "FG#")) %>% 
                 select(CellID = X1, CellSize = X2, X_centre = X3, Y_centre = X4) %>%
                 left_join(single.ID, by = "CellID") %>% 
                 mutate(Multi = NA, MultiPlot = NA)
    }
  
  cell.df <- cell.df %>% mutate(CellSize = as.numeric(CellSize), 
                                X_centre = as.numeric(X_centre), 
                                Y_centre = as.numeric(Y_centre),
                                data_source = file_path)
  
  well.df.multi <- well.df.multi %>% 
                   mutate(X1 = as.numeric(X1), 
                          X2 = as.numeric(X2)) %>% 
                   select(c(X1,X2))
  
  well.df.single <- well.df.single %>% 
                    mutate(X1 = as.numeric(X1), 
                           X2 = as.numeric(X2)) %>% 
                    select(c(X1,X2))

  # / return data frame
  return(list(cell.df = cell.df,
              well_single_df = well.df.single,
              well_multi_df = well.df.multi))
}


######## crop images using cell locations
######## support function for func_build_image_crops
func_get_cropped_image_vector <- function(df,
                                          image,
                                          crop_region,
                                          shift){
  if (dim(df)[1] > 0) {

    # Create a dummy image
    image_vector <- image_crop(image, crop_region)
    
    # Loop to find cells based on df cell content information
    for (j in 1:length(df$X_centre)) {
    
      geometry <- paste(crop_region, 
                        as.numeric(df$X_centre[j])-shift,
                        as.numeric(df$Y_centre[j])-shift, 
                        sep = "+") 
    
      cropped <- image_crop(image, geometry)
      image_vector <- c(image_vector, cropped)
    }
  
    # get rid of the first dummy image
    image_vector <- image_vector[-1]
    
  }else{
    image_vector <- c()
  }
  
  # return image vector
  return(image_vector)
}


######## image threshold function
######## support function for func_build_image_crops

func_filter_low_FI_image <- function(images,
                                     images_BF,
                                     image_size,
                                     ratio_thrd,
                                     n_clusters,
                                     FI_max_thrd, 
                                     FI_mean_thrd){
  
  func_area_ratio <- function(image_frame, thred){
    x <- image_frame > thred
    out <- sum(x)/(dim(x)[1]*dim(x)[2])
    return(out)
  }
  
  if (length(images) > 0) {
    
    cell_area_ratio <- seq_len(length(images))
    cell_FI_max <- seq_len(length(images))
    cell_FI_meam <- seq_len(length(images))
    image_index <- seq_len(length(images))
    
    for (i in seq_along(images)) {
    image_x <- as_EBImage(images[i])
    cell_area_ratio[i] <- func_area_ratio(image_frame = image_x, thred = otsu(image_x))
    cell_FI_max[i] <- max(image_x)
    cell_FI_meam[i] <- mean(image_x)
    }
    
    ## image clustering
    if (length(images) > 10) {

      clusters <- 1:length(images)
      ebi_images <- lapply(images, as_EBImage )
      images_full_crop_area <- sapply(ebi_images, length) == image_size
      ebi_images <- ebi_images[images_full_crop_area]
      images_df <- sapply(ebi_images, as.array) %>% t()
    
      pca <- prcomp(images_df, center = TRUE, scale = TRUE)
      temp_pca <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2])
    
      set.seed(99)
      cluster_pca <- kmeans(pca$x[, 1:10], n_clusters)
    
      clusters[images_full_crop_area] <- cluster_pca$cluster
      clusters[!images_full_crop_area] <- "incomplete image"
    
      cluster_images <- split(x = images, f = clusters)
      cluster_images_BF <- split(x = images_BF, f = clusters)
    
    }else{
      clusters <- rep("1", length(images))
      cluster_images <- list("1" = images)
      cluster_images_BF <- list("1" = images)
    }
    
    ## filter images
    image_index <- c(which(cell_area_ratio > ratio_thrd), 
                         which(cell_FI_max < FI_max_thrd), 
                         which(cell_FI_meam > FI_mean_thrd) ) %>% 
                   unique()
    
    ## create output data frame
    cell_call_df <- data.frame(max = cell_FI_max, 
                               mean = cell_FI_meam, 
                               ratio = cell_area_ratio,
                               cluster = clusters) %>% 
                    mutate(max_noise = max < FI_max_thrd,
                           mean_noise = mean > FI_mean_thrd,
                           ratio_noise = ratio > ratio_thrd)
    
  }else{
    cell_call_df <- data.frame(max = integer(0), 
                               mean = integer(0), 
                               ratio = integer(0),
                               cluster = integer(0),
                               max_noise = integer(0),
                               mean_noise = integer(0),
                               ratio_noise = integer(0))
    image_index <- integer(0)
    cluster_images <- integer(0)
    cluster_images_BF <- integer(0)
  }
  return(list(poor_image_index = image_index,
              cluster_images = cluster_images,
              cluster_images_BF = cluster_images_BF,
              cell_call_df = cell_call_df)
         )
}


######## create a pool of images of all cells
func_build_image_crops <- function(input_info){
  
  # / file list
  BF_layer = "BF.image.files"
  BF_image_path <- paste(input_info$image.files.path, 
                         input_info$BF.image.files, 
                         sep = "/")

  FG_layer = "FG.image.files"
  FG_image_path <- paste(input_info$image.files.path, 
                         input_info$FG.image.files, 
                         sep = "/")
  
  csv_path <- paste(input_info$csv.files.path,
                    input_info$csv.files,
                    sep = "/")
  
  # output object
  BF_single_vector <- c()
  FG_single_vector <- c()
  BF_multi_vector <- c()
  FG_multi_vector <- c()
  single_df <- c()
  multi_df <- c()
  single_well_pos_df <- c()
  multi_well_pos_df <- c()
  
  # / Get cropped images
  for (i in 1:input_info$file.length) {
    
    # / create cell location data frame
    all_df <- func_build_image_df(csv_path[i])
    cell_df <- all_df$cell.df
    
    single_cell <- cell_df %>% 
                   dplyr::filter(Single==TRUE & CellInWell==TRUE)
    
    multi_cell <- cell_df %>% 
                  dplyr::filter(MultiPlot==TRUE & CellInWell==TRUE)
    
    # / read image
    BF.image.load <- image_read(BF_image_path[i])
    FG.image.load <- image_read(FG_image_path[i])
    
    # get cropped image vector
    BF_temp_single <- func_get_cropped_image_vector(df = single_cell, 
                                                    image = BF.image.load,
                                                    crop_region = input_info$crop_size,
                                                    shift = input_info$crop_shift)
    
    FG_temp_single <- func_get_cropped_image_vector(df = single_cell, 
                                                    image = FG.image.load,
                                                    crop_region = input_info$crop_size,
                                                    shift = input_info$crop_shift)
    
    BF_temp_multi <- func_get_cropped_image_vector(df = multi_cell, 
                                                   image = BF.image.load,
                                                   crop_region = input_info$crop_size,
                                                   shift = input_info$crop_shift)
    
    FG_temp_multi <- func_get_cropped_image_vector(df = multi_cell, 
                                                   image = FG.image.load,
                                                   crop_region = input_info$crop_size,
                                                   shift = input_info$crop_shift)
    
    BF_single_vector <- append(BF_single_vector, BF_temp_single)
    FG_single_vector <- append(FG_single_vector, FG_temp_single)
    
    BF_multi_vector <- append(BF_multi_vector, BF_temp_multi)
    FG_multi_vector <- append(FG_multi_vector, FG_temp_multi)
    
    # get cell sizes
    single_df <- rbind(single_df, single_cell)
    multi_df <- rbind(multi_df, multi_cell)
    
    # get singlets and multiplets positions in the cartridge
    single_well_pos_df <- rbind(single_well_pos_df, all_df$well_single_df)
    multi_well_pos_df <- rbind(multi_well_pos_df, all_df$well_multi_df)
    
  }
  
  # / order images based on sizes (ascending)
  temp_single_order <- order(single_df$CellSize)
  BF_single_vector <- BF_single_vector[temp_single_order]
  FG_single_vector <- FG_single_vector[temp_single_order]
  single_df <- single_df[temp_single_order, ]
  single_well_pos_df <- single_well_pos_df[temp_single_order, ]
  
  temp_multi_order <- order(multi_df$CellSize)
  BF_multi_vector <- BF_multi_vector[temp_multi_order]
  FG_multi_vector <- FG_multi_vector[temp_multi_order]
  multi_df <- multi_df[temp_multi_order, ]
  multi_well_pos_df <- multi_well_pos_df[temp_multi_order, ]
  
  # / save this for manual images filtering using clusters number(s)
  cell_and_noise_images <- list(BF_single_and_noise_vector = BF_single_vector,
                                FG_single_and_noise_vector = FG_single_vector,
                                BF_multi_and_noise_vector = BF_multi_vector,
                                FG_multi_and_noise_vector = FG_multi_vector)
  
  # / Image QC filtering (using FG images)
  image_crop_dim <- input_info$crop_size %>% strsplit("x") %>% unlist() %>% as.numeric()
  image_crop_region <- image_crop_dim[1]*image_crop_dim[2]
  poor_single_cell_call <- func_filter_low_FI_image(FG_single_vector,
                                                    BF_single_vector,
                                                    image_size = image_crop_region,
                                                    ratio_thrd = input_info$fl_ratio,
                                                    n_clusters = input_info$n_clusters,
                                                    FI_max_thrd = input_info$fl_max,
                                                    FI_mean_thrd = input_info$fl_mean)

  poor_multi_cell_call <- func_filter_low_FI_image(FG_multi_vector,
                                                   BF_multi_vector,
                                                   image_size = image_crop_region,
                                                   ratio_thrd = input_info$fl_ratio,
                                                   n_clusters = input_info$n_clusters,
                                                   FI_max_thrd = input_info$fl_max,
                                                   FI_mean_thrd = input_info$fl_mean)
  
  poor_single_images_index <- poor_single_cell_call$poor_image_index
  poor_multi_images_index <- poor_multi_cell_call$poor_image_index
  
  BF_poor_single_vector <- BF_single_vector[poor_single_images_index]
  FG_poor_single_vector <- FG_single_vector[poor_single_images_index]
  
  BF_poor_multi_vector <- BF_multi_vector[poor_multi_images_index]
  FG_poor_multi_vector <- FG_multi_vector[poor_multi_images_index]
  
  BF_single_vector <- BF_single_vector[setdiff(seq_along(BF_single_vector), poor_single_images_index)]
  FG_single_vector <- FG_single_vector[setdiff(seq_along(FG_single_vector), poor_single_images_index)]
  
  BF_multi_vector <- BF_multi_vector[setdiff(seq_along(BF_multi_vector), poor_multi_images_index)]
  FG_multi_vector <- FG_multi_vector[setdiff(seq_along(FG_multi_vector), poor_multi_images_index)]
  
  single_df$FalseCall <- ifelse(seq_len(dim(single_df)[1]) %in% poor_single_images_index, TRUE, FALSE)
  multi_df$FalseCall <- ifelse(seq_len(dim(multi_df)[1]) %in% poor_multi_images_index, TRUE, FALSE)

  single_well_pos_df$FalseCall <- ifelse(seq_len(dim(single_df)[1]) %in% poor_single_images_index, TRUE, FALSE)
  multi_well_pos_df$FalseCall <- ifelse(seq_len(dim(multi_df)[1]) %in% poor_multi_images_index, TRUE, FALSE)
  
  # / return output
  return(list(BF_singlets_all = BF_single_vector,
              FG_singlets_all = FG_single_vector,
              BF_multiplets_all = BF_multi_vector,
              FG_multiplets_all = FG_multi_vector,
              cell_and_noise_images = cell_and_noise_images,
              single_cell_info = single_df,
              multi_cell_info = multi_df,
              single_pos_dist = single_well_pos_df,
              multi_pos_dist = multi_well_pos_df,
              BF_singlets_noise = BF_poor_single_vector,
              FG_singlets_noise = FG_poor_single_vector,
              BF_multiplets_noise = BF_poor_multi_vector,
              FG_multiplets_noise = FG_poor_multi_vector,
              poor_single_cell_call = poor_single_cell_call,
              poor_multi_cell_call = poor_multi_cell_call)
         )
}


######## to combine a list of images into one
######## supporting function for func_make_gifs_and_combine_images
func_images_combine <- function(image_vector){
  
  image.chunks <- split(image_vector,
                        ceiling(seq_along(image_vector)/floor(sqrt(length(image_vector)))))
  
  image.chunks.stackFalse <- lapply(image.chunks,magick::image_append, 
                                    stack = F)
  
  image.whole <- image.chunks.stackFalse[[1]]
  
  if (length(image.chunks.stackFalse) >1 ) {
      for (i in 1:(length(image.chunks.stackFalse)-1)) {
        image.whole <- image_append(c(image.whole,
                                      image.chunks.stackFalse[[i+1]]),
                                    stack = T)
      }
  }
  
  return(image.whole)
} 


######## create gif and combined images in different ranges of cell sizes
func_make_gifs_and_combine_images <- function(
                    cutpoint, 
                    cutpoint_labels,
                    image_vector,
                    zoom_region,
                    toCreateGIF, 
                    cell_sizes
  ){
  
  temp.cut <- cut(cell_sizes, 
                  cutpoint, 
                  cutpoint_labels)
  
  temp <- table(temp.cut)
  temp.labels <- names(temp[temp !=0])
  temp.image.vector.list <- as.list(temp.labels)
  names(temp.image.vector.list) <- temp.labels
  
  temp_index <- 1:length(image_vector)
  temp.image.vector.list <- lapply(temp.image.vector.list, 
                                   FUN = function(x)
                                     image_vector[temp_index[temp.cut == x]])
  
  temp.image.gif.list <- list()
  temp.image.combined.list <- list()
  
  # Combine images and/or make gif
  for (i in names(temp.image.vector.list)) {
      
      temp.image.combined.list[[i]] <- func_images_combine(temp.image.vector.list[[i]]) 
    
      if (toCreateGIF) {
      temp.image.gif.list[[i]] <- image_animate(image_scale(temp.image.vector.list[[i]], 
                                                            geometry = zoom_region),
                                                fps = 4, 
                                                dispose = "previous")
      }
  }
  
  # Return processed data
  return(list(image.combined = temp.image.combined.list, 
              image.gif = temp.image.gif.list))
}


######## function to blend bright field and fluorescence images
func_image_blending <- function(BF_img, FG_img){
  converted <- FG_img %>% 
    as_EBImage(.) %>%
    EBImage::channel(.,"asgreen") %>% 
    magick::image_read(.)
  
  blended <- image_composite(BF_img, 
                             converted, 
                             operator = "blend") 
  return(blended)
}


######## function to filter cells based on cluster number(s)
func_extra_cluster_filter <- function(results, 
                                      cluster_singlet = integer(0), 
                                      cluster_multiplet = integer(0)){
  
  out <- list(singlets_filtered = c(),
              multiplets_filtered = c())
  
  # singlet
  if (length(cluster_singlet) > 0) {
    temp_index <- data.frame(Type = results$processed_images_dataframe$single_cell_info$FalseCall,
                             clusters = results$processed_images_dataframe$poor_single_cell_call$cell_call_df$cluster) %>%
      dplyr::filter(Type == FALSE & clusters %in% cluster_singlet) %>%
      {which(.$Type == FALSE & .$clusters %in% cluster_singlet)}
    
    temp_filter <- lapply(results$processed_images_dataframe$cell_and_noise_images[c(1:2)], FUN = function(x) x[temp_index])
    img_out <- lapply(temp_filter, func_images_combine)
    out$singlets_filtered <- func_image_blending(BF_img = img_out[[1]], FG_img = img_out[[2]])
  }
  
  # multiplets
  if (length(cluster_multiplet) > 0){
    temp_index <- data.frame(Type = results$processed_images_dataframe$multi_cell_info$FalseCall,
                             clusters = results$processed_images_dataframe$poor_multi_cell_call$cell_call_df$cluster) %>%
      dplyr::filter(Type == FALSE & clusters %in% cluster_multiplet) %>%
      {which(.$Type == FALSE & .$clusters %in% cluster_multiplet)}
    
    temp_filter <- lapply(results$processed_images_dataframe$cell_and_noise_images[c(3:4)], FUN = function(x) x[temp_index])
    img_out <- lapply(temp_filter, func_images_combine)
    out$multiplets_filtered <- func_image_blending(BF_img = img_out[[1]], FG_img = img_out[[2]])
  }
  
  return(out)
}


######## workflow function
func_workflow <- function(input_info){
  
  # 1. get cell position tables and crop images
  img_df <- func_build_image_crops(input_info)
  img <- img_df[sapply(img_df, class) == "magick-image" & sapply(img_df, length) > 0]
  df <- img_df[sapply(img_df, class) == "data.frame"]
  
  # 2. combine all cell images
  img_combined <- lapply(img, func_images_combine)

  # 3. blending bright field and fluorescence images
  temp_odd <- img_combined[seq(1, length(img_combined),2)]
  temp_even <- img_combined[seq(2, length(img_combined),2)]
  
  img_blended <- mapply(func_image_blending, temp_odd, temp_even, SIMPLIFY = F)
  
  names(img_combined) <- paste("combined", names(img_combined), sep = "_")
  names(img_blended) <- paste("blended", names(img_blended), sep = "_")
  
  # 4. split cells with pre-set cutoffs (single cells only)
  BF_split_single <- func_make_gifs_and_combine_images(cutpoint = input_info$cuts,
                                                       cutpoint_labels = input_info$cut_labels, 
                                                       image_vector = img_df$BF_singlets_all,
                                                       zoom_region = input_info$zoom,
                                                       toCreateGIF = TRUE,
                                                       cell_sizes = img_df$single_cell_info %>% 
                                                                    dplyr::filter(FalseCall == FALSE) %>% 
                                                                    pull(CellSize))
  
  FG_split_single <- func_make_gifs_and_combine_images(cutpoint = input_info$cuts,
                                                       cutpoint_labels = input_info$cut_labels, 
                                                       image_vector = img_df$FG_singlets_all,
                                                       zoom_region = input_info$zoom,
                                                       toCreateGIF = TRUE,
                                                       cell_sizes = img_df$single_cell_info %>% 
                                                                    dplyr::filter(FalseCall == FALSE) %>%
                                                                    pull(CellSize))
  
  # 5. save images
  
  ## a) histogram
  temp_df_s <- data.frame(CellSize = img_df$single_cell_info$CellSize, 
                          Type = ifelse(img_df$single_cell_info$FalseCall, "Noise", "Singlet"))
  
  temp_df_m <- data.frame(CellSize = img_df$multi_cell_info$CellSize, 
                          Type = ifelse(img_df$multi_cell_info$FalseCall, "Noise", "Multiplet"))
  
  temp_df <- rbind(temp_df_s, temp_df_m)
  
  p_hist <- ggplot(data = temp_df, 
                   aes(x = CellSize, fill = Type)) + 
    geom_histogram(binwidth = 1,
                   position = "dodge", 
                   # alpha =0.8,
                   colour = "black",
                   lwd = 0.75,
                   linetype = 1) +
    scale_y_log10() + 
    xlab("Cell Size (μm)") + 
    ylab("Count") + 
    theme(axis.text = element_text(size = 12)) + 
    theme(axis.title = element_text(size = 12)) +
    theme(plot.title = element_text(size = 14)) +
    ggtitle(paste("Histogram of sample -", 
                  folder_name))
  
  ggsave(filename = "Cell Size Distribution (log10_Y).png",
          dpi = 300,
          plot = p_hist, 
          path = "./Image Output")
  
  p_hist2 <- ggplot(data = temp_df, 
                   aes(x = CellSize, fill = Type)) + 
    geom_histogram(binwidth = 1,
                   position = "dodge", 
                   # alpha =0.8,
                   colour = "black",
                   lwd = 0.75,
                   linetype = 1) +
    # scale_y_log10() + 
    xlab("Cell Size (μm)") + 
    ylab("Count") + 
    theme(axis.text = element_text(size = 12)) + 
    theme(axis.title = element_text(size = 12)) +
    theme(plot.title = element_text(size = 14)) +
    ggtitle(paste("Histogram of sample -", 
                  folder_name))
  
  ggsave(filename = "Cell Size Distribution.png",
         dpi = 300,
         plot = p_hist2, 
         path = "./Image Output")
  
  ## b) images
  #### all cells 
  lapply(names(img_combined), FUN=function(x) image_write(img_combined[[x]], path = paste0("./Image Output/", x, ".png")))
  lapply(names(img_blended), FUN=function(x) image_write(img_blended[[x]], path = paste0("./Image Output/", x, ".png")))
  
  #### split cells
  temp_name <- names(BF_split_single$image.combined)
  if (!is.null(temp_name)) {
    lapply(temp_name, FUN = function(x) image_write(BF_split_single$image.combined[[x]],
                                                    path = paste("./Image Output", paste0(x, "_singlet_BF.png"), 
                                                                 sep = "/")))
  }
  
  temp_name <- names(BF_split_single$image.gif)
  if (!is.null(temp_name)) {
    lapply(temp_name, FUN = function(x) image_write(BF_split_single$image.gif[[x]],
                                                    path = paste("./Image Output", paste0(x, "_singlet_BF.gif"), 
                                                                 sep = "/")))
  }
  
  temp_name <- names(FG_split_single$image.combined)
  if (!is.null(temp_name)) {
    lapply(temp_name, FUN = function(x) image_write(FG_split_single$image.combined[[x]],
                                                    path = paste("./Image Output", paste0(x, "_singlet_FG.png"), 
                                                                 sep = "/")))
  }
  
  temp_name <- names(FG_split_single$image.gif)
  if (!is.null(temp_name)) {
    lapply(temp_name, FUN = function(x) image_write(FG_split_single$image.gif[[x]],
                                                    path = paste("./Image Output", paste0(x, "_singlet_FG.gif"), 
                                                                 sep = "/")))
  }
  
  #### cluster cells
  # check whether clustering step is skipped
  cluster_number <- length(img_df$poor_single_cell_call$cluster_images)
  if (cluster_number > 1) {
    each_cluster_combined <- lapply(img_df$poor_single_cell_call$cluster_images, func_images_combine)
    each_cluster_combined_BF <- lapply(img_df$poor_single_cell_call$cluster_images_BF, func_images_combine)
    img_blended <- mapply(func_image_blending, each_cluster_combined_BF, each_cluster_combined, SIMPLIFY = F)
    temp_cluster_name <- names(img_blended)
    lapply(temp_cluster_name, FUN = function(x) image_write(img_blended[[x]],
                                                    path = paste("./Image Output", paste0("cluster_", x, "_singlet_FG.png"), 
                                                                 sep = "/")))
  }
  
  cluster_number <- length(img_df$poor_multi_cell_call$cluster_images)
  if (cluster_number > 1) {
    each_cluster_combined <- lapply(img_df$poor_multi_cell_call$cluster_images, func_images_combine)
    each_cluster_combined_BF <- lapply(img_df$poor_multi_cell_call$cluster_images_BF, func_images_combine)
    img_blended <- mapply(func_image_blending, each_cluster_combined_BF, each_cluster_combined, SIMPLIFY = F)
    temp_cluster_name <- names(img_blended)
    lapply(temp_cluster_name, FUN = function(x) image_write(img_blended[[x]],
                                                            path = paste("./Image Output", paste0("cluster_", x, "_multiplet_FG.png"), 
                                                                         sep = "/")))
  }
  
  ## c) cell positions in Cartridge
  pos_single <- img_df$single_pos_dist %>% 
                mutate(Type = ifelse(img_df$single_cell_info$FalseCall, "Noise", "Singlet"))
  
  pos_multi <- img_df$multi_pos_dist %>% 
               mutate(Type = ifelse(img_df$multi_cell_info$FalseCall, "Noise", "Multiplet"))
  
  if (dim(pos_multi)[1] > 10) {
    pos <- rbind(pos_single, pos_multi)
  }else{
    pos <- pos_single
  }
  
  
  p_pos <- ggplot(pos, aes(x=X1, y=X2, colour=Type)) +
           geom_point(shape = 20) + 
           geom_density_2d() + 
           geom_rug(position="identity", linewidth = 0.1) + 
           theme_classic() +
           xlab("Well Position (X)") + 
           ylab("Well Position (Y)") + 
           theme(axis.text = element_text(size = 12)) + 
           theme(axis.title = element_text(size = 12)) +
           theme(plot.title = element_text(size = 14)) +
           theme(legend.position = "bottom") + 
           ggtitle(paste("Cell positions in cartridge (contour) -", 
                   folder_name))
  
  p_pos_out <- ggMarginal(p_pos, type='density', groupColour = T, groupFill = T)
  
  ggsave(filename = "Cell positions in cartridge (contour).png",
         dpi = 300, 
         width = 24, 
         height = 16,
         units = "cm",
         plot = p_pos_out, 
         path = "./Image Output")
  
  p_pos2 <- ggplot(pos, aes(x=X1, y=X2, colour=Type)) +
            geom_point(shape = 20) + 
            # geom_density_2d() + 
            geom_rug(position="identity", linewidth = 0.1) + 
            theme_classic() +
            xlab("Well Position (X)") + 
            ylab("Well Position (Y)") + 
            theme(axis.text = element_text(size = 12)) + 
            theme(axis.title = element_text(size = 12)) +
            theme(plot.title = element_text(size = 14)) +
            theme(legend.position = "bottom") + 
            ggtitle(paste("Cell positions in cartridge -", 
                    folder_name))
  
  p_pos2_out <- ggMarginal(p_pos2, type='density', groupColour = T, groupFill = T)
  
  ggsave(filename = "Cell positions in cartridge.png",
         dpi = 300,
         width = 24, 
         height = 16,
         units = "cm",
         plot = p_pos2_out, 
         path = "./Image Output")
  
  
  ## d) filtering cells
  call_df <- img_df$poor_single_cell_call$cell_call_df %>% 
             mutate(noise_score = max_noise + mean_noise + ratio_noise,
                    Type = ifelse(noise_score == 0, "cell", "noise"))
  
  p_call_1 <- ggplot(call_df, 
                   aes(mean, max, colour = cluster, group = Type)) +
             geom_point(aes(shape = Type), alpha = 0.5) +
             scale_shape_manual(values = c(3,16)) +
             theme_classic() +
             xlab("Mean intensity") + 
             ylab("Max intensity") + 
             theme(axis.text = element_text(size = 12)) + 
             theme(axis.title = element_text(size = 12)) +
             theme(plot.title = element_text(size = 14)) +
             theme(legend.position = "bottom")
  
  p_call_2 <- ggplot(call_df, 
                     aes(ratio, max, colour = cluster, group = Type)) +
              geom_point(aes(shape = Type), alpha = 0.5) +
              scale_shape_manual(values = c(3,16)) +
              theme_classic() +
              xlab("Fluorescence region (%)") + 
              ylab("Max intensity") + 
              theme(axis.text = element_text(size = 12)) + 
              theme(axis.title = element_text(size = 12)) +
              theme(plot.title = element_text(size = 14)) +
              theme(legend.position = "bottom")
  
  p_call_3 <- ggplot(call_df, 
                     aes(ratio, mean, colour = cluster, group = Type)) +
              geom_point(aes(shape = Type), alpha = 0.5) +
              scale_shape_manual(values = c(3,16)) +
              theme_classic() +
              xlab("Fluorescence region (%)") + 
              ylab("Mean intensity") + 
              theme(axis.text = element_text(size = 12)) + 
              theme(axis.title = element_text(size = 12)) +
              theme(plot.title = element_text(size = 14)) +
              theme(legend.position = "bottom")
  
  p_call <- (p_call_1 + p_call_2 + p_call_3) + 
            plot_layout(guides = 'collect') + 
            plot_annotation(title = "Identify fluorescence noise") & 
            theme(legend.position = "bottom")
  
  ggsave(filename = "Cell and noise.png",
         dpi = 300,
         width = 24, 
         height = 10,
         units = "cm",
         plot = p_call, 
         path = "./Image Output")
  
  ## Sankey plot
  # / singlets
  temp_df <- data.frame(Type = ifelse(img_df$single_cell_info$FalseCall, "noise", "singlet"), 
                        Cluster = img_df$poor_single_cell_call$cell_call_df$cluster) %>% 
                        ggsankey::make_long(Type, Cluster)
  
  p_sankey_single <- ggplot(temp_df, aes(x = x, next_x = next_x, 
                                  node=node, next_node=next_node, 
                                  fill=factor(node), label = node)) +   
    geom_sankey(flow.alpha = .6, node.color = "gray30") +
    geom_sankey_label(size = 3, color = "white", fill = "gray40") +
    scale_fill_viridis_d() +
    theme_sankey(base_size = 18) +
    labs(x = NULL) +
    theme(legend.position = "none", 
          plot.title = element_text(hjust = .5)) +
    ggtitle("Singlets")
  
  ggsave(filename = "Sankey Singlet.png",
         dpi = 300,
         width = 15, 
         height = 10,
         units = "cm",
         plot = p_sankey_single, 
         path = "./Image Output")
  
  # / multiplets
  temp_df <- data.frame(Type = ifelse(img_df$multi_cell_info$FalseCall, "noise", "multiplet"), 
                        Cluster = img_df$poor_multi_cell_call$cell_call_df$cluster) %>% 
                        ggsankey::make_long(Type, Cluster)
  
  p_sankey_multi <- ggplot(temp_df, aes(x = x, next_x = next_x, 
                                         node=node, next_node=next_node, 
                                         fill=factor(node), label = node)) +   
    geom_sankey(flow.alpha = .6, node.color = "gray30") +
    geom_sankey_label(size = 3, color = "white", fill = "gray40") +
    scale_fill_viridis_d() +
    theme_sankey(base_size = 18) +
    labs(x = NULL) +
    theme(legend.position = "none", 
          plot.title = element_text(hjust = .5)) +
    ggtitle("Multiplets")
  
  ggsave(filename = "Sankey Multiplet.png",
         dpi = 300,
         width = 15, 
         height = 10,
         units = "cm",
         plot = p_sankey_multi, 
         path = "./Image Output")
  
  
  # 5. Finish
  return(list(processed_images_dataframe = img_df,
              grouped_BF_images_singlets = BF_split_single,
              grouped_FG_images_singlets = FG_split_single
              )
        )
 }

# /// ------------ Run the whole pipeline ------------ 
results <- func_workflow(input_info)

# /// ------------ Data frames to save ------------ 
results_df <- results$processed_images_dataframe[sapply(results$processed_images_dataframe, class) == "data.frame"]
results_df[["input_info"]] <- input_info
save(input_info, results_df, file = paste0("./Image Output/Out_", folder_name, ".RData"))


 