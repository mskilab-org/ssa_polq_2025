library(Flow)
library(dplyr)
library(tidyr)
library(skitools)
library(RColorBrewer)
library(ggnewscale)
library(skitools)
library(ggpubr)
library(pbmcapply)
library(data.table)
library(scales)
library(ComplexHeatmap)
library(circlize)
library(pals)
library(colorspace)
#library(scico)
#library(viridis)

add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha))  
}

cividis_hcl <- function(n) {
  i <- seq(1, 0, length.out = n)
  hex(polarLUV(
    L = 92 - (92 - 13) * i,
    C = approx(c(1, 0.9, 0.5, 0), c(30, 50, 0, 95), xout = i)$y,
    H = c(255, 75)[1 + (i < 0.5)]
  ), fix = TRUE)
}


d <- readRDS("/gpfs/commons/home/sbrylka/Projects/homeology/final_data/prediction.rds")

d <- readRDS("~/Projects/homeology/stash/entire_amplicon.rds")


#homology complexheatmap
#all regions
data_matrix <- merge(expand.grid(
  position = -250:0,
  deletion_size = 1:250
),d, by = c("position", "deletion_size"), all.x = TRUE)

## data_matrix <- merge(expand.grid(
##   position = 1:503,
##   deletion_size = 1:503
## ),d, by = c("position", "deletion_size"), all.x = TRUE)


data_matrix$position <- as.numeric(data_matrix$position)
data_matrix$deletion_size<- as.numeric(data_matrix$deletion_size)

total_data_matrix = data_matrix %>%
  rowwise() %>%
  mutate(homology = as.numeric(case_when(
           !is.na(homology) ~ homology,
           abs(position) - deletion_size > 0 ~ NA,
           T ~ 0
         ))) %>%
  #dplyr::select(-c("homeology")) %>%
  as.data.table()

small_data_matrix = dt1 %>%
  mutate(position = position - 251) %>%
  filter(position)


  filter(position %in% c(-69:0),
         deletion_size %in% c(1:100))



#' small data table
small_data_matrix = total_data_matrix %>%
v  filter(position %in% c(-69:0),
         deletion_size %in% c(1:100))
medium_data_matrix = total_data_matrix %>%
  filter(position %in% c(-200:0),
         deletion_size %in% c(1:200))
#' large data table
large_data_matrix = total_data_matrix %>%
  filter(position %in% c(-250:0),
         deletion_size %in% c(1:250))

#' small homology
small.homology = small_data_matrix %>%
  dplyr::select(-c("homeology")) %>%
  pivot_wider(
    names_from = "deletion_size",
    values_from = "homology",
    names_sort = T
  ) %>%
  as.data.table %>%
  setkey(position) %>%
  as.matrix(rownames = T)
#' small homeology
small.homeology = small_data_matrix %>%
  dplyr::select(-c("homology")) %>%
  pivot_wider(
    names_from = "deletion_size",
    values_from = "homeology",
    names_sort = T
  ) %>%
  as.data.table %>%
  setkey(position) %>%
  as.matrix(rownames = T)

#' medium homology
medium.homology = medium_data_matrix %>%
  dplyr::select(-c("homeology")) %>%
  pivot_wider(
    names_from = "deletion_size",
    values_from = "homology",
    names_sort = T
  ) %>%
  as.data.table %>%
  setkey(position) %>%
  as.matrix(rownames = T)
#' large homeology
medium.homeology = medium_data_matrix %>%
  dplyr::select(-c("homology")) %>%
  pivot_wider(
    names_from = "deletion_size",
    values_from = "homeology",
    names_sort = T
  ) %>%
  as.data.table %>%
  setkey(position) %>%
  as.matrix(rownames = T)

#' large homology
large.homology = large_data_matrix %>%
  dplyr::select(-c("homeology")) %>%
  pivot_wider(
    names_from = "deletion_size",
    values_from = "homology",
    names_sort = T
  ) %>%
  as.data.table %>%
  setkey(position) %>%
  as.matrix(rownames = T)
#' large homeology
large.homeology = large_data_matrix %>%
  dplyr::select(-c("homology")) %>%
  pivot_wider(
    names_from = "deletion_size",
    values_from = "homeology",
    names_sort = T
  ) %>%
  as.data.table %>%
  setkey(position) %>%
  as.matrix(rownames = T)

#' jrafailov Tuesday, Jul 30, 2024 10:37:59 AM
#' lets make the small matrix first
  

small.homology.breaks <- c(0:max(small.homology, na.rm = T))
small.homeology.breaks <- c(0:max(small.homeology, na.rm = T))
n.breaks <- length(small.homology.breaks) - 2
breaks.left = length(small.homeology.breaks) - 2 - n.breaks
total.breaks = n.breaks + breaks.left + 2

color0 = "#e3e3e3"

color0 = "#f9f9f9"

col = rev(sequential_hcl(n.breaks + 3, palette = "Mako")[-c(1,2, n.breaks + 3)])
col2 = colorRamp2(
  breaks = c(max(small.homeology.breaks) - breaks.left,
             max(small.homeology.breaks)),
  colors = c(col[length(col)], "palevioletred"))
col2 = col2(c((max(small.homeology.breaks) - breaks.left + 1):
                max(small.homeology.breaks)))
col_fun_prop=colorRamp2(small.homeology.breaks,c(color0,color0, col, col2), space = "RGB")

saveRDS(col_fun_prop, "/gpfs/commons/home/sbrylka/Projects/homeology/jr_stash/col_fun_prop/small.rds")

col_labels <- function(matrix){
  colnames = as.numeric(colnames(matrix))
  pos_of_min = which(abs(colnames) == min(abs(colnames)))
  pos_of_5_mod = which(colnames %% 5 == 0)
  keep = c(pos_of_min, pos_of_5_mod)
  colnames[-keep] <- ""
  return(colnames)
}
row_labels <- function(matrix){
  rownames = as.numeric(rownames(matrix))
  pos_of_max = which(abs(rownames) == max(abs(rownames)))
  pos_of_5_mod = which(rownames %% 5 == 0)
  keep = c(pos_of_max, pos_of_5_mod)
  rownames[-keep] <- ""
  return(rownames)
}
                                        #[-c(1,2,3,4)])

                                        #col = rev(cividis_hcl(9))

#sequential_hcl(n.breaks - 2, h = c(260, 60), c = 60, l = c(40, 95), power = 1)

#[-c(10,11,12,13,14,15)]

#[-c(8, 9, 10, 11, 12)]

#col_fun_prop=colorRamp2(small.homeology.breaks,c(color0,color0, col, col2), space = "RGB")

pdf(file=paste0("~/public_html/plots/homology/indel_heatmaps/", "HOMOLOGY_small_amplicon_prediction.pdf"), width = unit(6, "in"), height = unit(6, "in"))
lgd = Legend(
  at = small.homology.breaks,
  #col_fun = col_fun_prop,
  legend_gp = gpar(fill = col_fun_prop(small.homology.breaks)),
  #breaks = c(0:11),
  #at = c(0,  1e-6 - 1e-7,5e-6,1e-5,5e-5,1e-4,5e-4,1e-3,5e-3,1e-2,5e-2,0.1,0.5),
  #labels = c("0",  "1e-6","5e-6","1e-5","5e-5","1e-4","5e-4","1e-3","5e-3","1e-2","5e-2","0.1","0.5"),
  title = paste0("Homology"),
  title_position = "lefttop-rot",
  legend_height = unit(0.5, "in"))
Heatmap(
  small.homology,
  row_order = rev(rownames(small.homology)),
  row_title = "Distance from cut site",
  row_title_side = "right",
  col = col_fun_prop,
  #top_annotation = deletion_histogram,
  #right_annotation = row_histogram,
  column_title_side = "top",
  column_title = "Deletion size",
  rect_gp = gpar(col= "white", lwd = 0.01),
  column_order = colnames(small.homology),
  show_column_names = T,
  show_row_names = T,
  column_labels = col_labels(small.homology),
  row_labels = row_labels(small.homology),
  column_names_side = "top",
  na_col = "white",
  width = ncol(small.homology)*unit(0.8, "mm"),
  height = nrow(small.homology)*unit(0.8, "mm"),
  show_heatmap_legend = F,
  row_names_gp = gpar(fontsize = 8),
  column_names_gp = gpar(fontsize = 8)
) -> hm
draw(hm, heatmap_legend_list = lgd, heatmap_legend_side = "left")
dev.off()

pdf(file=paste0("~/public_html/plots/homology/indel_heatmaps/", "HOMEOLOGY_small_amplicon_prediction.pdf"), width = unit(6, "in"), height = unit(6, "in"))
lgd = Legend(
  at = small.homeology.breaks,
  #col_fun = col_fun_prop,
  legend_gp = gpar(fill = col_fun_prop(small.homeology.breaks)),
  #breaks = c(0:11),
  #at = c(0,  1e-6 - 1e-7,5e-6,1e-5,5e-5,1e-4,5e-4,1e-3,5e-3,1e-2,5e-2,0.1,0.5),
  #labels = c("0",  "1e-6","5e-6","1e-5","5e-5","1e-4","5e-4","1e-3","5e-3","1e-2","5e-2","0.1","0.5"),
  title = paste0("Homeology"),
  title_position = "lefttop-rot",
  legend_height = unit(0.5, "in"))
Heatmap(
  small.homeology,
  row_order = rev(rownames(small.homeology)),
  row_title = "Distance from cut site",
  row_title_side = "right",
  col = col_fun_prop,
  rect_gp = gpar(col = "white", lwd = 0.01),
  #top_annotation = deletion_histogram,
  #right_annotation = row_histogram,
  column_title_side = "top",
  column_title = "Deletion size",
  column_order = colnames(small.homeology),
  show_column_names = T,
  show_row_names = T,
  column_labels = col_labels(small.homeology),
  row_labels = row_labels(small.homeology),
  column_names_side = "top",
  na_col = "white",
  width = ncol(small.homeology)*unit(0.8, "mm"),
  height = nrow(small.homeology)*unit(0.8, "mm"),
  show_heatmap_legend = F,
  row_names_gp = gpar(fontsize = 8),
  column_names_gp = gpar(fontsize = 8)
) -> hm
draw(hm, heatmap_legend_list = lgd, heatmap_legend_side = "left")
dev.off()

medium.homology.breaks <- c(0:max(medium.homology, na.rm = T))
medium.homeology.breaks <- c(0:max(medium.homeology, na.rm = T))
n.breaks <- length(medium.homology.breaks) - 2
breaks.left = length(medium.homeology.breaks) - 2 - n.breaks
total.breaks = n.breaks + breaks.left + 2

color0 = "#e3e3e3"

color0 = "#f2f2f2"

col = rev(sequential_hcl(n.breaks + 3, palette = "Mako")[-c(1,2, n.breaks + 3)])
col2 = colorRamp2(
  breaks = c(max(medium.homeology.breaks) - breaks.left,
             max(medium.homeology.breaks)),
  colors = c(col[length(col)], "palevioletred"))
col2 = col2(c((max(medium.homeology.breaks) - breaks.left + 1):
                max(medium.homeology.breaks)))
col_fun_prop=colorRamp2(medium.homeology.breaks,c(color0,color0, col, col2), space = "RGB")
saveRDS(col_fun_prop, "/gpfs/commons/home/sbrylka/Projects/homeology/jr_stash/col_fun_prop/medium.rds")


col_labels <- function(matrix){
  colnames = as.numeric(colnames(matrix))
  pos_of_min = which(abs(colnames) == min(abs(colnames)))
  pos_of_mod = which(colnames %% 10 == 0)
  keep = c(pos_of_min, pos_of_mod)
  colnames[-keep] <- ""
  return(colnames)
}
row_labels <- function(matrix){
  rownames = as.numeric(rownames(matrix))
  pos_of_max = which(abs(rownames) == max(abs(rownames)))
  pos_of_mod = which(rownames %% 10 == 0)
  keep = c(pos_of_max, pos_of_mod)
  rownames[-keep] <- ""
  return(rownames)
}

                                        #[-c(1,2,3,4)])

                                        #col = rev(cividis_hcl(9))

#sequential_hcl(n.breaks - 2, h = c(260, 60), c = 60, l = c(40, 95), power = 1)

#[-c(10,11,12,13,14,15)]

#[-c(8, 9, 10, 11, 12)]

#col_fun_prop=colorRamp2(small.homeology.breaks,c(color0,color0, col, col2), space = "RGB")

pdf(file=paste0("~/public_html/plots/homology/indel_heatmaps/", "HOMOLOGY_medium_amplicon_prediction.pdf"), width = unit(6, "in"), height = unit(6, "in"))
lgd = Legend(
  at = medium.homology.breaks,
  #col_fun = col_fun_prop,
  legend_gp = gpar(fill = col_fun_prop(medium.homology.breaks)),
  #breaks = c(0:11),
  #at = c(0,  1e-6 - 1e-7,5e-6,1e-5,5e-5,1e-4,5e-4,1e-3,5e-3,1e-2,5e-2,0.1,0.5),
  #labels = c("0",  "1e-6","5e-6","1e-5","5e-5","1e-4","5e-4","1e-3","5e-3","1e-2","5e-2","0.1","0.5"),
  title = paste0("Homology"),
  title_position = "lefttop-rot",
  legend_height = unit(0.5, "in"))
Heatmap(
  medium.homology,
  row_order = rev(rownames(medium.homology)),
  row_title = "Distance from cut site",
  row_title_side = "right",
  col = col_fun_prop,
  #top_annotation = deletion_histogram,
  #right_annotation = row_histogram,
  column_title_side = "top",
  column_title = "Deletion size",
  #rect_gp = gpar(col= "white", lwd = 0.01),
  column_order = colnames(medium.homology),
  show_column_names = T,
  show_row_names = T,
  column_labels = col_labels(medium.homology),
  row_labels = row_labels(medium.homology),
  column_names_side = "top",
  na_col = "white",
  width = ncol(medium.homology)*unit(0.43, "mm"),
  height = nrow(medium.homology)*unit(0.43, "mm"),
  show_heatmap_legend = F,
  row_names_gp = gpar(fontsize = 8),
  column_names_gp = gpar(fontsize = 8)
) -> hm
draw(hm, heatmap_legend_list = lgd, heatmap_legend_side = "left")
dev.off()

pdf(file=paste0("~/public_html/plots/homology/indel_heatmaps/", "HOMEOLOGY_medium_amplicon_prediction.pdf"), width = unit(6, "in"), height = unit(6, "in"))
lgd = Legend(
  at = medium.homeology.breaks,
  #col_fun = col_fun_prop,
  legend_gp = gpar(fill = col_fun_prop(medium.homeology.breaks)),
  #breaks = c(0:11),
  #at = c(0,  1e-6 - 1e-7,5e-6,1e-5,5e-5,1e-4,5e-4,1e-3,5e-3,1e-2,5e-2,0.1,0.5),
  #labels = c("0",  "1e-6","5e-6","1e-5","5e-5","1e-4","5e-4","1e-3","5e-3","1e-2","5e-2","0.1","0.5"),
  title = paste0("Homeology"),
  title_position = "lefttop-rot",
  legend_height = unit(0.5, "in"))
Heatmap(
  medium.homeology,
  row_order = rev(rownames(medium.homeology)),
  row_title = "Distance from cut site",
  row_title_side = "right",
  col = col_fun_prop,
  #rect_gp = gpar(col = "white", lwd = 0.01),
  #top_annotation = deletion_histogram,
  #right_annotation = row_histogram,
  column_title_side = "top",
  column_title = "Deletion size",
  column_order = colnames(medium.homeology),
  show_column_names = T,
  show_row_names = T,
  column_labels = col_labels(medium.homeology),
  row_labels = row_labels(medium.homeology),
  column_names_side = "top",
  na_col = "white",
  width = ncol(medium.homeology)*unit(0.43, "mm"),
  height = nrow(medium.homeology)*unit(0.43, "mm"),
  show_heatmap_legend = F,
  row_names_gp = gpar(fontsize = 8),
  column_names_gp = gpar(fontsize = 8)
) -> hm
draw(hm, heatmap_legend_list = lgd, heatmap_legend_side = "left")
dev.off()




#' jrafailov Tuesday, Jul 30, 2024 10:37:59 AM
#' lets make the small matrix first
  

large.homology.breaks <- c(0:max(large.homology, na.rm = T))
large.homeology.breaks <- c(0:max(large.homeology, na.rm = T))
n.breaks = length(large.homology.breaks) - 2
breaks.left = length(large.homeology.breaks) - 2 - n.breaks
total.breaks = n.breaks + breaks.left + 2


color0 = "#f2f2f2"

color0 = "#e3e3e3"
col = rev(sequential_hcl(n.breaks + 3, palette = "Mako")[-c(1,2, n.breaks + 3)])
col2 = colorRamp2(
  breaks = c(max(large.homeology.breaks) - breaks.left,
             max(large.homeology.breaks)),
  colors = c(col[length(col)], "palevioletred"))
col2 = col2(c((max(large.homeology.breaks) - breaks.left + 1):
              max(large.homeology.breaks)))
col_fun_prop=colorRamp2(large.homeology.breaks,c(color0,color0, col, col2), space = "RGB")

saveRDS(col_fun_prop, "/gpfs/commons/home/sbrylka/Projects/homeology/jr_stash/col_fun_prop/large.rds")

col_labels <- function(matrix){
  colnames = as.numeric(colnames(matrix))
  pos_of_min = which(abs(colnames) == min(abs(colnames)))
  pos_of_mod = which(colnames %% 10 == 0)
  keep = c(pos_of_min, pos_of_mod)
  colnames[-keep] <- ""
  return(colnames)
}
row_labels <- function(matrix){
  rownames = as.numeric(rownames(matrix))
  pos_of_max = which(abs(rownames) == max(abs(rownames)))
  pos_of_mod = which(rownames %% 10 == 0)
  keep = c(pos_of_max, pos_of_mod)
  rownames[-keep] <- ""
  return(rownames)
}

pdf(file=paste0("~/public_html/plots/homology/indel_heatmaps/", "HOMOLOGY_large_amplicon_prediction.pdf"), width = unit(6, "in"), height = unit(6, "in"))
lgd = Legend(
  at = large.homology.breaks,
  #col_fun = col_fun_prop,
  legend_gp = gpar(fill = col_fun_prop(large.homology.breaks)),
  #breaks = c(0:11),
  #at = c(0,  1e-6 - 1e-7,5e-6,1e-5,5e-5,1e-4,5e-4,1e-3,5e-3,1e-2,5e-2,0.1,0.5),
  #labels = c("0",  "1e-6","5e-6","1e-5","5e-5","1e-4","5e-4","1e-3","5e-3","1e-2","5e-2","0.1","0.5"),
  title = paste0("Homology"),
  title_position = "lefttop-rot",
  legend_height = unit(0.5, "in"))
Heatmap(
  large.homology,
  #rect_gp = gpar(col = "white", lwd = 0.1),
  row_order = rev(rownames(large.homology)),
  row_title = "Distance from cut site",
  row_title_side = "right",
  col = col_fun_prop,
  #top_annotation = deletion_histogram,
  #right_annotation = row_histogram,
  column_title_side = "top",
  column_title = "Deletion size",
  column_order = colnames(large.homology),
  show_column_names = T,
  show_row_names = T,
  column_labels = col_labels(large.homology),
  row_labels = row_labels(large.homology),
  column_names_side = "top",
  na_col = "white",
  width = ncol(large.homology)*unit(0.35, "mm"),
  height = nrow(large.homology)*unit(0.35, "mm"),
  show_heatmap_legend = F,
  row_names_gp = gpar(fontsize = 6),
  column_names_gp = gpar(fontsize = 6)
) -> hm
draw(hm, heatmap_legend_list = lgd, heatmap_legend_side = "left")
dev.off()

pdf(file=paste0("~/public_html/plots/homology/indel_heatmaps/", "HOMEOLOGY_large_amplicon_prediction.pdf"), width = unit(6, "in"), height = unit(6, "in"))
lgd = Legend(
  at = large.homeology.breaks,
  #col_fun = col_fun_prop,
  legend_gp = gpar(fill = col_fun_prop(large.homeology.breaks)),
  #breaks = c(0:11),
  #at = c(0,  1e-6 - 1e-7,5e-6,1e-5,5e-5,1e-4,5e-4,1e-3,5e-3,1e-2,5e-2,0.1,0.5),
  #labels = c("0",  "1e-6","5e-6","1e-5","5e-5","1e-4","5e-4","1e-3","5e-3","1e-2","5e-2","0.1","0.5"),
  title = paste0("bp of hom(e)ology"),
  title_position = "lefttop-rot",
  legend_height = unit(0.5, "in"))
Heatmap(
  large.homeology,
  #rect_gp = gpar(col = "white", lwd = 0.1),
  row_order = rev(rownames(large.homeology)),
  row_title = "Distance from cut site",
  row_title_side = "right",
  col = col_fun_prop,
  #top_annotation = deletion_histogram,
  #right_annotation = row_histogram,
  column_title_side = "top",
  column_title = "Deletion size",
  column_order = colnames(large.homeology),
  show_column_names = T,
  show_row_names = T,
  column_labels = col_labels(large.homeology),
  row_labels = row_labels(large.homeology),
  column_names_side = "top",
  na_col = "white",
  width = ncol(large.homeology)*unit(0.35, "mm"),
  height = nrow(large.homeology)*unit(0.35, "mm"),
  show_heatmap_legend = F,
  row_names_gp = gpar(fontsize = 6),
  column_names_gp = gpar(fontsize = 6)
) -> hm
draw(hm, heatmap_legend_list = lgd, heatmap_legend_side = "left")
dev.off()




#' jrafailov Monday, Aug 12, 2024 03:15:32 PM
#' custom col_fun_prop for the graph's manisha wanted

#' small

small.homology.breaks <- c(0:max(small.homology, na.rm = T))
small.homeology.breaks <- c(0:max(small.homeology, na.rm = T))
n.breaks = length(small.homology.breaks) - 2
breaks.left = length(small.homeology.breaks) - 2 - n.breaks
total.breaks = n.breaks + breaks.left + 2

color0 = "#e3e3e3"
col = rev(sequential_hcl(n.breaks + 3, palette = "Mako")[-c(1,2, n.breaks + 3)])
color1 = col[1]
color2 = col[length(col)]

col_fun_prop=colorRamp2(small.homeology.breaks,
                        c(rep(color0, 2),
                          rep(color1, 3),
                          rep(color2, 6)), space = "RGB")


#' large

large.homology.breaks <- c(0:max(large.homology, na.rm = T))
large.homeology.breaks <- c(0:max(large.homeology, na.rm = T))
n.breaks = length(large.homology.breaks) - 2
breaks.left = length(large.homeology.breaks) - 2 - n.breaks
total.breaks = n.breaks + breaks.left + 2


color0 = "#e3e3e3"
col = rev(sequential_hcl(n.breaks + 3, palette = "Mako")[-c(1,2, n.breaks + 3)])
color1 = col[1]
color2 = col[length(col)]

## col = rev(sequential_hcl(n.breaks + 3, palette = "Mako")[-c(1,2, n.breaks + 3)])
## col2 = colorRamp2(
##   breaks = c(max(large.homeology.breaks) - breaks.left,
##              max(large.homeology.breaks)),
##   colors = c(col[length(col)], "palevioletred"))
## col2 = col2(c((max(large.homeology.breaks) - breaks.left + 1):
##               max(large.homeology.breaks)))

col_fun_prop=colorRamp2(large.homeology.breaks,c(rep(color0, 2), rep(color1, 3), rep(color2, 7)), space = "RGB")
col_fun_prop=colorRamp2(large.homeology.breaks,c(rep(color0, 5), rep(color2, 7)), space = "RGB")


saveRDS(col_fun_prop, "/gpfs/commons/home/sbrylka/Projects/homeology/jr_stash/col_fun_prop/custom.rds")


heatmap.fn <- function(table,
                       min, max,
                       order = NULL,
                       min.homo = 2,
                       stacked = F,
                       homo_homeo = "homeo"){
  n_geno = unique(table$genotype)

  type.n <-case_when(
    homo_homeo == "homeo" ~ "homeology",
    T ~ "homology")

  table.filtered <- table %>%
    filter(type == type.n,
           deletion_size >= min,
           deletion_size <= max,
           bp > min.homo) %>%
    select(genotype, bp, counts, deletion_density, replicate) %>%
    group_by(genotype,replicate, bp) %>%
    #filter(any(counts > 0)) %>%
    summarize(counts = sum(counts),
              deletion_density = sum(deletion_density)) %>%
    ungroup %>%
    mutate(score = deletion_density * bp) %>%
    as.data.table()

  scores.per.replicate <- table.filtered %>%
    group_by(genotype, replicate) %>%
    summarize(score = sum(score),
              deletion_density = sum(deletion_density)) %>%
    mutate(sd = sd(score),
           mean = mean(score),
           se = sd / sqrt(n()))

  table %>%
    filter(type == type.n,
           deletion_size >= min,
           deletion_size <= max,
           bp > min.homo) %>%
    select(genotype, bp, counts, total_dels, replicate, deletion_density_gt) %>%
    group_by(genotype, bp) %>%
    summarize(deletion_density_gt = sum(deletion_density_gt)) %>%
    as.data.table() -> bp.to.plot.all

  table %>%
    filter(type == type.n,
           deletion_size >= min,
           deletion_size <= max,
           bp > min.homo) %>%
    select(genotype, bp, counts, total_dels, replicate, deletion_density_gt, deletion_density) %>%
    group_by(genotype, bp, replicate) %>%
    summarize(deletion_density_gt = sum(deletion_density_gt),
              deletion_density = sum(deletion_density)) %>%
    as.data.table() -> bp.to.plot.rep

  
  if(stacked == T){

    stats = scores.per.replicate %>%
      group_by(genotype) %>%
      t_test(genotype ~ score)
    
  } else if (stacked == F){
    bar <- geom_bar(aes(fill = factor(bp), color = factor(bp)), position = "dodge", stat = "identity")
  }

  ggplot(bp.to.plot, aes(x = genotype, y = deletion_density_gt)) +
    geom_bar(aes(fill = factor(bp), color = factor(bp)), stat = "identity") +
    #stat_pvalue_manual(stats, label = "p.adj.signif") +
    scale_fill_manual(values = col_fun_prop(unique(bp.to.plot$bp))) +
    scale_color_manual(values = rep("black", length(unique(table.filtered$bp)))) +
    theme_bw() +
    theme(panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.major.y = element_line(size = .1, color = "black"),
          panel.grid.minor.y = element_blank(),
          strip.background = element_blank(),
          strip.text = element_text(size = 8)
          #plot.margin = unit(c(0, 0, 0, 0), "cm")
          ) +
    labs(title = paste0("Deletion sizes ", min, "-", max, ", inclusive")) +
    ylab("% of deletions") -> all

  ggplot(bp.to.plot.rep, aes(x = replicate, y = deletion_density)) +
    geom_bar(aes(fill = factor(bp), color = factor(bp)), stat = "identity") +
                                        #stat_pvalue_manual(stats, label = "p.adj.signif") +
    scale_fill_manual(values = col_fun_prop(unique(bp.to.plot$bp))) +
    scale_color_manual(values = rep("black", length(unique(table.filtered$bp)))) +
    facet_wrap( ~ genotype, scales = "free_x", ncol = length(unique(bp.to.plot.rep$genotype))) +
    theme_bw() +
    theme(panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.major.y = element_line(size = .1, color = "black"),
          panel.grid.minor.y = element_blank(),
          strip.background = element_blank(),
          strip.text = element_text(size = 8)
                                        #plot.margin = unit(c(0, 0, 0, 0), "cm")
          ) +
    labs(title = paste0("Deletion sizes ", min, "-", max, ", inclusive")) +
    ylab("% of deletions") -> all.rep

  ggplot(scores.per.replicate, aes(x = genotype, y = score)) +
    geom_bar(stat = "summary", color = "black", fill = "skyblue", width = 0.7) +
    geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.4) +
    stat_pvalue_manual(stats, label = "p.adj.signif") +
    theme_bw() +
    theme(panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.major.y = element_line(size = .1, color = "black"),
          panel.grid.minor.y = element_blank(),
          strip.background = element_blank(),
          strip.text = element_text(size = 8)
                                        #plot.margin = unit(c(0, 0, 0, 0), "cm")
          ) +
    labs(title = paste0("Deletion sizes ", min, "-", max, ", inclusive")) +
    ylab("Score") -> all.2

  skitools::ppdf(print(all))
  skitools::ppdf(print(all.2))
  skitools::ppdf(print(all.rep))
  
  
  return(all)
  
}



















