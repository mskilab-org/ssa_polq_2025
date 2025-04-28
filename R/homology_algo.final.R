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

#' need to generate the new data...
seq1 <- "TTACCTCTCTAGTCTGTGCTAGCTCttccagccccctgtcatggcatcttccaggggtccgagagctcagctagtcttcttcctccaacccgggcccctatgtccacttcaggacagcatgtttgctgcctccagggatcctgtgtccccgagctgggaccaccttatattcccagggccggttaatgtggctctggttctgggtacttttatctgtcccctccaccccacagtggggccactagggacaggattggtgacagaaaagccccatccttaggcctcctccttcctagtctcctgatattgggtctaacccccacctcctgttaggcagattccttatctggtgacacacccccatttcctggagccatctctctccttgccagaacctctaaggtttgcttacgatggagccagagaggatcctgggagggagagcttggcagggggtgggagggaagggggggatgcgtgacctGCCCGGTTCTCAGTGGCCA" %>% toupper()
nchar(seq1)

seq2 <- "ttaatgtggctctggttctgggtacttttatctgtcccctccaccccacagtGGGGCCACTAGGGACAGGATTggtgacagaaaagccccatccttaggcctcctccttcctagtctcctgatattgggtctaacccccacctcctgttaggcagattccttatctggtgacacacccccatttcctggagccatctctctccttgccagaacctctaaggtttgcttacgatggagccagagag" %>% toupper()
nchar(seq2)

## Functions
find_homology <- function(position, deletion, seq, debug = F) {
  if (position + deletion > 503) {
    return(NA)
  }
  count.x <- 0
  if (debug) print("5' retained sequence vs 3' deleted sequence")
  if (debug) {
    print(paste0("starting 5' retained sequence: ", substr(seq, position - count.x - 1, position - count.x - 1)))
    print(paste0("starting 3' deleted sequence: ", substr(seq, position + deletion - count.x - 1, position + deletion - count.x - 1)))
  }
  while (substr(seq, position - count.x - 1, position - count.x - 1) == substr(seq, position + deletion - count.x - 1, position + deletion - count.x - 1)) {
    if (debug) print(paste0(substr(seq, position - count.x - 1, position - count.x - 1), " ", substr(seq, position + deletion - count.x - 1, position + deletion - count.x - 1), " count:", count.x + 1))
    count.x <- count.x + 1
  }
  count.y <- 0 #' bidirectional search!
  if (debug) print("5' deleted sequence vs 3' retained sequence")
  if (debug) {
    print(paste0("starting 5' deleted sequence: ", substr(seq, position + deletion + count.y, position + deletion + count.y)))
    print(paste0("starting 3' retained sequence: ", substr(seq, position + count.y, position + count.y)))
  }
  while (substr(seq, position + deletion + count.y, position + deletion + count.y) == substr(seq, position + count.y, position + count.y)) {
    if (debug) print(paste0(substr(seq, position + deletion + count.y, position + deletion + count.y), " ", substr(seq, position + count.y, position + count.y), " count:", count.y + 1))
    count.y <- count.y + 1
  }
  return(max(count.x, count.y))
}

find_homeology <- function(position, deletion, seq, debug = F, thresh) {
  if (position + deletion > 503) {
    return(NA)
  }
  max.x <- 0
  count.x <- 1
  match.x <- 0
  x.s <- ""
  x.t <- ""
  if (debug) print("5' retained sequence vs 3' deleted sequence")
  while (count.x <= deletion) {
    s1 <- substr(seq, position - count.x, position - count.x)
    s2 <- substr(seq, position - count.x + deletion, position - count.x + deletion)
    # if(debug) print(paste0(s1, " ", s2))
    match.x <- match.x + mapply(function(x, y) sum(x == y), strsplit(s1, ""), strsplit(s2, ""))
    if (match.x / count.x >= thresh & s1 == s2) {
      max.x <- count.x
    }
    count.x <- count.x + 1
    x.s <- paste0(s1, x.s)
    x.t <- paste0(s2, x.t)
  }
  if (debug) {
    print(paste0(
      stri_reverse(substring(stri_reverse(x.s), 1, 12)),
      "  ",
      stri_reverse(substring(stri_reverse(x.t), 1, 12))
    ))
  }
  max.y <- 0
  count.y <- 1
  match.y <- 0
  y.s <- ""
  y.t <- ""
  if (debug) print("5' deleted sequence vs 3' retained sequence")
  while (count.y <= deletion) {
    s1 <- substr(seq, position + count.y - 1, position + count.y - 1)
    s2 <- substr(seq, position + count.y + deletion - 1, position + count.y + deletion - 1)
    # if(debug) print(paste0(s1, " ", s2))
    match.y <- match.y + mapply(function(x, y) sum(x == y), strsplit(s1, ""), strsplit(s2, ""))
    if (match.y / count.y >= thresh & s1 == s2) {
      max.y <- count.y
    }
    count.y <- count.y + 1
    y.s <- paste0(y.s, s1)
    y.t <- paste0(y.t, s2)
  }
  if (debug) {
    print(paste0(
      substring(y.s, 1, 12),
      "  ",
      substring(y.t, 1, 12)
    ))
  }
  if (debug) print(paste0(max.x, " vs. ", max.y))
  return(max(max.x, max.y))
}

hbond <- function(nucleotide) {
  nucleotide <- nucleotide %>% toupper()
  if (nucleotide %in% c("A", "T")) {
    return(2)
  } else if (nucleotide %in% c("C", "G")) {
    return(3)
  } else {
    return(0)
  }
}

find_homology.hbonds <- function(position, deletion, seq, debug = F) {
  if (position + deletion > 503) {
    return(NA)
  }
  count.x <- 0
  if (debug) print("5' retained sequence vs 3' deleted sequence")
  if (debug) {
    print(paste0("starting 5' retained sequence: ", substr(seq, position - count.x - 1, position - count.x - 1)))
    print(paste0("starting 3' deleted sequence: ", substr(seq, position + deletion - count.x - 1, position + deletion - count.x - 1)))
  }
  left.pos.x <- substr(seq, position - count.x - 1, position - count.x - 1)
  right.pos.x <- substr(seq, position + deletion - count.x - 1, position + deletion - count.x - 1)
  hbonds.x <- 0
  while (left.pos.x == right.pos.x) {
    if (debug) print(paste0(left.pos.x, " ", right.pos.x, " count: ", hbonds.x + hbond(left.pos.x)))
    count.x <- count.x + 1
    left.pos.x <- substr(seq, position - count.x - 1, position - count.x - 1)
    right.pos.x <- substr(seq, position + deletion - count.x - 1, position + deletion - count.x - 1)
    hbonds.x <- hbonds.x + hbond(left.pos.x)
  }
  if (debug) print("5' deleted sequence vs 3' retained sequence")
  count.y <- 0 #' bidirectional search!
  if (debug) {
    print(paste0("starting 5' deleted sequence: ", substr(seq, position + deletion + count.y, position + deletion + count.y)))
    print(paste0("starting 3' retained sequence: ", substr(seq, position + count.y, position + count.y)))
  }
  left.pos.y <- substr(seq, position + deletion + count.y, position + deletion + count.y)
  right.pos.y <- substr(seq, position + count.y, position + count.y)
  hbonds.y <- 0
  while (left.pos.y == right.pos.y) {
    if (debug) print(paste0(left.pos.y, " ", right.pos.y, " count: ", hbonds.y + hbond(left.pos.y)))
    count.y <- count.y + 1
    hbonds.y <- hbonds.y + hbond(left.pos.y)
    left.pos.y <- substr(seq, position + deletion + count.y, position + deletion + count.y)
    right.pos.y <- substr(seq, position + count.y, position + count.y)
  }
  return(max(hbonds.x, hbonds.y))
}

find_homeology.hbonds <- function(position, deletion, seq, debug = F, thresh = 0.80) {
  if (position + deletion > 503) {
    return(NA)
  }
  max.x <- 0
  count.x <- 1
  match.x <- 0
  hbonds.x <- 0
  hbonds.evaluated.x <- 0
  x.s <- ""
  x.t <- ""
  if (debug) print("5' retained sequence vs 3' deleted sequence")
  while (count.x <= deletion) {
    s1 <- substr(seq, position - count.x, position - count.x)
    s2 <- substr(seq, position - count.x + deletion, position - count.x + deletion)
    if (debug) print(paste0(s1, " ", s2))
    match.x <- match.x + mapply(function(x, y) sum(x == y), strsplit(s1, ""), strsplit(s2, ""))
    if (match.x / count.x >= thresh & s1 == s2) {
      # hbonds.x <- hbonds.x + hbond(s1)
      # max.x <- count.x
      hbonds.x <- hbonds.x + hbond(s1)
      hbonds.evaluated.x <- hbonds.evaluated.x + hbond(s1)
    }
    count.x <- count.x + 1
    x.s <- paste0(s1, x.s)
    x.t <- paste0(s2, x.t)
  }
  if (debug) {
    print(paste0(
      stri_reverse(substring(stri_reverse(x.s), 1, 12)),
      "  ",
      stri_reverse(substring(stri_reverse(x.t), 1, 12))
    ))
  }
  max.y <- 0
  count.y <- 1
  match.y <- 0
  hbonds.y <- 0
  hbonds.evaluated.y <- 0
  y.s <- ""
  y.t <- ""
  if (debug) print("5' deleted sequence vs 3' retained sequence")
  while (count.y <= deletion) {
    s1 <- substr(seq, position + count.y - 1, position + count.y - 1)
    s2 <- substr(seq, position + count.y + deletion - 1, position + count.y + deletion - 1)
    if (debug) print(paste0(s1, " ", s2))
    match.y <- match.y + mapply(function(x, y) sum(x == y), strsplit(s1, ""), strsplit(s2, ""))
    if (match.y / count.y >= thresh & s1 == s2) {
      # max.y <- hbonds.y + hbond(s1)
      hbonds.y <- hbonds.y + hbond(s1)
      hbonds.evaluated.y <- hbonds.evaluated.y + hbond(s1)
      # max.y = count.y
    }
    count.y <- count.y + 1
    y.s <- paste0(y.s, s1)
    y.t <- paste0(y.t, s2)
  }
  if (debug) {
    print(paste0(
      substring(y.s, 1, 12),
      "  ",
      substring(y.t, 1, 12)
    ))
  }
  if (debug) print(paste0(hbonds.x, " vs. ", hbonds.y))
  return(max(hbonds.x, hbonds.y))
}

#' homology 503bp
dt1 <- data.table(
  position = 1:nchar(seq1),
  deletion_size = 1:nchar(seq1)
) %>%
  expand.grid() %>%
  filter(position + deletion_size < nchar(seq1)) %>%
  as.data.table()

dt1$homology <- pbmclapply(1:nrow(dt1), function(i) {
  row <- dt1[i, ]
  x <- find_homology(row$position, row$deletion_size, seq1)
  return(x)
}, mc.cores = 40) %>% unlist()

dt1$homology.hbonds <- mclapply(1:nrow(dt1), function(i) {
  row <- dt1[i, ]
  x <- find_homology.hbonds(row$position, row$deletion_size, seq1)
  return(x)
}, mc.cores = 40) %>% unlist()



dt1$homeology0.75 <- pbmclapply(1:nrow(dt1), function(i) {
  row <- dt1[i, ]
  x <- find_homeology(row$position, row$deletion_size, seq1, thresh = 0.75)
  return(x)
}, mc.cores = 80) %>% unlist()

dt1$homeology0.8 <- pbmclapply(1:nrow(dt1), function(i) {
  row <- dt1[i, ]
  x <- find_homeology(row$position, row$deletion_size, seq1, thresh = 0.80)
  return(x)
}, mc.cores = 80) %>% unlist()

dt1$homeology0.5 <- pbmclapply(1:nrow(dt1), function(i) {
  row <- dt1[i, ]
  x <- find_homeology(row$position, row$deletion_size, seq1, thresh = 0.50)
  return(x)
}, mc.cores = 100) %>% unlist()


#' homeology 245bp
dt2 <- data.table(
  position = 1:nchar(seq2),
  deletion_size = 1:nchar(seq2)
) %>%
  expand.grid() %>%
  filter(position + deletion_size < nchar(seq2))

dt2$homology <- pbmclapply(1:nrow(dt2), function(i) {
  row <- dt2[i, ]
  x <- find_homology(row$position, row$deletion_size, seq2)
  return(x)
}, mc.cores = 40) %>% unlist()

dt2$homeology <- pbmclapply(1:nrow(dt2), function(i) {
  row <- dt2[i, ]
  x <- find_homeology(row$position, row$deletion_size, seq2, thresh = 0.50)
  return(x)
}, mc.cores = 80) %>% unlist()


#' jrafailov Monday, Oct 07, 2024 11:07:47 AM
#saveRDS(dt2, "/gpfs/commons/home/sbrylka/Projects/homeology/jr_stash/final_prediction_matrix.SMALL.50thres.rds")
#saveRDS(dt1, "/gpfs/commons/home/sbrylka/Projects/homeology/jr_stash/final_prediction_matrix.LARGE.50thres.rds")

dt2 <- readRDS"./data/final_prediction_matrix.SMALL.rds")
dt1 <- readRDS("./data/final_prediction_matrix.LARGE.rds")

small_data_matrix <- dt2 %>%
  transmute(
    position = position - 70,
    deletion_size,
    homology = case_when(
      deletion_size < abs(position) ~ NA,
      T ~ homology
    ),
    homeology = case_when(
      deletion_size < abs(position) ~ NA,
      T ~ homeology
    )
  ) %>%
  filter(
    !is.na(homology),
    position <= 0,
    deletion_size <= 100
  ) %>%
  as.data.table()


medium_data_matrix <- dt1 %>%
  transmute(
    position = position - 252,
    deletion_size,
    homology = case_when(
      deletion_size < abs(position) ~ NA,
      T ~ homology
    ),
    homeology = case_when(
      deletion_size < abs(position) ~ NA,
      T ~ homeology
    )
  ) %>%
  filter(
    position <= 0,
    deletion_size <= 200,
    position >= -200
  )

large_data_matrix <- dt1 %>%
  transmute(
    position = position - 252,
    deletion_size,
    homology = case_when(
      deletion_size < abs(position) ~ NA,
      T ~ homology
    ),
    homeology = case_when(
      deletion_size < abs(position) ~ NA,
      T ~ homeology
    )
  ) %>%
  filter(
    position <= 0,
    deletion_size <= 251,
    position >= -251
  )

#' small homology
small.homology <- small_data_matrix %>%
  dplyr::select(-c("homeology")) %>%
  pivot_wider(
    names_from = "deletion_size",
    values_from = "homology",
    names_sort = T
  ) %>%
  as.data.table() %>%
  setkey(position) %>%
  as.matrix(rownames = T)

#' small homeology
small.homeology <- small_data_matrix %>%
  dplyr::select(-c("homology")) %>%
  pivot_wider(
    names_from = "deletion_size",
    values_from = "homeology",
    names_sort = T
  ) %>%
  as.data.table() %>%
  setkey(position) %>%
  as.matrix(rownames = T)

#' medium homology
medium.homology <- medium_data_matrix %>%
  dplyr::select(-c("homeology")) %>%
  pivot_wider(
    names_from = "deletion_size",
    values_from = "homology",
    names_sort = T
  ) %>%
  as.data.table() %>%
  setkey(position) %>%
  as.matrix(rownames = T)
#' large homeology
medium.homeology <- medium_data_matrix %>%
  dplyr::select(-c("homology")) %>%
  pivot_wider(
    names_from = "deletion_size",
    values_from = "homeology",
    names_sort = T
  ) %>%
  as.data.table() %>%
  setkey(position) %>%
  as.matrix(rownames = T)

#' large homology
large.homology <- large_data_matrix %>%
  dplyr::select(-c("homeology")) %>%
  pivot_wider(
    names_from = "deletion_size",
    values_from = "homology",
    names_sort = T
  ) %>%
  as.data.table() %>%
  setkey(position) %>%
  as.matrix(rownames = T)
#' large homeology
large.homeology <- large_data_matrix %>%
  dplyr::select(-c("homology")) %>%
  pivot_wider(
    names_from = "deletion_size",
    values_from = "homeology",
    names_sort = T
  ) %>%
  as.data.table() %>%
  setkey(position) %>%
  as.matrix(rownames = T)

small.homology.breaks <- c(0:max(small.homology, na.rm = T))
small.homeology.breaks <- c(0:max(small.homeology, na.rm = T))

medium.homology.breaks <- c(0:max(medium.homology, na.rm = T))
medium.homeology.breaks <- c(0:max(medium.homeology, na.rm = T))

large.homology.breaks <- c(0:max(large.homology, na.rm = T))
large.homeology.breaks <- c(0:max(large.homeology, na.rm = T))

small.homology.breaks
small.homeology.breaks
medium.homology.breaks
medium.homeology.breaks
large.homology.breaks
large.homeology.breaks

n.breaks <- length(large.homology.breaks) - 1
breaks.left <- length(large.homeology.breaks) - 1 - n.breaks
total.breaks <- length(large.homeology.breaks)

n.breaks + breaks.left - 1

color0 <- "#f9f9f9"

color0 <- "#e3e3e3"
col <- rev(sequential_hcl(n.breaks + 1, palette = "Mako")[-c(1)])
col_fun_prop <- colorRamp2(large.homeology.breaks, c(color0, color0, col), space = "RGB")

color0 <- "#e3e3e3"

## color0 = "#f9f9f9"
## darkgray = "#e3e3e3"

## ### two gradients

## rev(sequential_hcl(length(small.homeology.breaks) + 3, palette = "Mako"))

col <- rev(sequential_hcl(n.breaks + 3, palette = "Mako")[-c(1, 2, n.breaks + 3)])
col2 <- colorRamp2(
  breaks = c(
    max(large.homeology.breaks) - breaks.left,
    max(large.homeology.breaks)
  ),
  colors = c(col[length(col)], "palevioletred")
)
col2 <- col2(c((max(large.homeology.breaks) - breaks.left + 1):
max(large.homeology.breaks)))
col_fun_prop <- colorRamp2(large.homeology.breaks, c(color0, color0, col, col2), space = "RGB")

#### single gradient
col <- rev(sequential_hcl(n.breaks + 3, palette = "Mako")[-c(1, 2, n.breaks + 3)])

col <- rev(sequential_hcl(length(small.homeology.breaks) - 2, palette = "Mako"))
col_fun_prop <- colorRamp2(small.homeology.breaks, c(color0, color0, col), space = "RGB")

col <- rev(sequential_hcl(n.breaks + 2, palette = "Mako")[-c(1, 2, n.breaks + 2)])
col <- rev(sequential_hcl(length(large.homeology.breaks) - 1, palette = "Mako"))

col_fun_prop <- colorRamp2(large.homeology.breaks, c("white", col), space = "RGB")

breaks.left <- length(large.homeology.breaks)

## #### requested gradient
col1 <- rev(sequential_hcl(breaks.left + 3, palette = "Mako")[-c(1, breaks.left + 3)])[1]
col3 <- rev(sequential_hcl(breaks.left + 3, palette = "Mako")[-c(1, breaks.left + 3)])[c(2:(1 + breaks.left))]
col_fun_prop <- colorRamp2(large.homeology.breaks, c(color0, color0, rep(col1, n.breaks), col3), space = "RGB")






saveRDS(col_fun_prop, "~/projects/polq/jr_stash/col_fun_prop/ALL_round2.rds")

col_labels <- function(matrix) {
  colnames <- as.numeric(colnames(matrix))
  pos_of_min <- which(abs(colnames) == min(abs(colnames)))
  pos_of_5_mod <- which(colnames %% 5 == 0)
  keep <- c(pos_of_min, pos_of_5_mod)
  colnames[-keep] <- ""
  return(colnames)
}
row_labels <- function(matrix) {
  rownames <- as.numeric(rownames(matrix))
  pos_of_max <- which(abs(rownames) == max(abs(rownames)))
  pos_of_5_mod <- which(rownames %% 5 == 0)
  keep <- c(pos_of_max, pos_of_5_mod)
  rownames[-keep] <- ""
  return(rownames)
}


pdf(file = paste0("~/public_html/plots/homology/indel_heatmaps/", "HOMOLOGY_NEW_small_amplicon_prediction.SINGLEGRADIENT.pdf"), width = unit(4, "in"), height = unit(3.5, "in"))
lgd <- Legend(
  at = small.homology.breaks,
  # col_fun = col_fun_prop,
  legend_gp = gpar(fill = col_fun_prop(small.homology.breaks)),
  # breaks = c(0:11),
  # at = c(0,  1e-6 - 1e-7,5e-6,1e-5,5e-5,1e-4,5e-4,1e-3,5e-3,1e-2,5e-2,0.1,0.5),
  # labels = c("0",  "1e-6","5e-6","1e-5","5e-5","1e-4","5e-4","1e-3","5e-3","1e-2","5e-2","0.1","0.5"),
  title = paste0("Homology"),
  title_position = "lefttop-rot",
  row_gap = unit(0, "cm"),
  grid_width = unit(3, "mm"),
  grid_height = unit(3, "mm"),
  column_gap = unit(0, "cm"),
  #  legend_height = unit(1, "in"),
  title_gp = gpar(fontsize = 8, fontface = 4),
  labels_gp = gpar(fontsize = 6)
)

Heatmap(
  small.homology,
  row_order = rev(rownames(small.homology)),
  row_title = "Distance from cut site",
  row_title_side = "right",
  row_title_gp = gpar(fontsize = 8),
  col = col_fun_prop,
  rect_gp = gpar(col = "white", lwd = 0.001),
  # top_annotation = deletion_histogram,
  # right_annotation = row_histogram,
  column_title_side = "top",
  column_title = "Deletion size",
  column_title_gp = gpar(fontsize = 8),
  column_order = colnames(small.homology),
  show_column_names = T,
  show_row_names = T,
  column_labels = col_labels(small.homology),
  row_labels = row_labels(small.homology),
  column_names_side = "top",
  na_col = "white",
  width = ncol(small.homology) * unit(0.8, "mm"),
  height = nrow(small.homology) * unit(0.8, "mm"),
  show_heatmap_legend = F,
  row_names_gp = gpar(fontsize = 5, fontfamily = "mono", fontface = "bold"),
  column_names_gp = gpar(fontsize = 5, fontfamily = "mono", fontface = "bold")
) -> hm

draw(hm)
draw(lgd, x = unit(0.75, "in"), y = unit(1.2, "in"))
dev.off()



pdf(file = paste0("~/public_html/plots/homology/indel_heatmaps/", "HOMEOLOGY_NEW_small_amplicon_prediction.SINGLEGRADIENT.50percentthresh.pdf"), width = unit(4, "in"), height = unit(3.5, "in"))
lgd <- Legend(
  at = small.homeology.breaks,
  # col_fun = col_fun_prop,
  legend_gp = gpar(fill = col_fun_prop(small.homeology.breaks)),
  # breaks = c(0:11),
  # at = c(0,  1e-6 - 1e-7,5e-6,1e-5,5e-5,1e-4,5e-4,1e-3,5e-3,1e-2,5e-2,0.1,0.5),
  # labels = c("0",  "1e-6","5e-6","1e-5","5e-5","1e-4","5e-4","1e-3","5e-3","1e-2","5e-2","0.1","0.5"),
  title = paste0("Homeology"),
  title_position = "lefttop-rot",
  row_gap = unit(0, "cm"),
  grid_width = unit(3, "mm"),
  grid_height = unit(3, "mm"),
  column_gap = unit(0, "cm"),
  #  legend_height = unit(1, "in"),
  title_gp = gpar(fontsize = 8, fontface = 4),
  labels_gp = gpar(fontsize = 6)
)
Heatmap(
  small.homeology,
  row_order = rev(rownames(small.homeology)),
  row_title = "Distance from cut site",
  row_title_side = "right",
  row_title_gp = gpar(fontsize = 8),
  col = col_fun_prop,
  rect_gp = gpar(col = "white", lwd = 0.001),
  # top_annotation = deletion_histogram,
  # right_annotation = row_histogram,
  column_title_side = "top",
  column_title = "Deletion size",
  column_title_gp = gpar(fontsize = 8),
  column_order = colnames(small.homeology),
  show_column_names = T,
  show_row_names = T,
  column_labels = col_labels(small.homeology),
  row_labels = row_labels(small.homeology),
  column_names_side = "top",
  na_col = "white",
  width = ncol(small.homeology) * unit(0.8, "mm"),
  height = nrow(small.homeology) * unit(0.8, "mm"),
  show_heatmap_legend = F,
  row_names_gp = gpar(fontsize = 5, fontfamily = "mono", fontface = "bold"),
  column_names_gp = gpar(fontsize = 5, fontfamily = "mono", fontface = "bold")
) -> hm
draw(hm)
draw(lgd, x = unit(0.75, "in"), y = unit(1.2, "in"))
dev.off()

pdf(file = paste0("~/public_html/plots/homology/indel_heatmaps/", "HOMOLOGY_NEW_small_amplicon_prediction.REQUESTEDGRADIENT.pdf"), width = unit(4, "in"), height = unit(3.5, "in"))
lgd <- Legend(
  at = small.homology.breaks,
  # col_fun = col_fun_prop,
  legend_gp = gpar(fill = col_fun_prop(small.homology.breaks)),
  # breaks = c(0:11),
  # at = c(0,  1e-6 - 1e-7,5e-6,1e-5,5e-5,1e-4,5e-4,1e-3,5e-3,1e-2,5e-2,0.1,0.5),
  # labels = c("0",  "1e-6","5e-6","1e-5","5e-5","1e-4","5e-4","1e-3","5e-3","1e-2","5e-2","0.1","0.5"),
  title = paste0("Homology"),
  title_position = "lefttop-rot",
  row_gap = unit(0, "cm"),
  grid_width = unit(3, "mm"),
  grid_height = unit(3, "mm"),
  column_gap = unit(0, "cm"),
  #  legend_height = unit(1, "in"),
  title_gp = gpar(fontsize = 8, fontface = 4),
  labels_gp = gpar(fontsize = 6)
)
Heatmap(
  small.homology,
  row_order = rev(rownames(small.homology)),
  row_title = "Distance from cut site",
  row_title_side = "right",
  row_title_gp = gpar(fontsize = 8),
  col = col_fun_prop,
  rect_gp = gpar(col = "white", lwd = 0.001),
  # top_annotation = deletion_histogram,
  # right_annotation = row_histogram,
  column_title_side = "top",
  column_title = "Deletion size",
  column_title_gp = gpar(fontsize = 8),
  column_order = colnames(small.homology),
  show_column_names = T,
  show_row_names = T,
  column_labels = col_labels(small.homology),
  row_labels = row_labels(small.homology),
  column_names_side = "top",
  na_col = "white",
  width = ncol(small.homology) * unit(0.8, "mm"),
  height = nrow(small.homology) * unit(0.8, "mm"),
  show_heatmap_legend = F,
  row_names_gp = gpar(fontsize = 5, fontfamily = "mono", fontface = "bold"),
  column_names_gp = gpar(fontsize = 5, fontfamily = "mono", fontface = "bold")
) -> hm
draw(hm)
draw(lgd, x = unit(0.75, "in"), y = unit(1.2, "in"))
dev.off()



pdf(file = paste0("~/public_html/plots/homology/indel_heatmaps/", "HOMEOLOGY_NEW_small_amplicon_prediction.REQUESTEDGRADIENT.75percentthresh.pdf"), width = unit(4, "in"), height = unit(3.5, "in"))
lgd <- Legend(
  at = small.homeology.breaks,
  # col_fun = col_fun_prop,
  legend_gp = gpar(fill = col_fun_prop(small.homeology.breaks)),
  # breaks = c(0:11),
  # at = c(0,  1e-6 - 1e-7,5e-6,1e-5,5e-5,1e-4,5e-4,1e-3,5e-3,1e-2,5e-2,0.1,0.5),
  # labels = c("0",  "1e-6","5e-6","1e-5","5e-5","1e-4","5e-4","1e-3","5e-3","1e-2","5e-2","0.1","0.5"),
  title = paste0("Homeology"),
  title_position = "lefttop-rot",
  row_gap = unit(0, "cm"),
  grid_width = unit(3, "mm"),
  grid_height = unit(3, "mm"),
  column_gap = unit(0, "cm"),
  #  legend_height = unit(1, "in"),
  title_gp = gpar(fontsize = 8, fontface = 4),
  labels_gp = gpar(fontsize = 6)
)
Heatmap(
  small.homeology,
  row_order = rev(rownames(small.homeology)),
  row_title = "Distance from cut site",
  row_title_side = "right",
  row_title_gp = gpar(fontsize = 8),
  col = col_fun_prop,
  rect_gp = gpar(col = "white", lwd = 0.001),
  # top_annotation = deletion_histogram,
  # right_annotation = row_histogram,
  column_title_side = "top",
  column_title = "Deletion size",
  column_title_gp = gpar(fontsize = 8),
  column_order = colnames(small.homeology),
  show_column_names = T,
  show_row_names = T,
  column_labels = col_labels(small.homeology),
  row_labels = row_labels(small.homeology),
  column_names_side = "top",
  na_col = "white",
  width = ncol(small.homeology) * unit(0.8, "mm"),
  height = nrow(small.homeology) * unit(0.8, "mm"),
  show_heatmap_legend = F,
  row_names_gp = gpar(fontsize = 5, fontfamily = "mono", fontface = "bold"),
  column_names_gp = gpar(fontsize = 5, fontfamily = "mono", fontface = "bold")
) -> hm
draw(hm)
draw(lgd, x = unit(0.75, "in"), y = unit(1.2, "in"))
dev.off()

pdf(file = paste0("~/public_html/plots/homology/indel_heatmaps/", "HOMOLOGY_NEW_small_amplicon_prediction.TWOGRADIENTS.pdf"), width = unit(4, "in"), height = unit(3.5, "in"))
lgd <- Legend(
  at = small.homology.breaks,
  # col_fun = col_fun_prop,
  legend_gp = gpar(fill = col_fun_prop(small.homology.breaks)),
  # breaks = c(0:11),
  # at = c(0,  1e-6 - 1e-7,5e-6,1e-5,5e-5,1e-4,5e-4,1e-3,5e-3,1e-2,5e-2,0.1,0.5),
  # labels = c("0",  "1e-6","5e-6","1e-5","5e-5","1e-4","5e-4","1e-3","5e-3","1e-2","5e-2","0.1","0.5"),
  title = paste0("Homology"),
  title_position = "lefttop-rot",
  row_gap = unit(0, "cm"),
  grid_width = unit(3, "mm"),
  grid_height = unit(3, "mm"),
  column_gap = unit(0, "cm"),
  #  legend_height = unit(1, "in"),
  title_gp = gpar(fontsize = 8, fontface = 4),
  labels_gp = gpar(fontsize = 6)
)
Heatmap(
  small.homology,
  row_order = rev(rownames(small.homology)),
  row_title = "Distance from cut site",
  row_title_side = "right",
  row_title_gp = gpar(fontsize = 8),
  col = col_fun_prop,
  rect_gp = gpar(col = "white", lwd = 0.001),
  # top_annotation = deletion_histogram,
  # right_annotation = row_histogram,
  column_title_side = "top",
  column_title = "Deletion size",
  column_title_gp = gpar(fontsize = 8),
  column_order = colnames(small.homology),
  show_column_names = T,
  show_row_names = T,
  column_labels = col_labels(small.homology),
  row_labels = row_labels(small.homology),
  column_names_side = "top",
  na_col = "white",
  width = ncol(small.homology) * unit(0.8, "mm"),
  height = nrow(small.homology) * unit(0.8, "mm"),
  show_heatmap_legend = F,
  row_names_gp = gpar(fontsize = 5, fontfamily = "mono", fontface = "bold"),
  column_names_gp = gpar(fontsize = 5, fontfamily = "mono", fontface = "bold")
) -> hm
draw(hm)
draw(lgd, x = unit(0.75, "in"), y = unit(1.2, "in"))
dev.off()



pdf(file = paste0("~/public_html/plots/homology/indel_heatmaps/", "HOMEOLOGY_NEW_small_amplicon_prediction.TWOGRADIENTS.pdf"), width = unit(4, "in"), height = unit(3.5, "in"))
lgd <- Legend(
  at = small.homeology.breaks,
  # col_fun = col_fun_prop,
  legend_gp = gpar(fill = col_fun_prop(small.homeology.breaks)),
  # breaks = c(0:11),
  # at = c(0,  1e-6 - 1e-7,5e-6,1e-5,5e-5,1e-4,5e-4,1e-3,5e-3,1e-2,5e-2,0.1,0.5),
  # labels = c("0",  "1e-6","5e-6","1e-5","5e-5","1e-4","5e-4","1e-3","5e-3","1e-2","5e-2","0.1","0.5"),
  title = paste0("Homeology"),
  title_position = "lefttop-rot",
  row_gap = unit(0, "cm"),
  grid_width = unit(3, "mm"),
  grid_height = unit(3, "mm"),
  column_gap = unit(0, "cm"),
  #  legend_height = unit(1, "in"),
  title_gp = gpar(fontsize = 8, fontface = 4),
  labels_gp = gpar(fontsize = 6)
)
Heatmap(
  small.homeology,
  row_order = rev(rownames(small.homeology)),
  row_title = "Distance from cut site",
  row_title_side = "right",
  row_title_gp = gpar(fontsize = 8),
  col = col_fun_prop,
  rect_gp = gpar(col = "white", lwd = 0.001),
  # top_annotation = deletion_histogram,
  # right_annotation = row_histogram,
  column_title_side = "top",
  column_title = "Deletion size",
  column_title_gp = gpar(fontsize = 8),
  column_order = colnames(small.homeology),
  show_column_names = T,
  show_row_names = T,
  column_labels = col_labels(small.homeology),
  row_labels = row_labels(small.homeology),
  column_names_side = "top",
  na_col = "white",
  width = ncol(small.homeology) * unit(0.8, "mm"),
  height = nrow(small.homeology) * unit(0.8, "mm"),
  show_heatmap_legend = F,
  row_names_gp = gpar(fontsize = 5, fontfamily = "mono", fontface = "bold"),
  column_names_gp = gpar(fontsize = 5, fontfamily = "mono", fontface = "bold")
) -> hm
draw(hm)
draw(lgd, x = unit(0.75, "in"), y = unit(1.2, "in"))
dev.off()



col_labels <- function(matrix) {
  colnames <- as.numeric(colnames(matrix))
  pos_of_min <- which(abs(colnames) == min(abs(colnames)))
  pos_of_mod <- which(colnames %% 10 == 0)
  keep <- c(pos_of_min, pos_of_mod)
  colnames[-keep] <- ""
  return(colnames)
}
row_labels <- function(matrix) {
  rownames <- as.numeric(rownames(matrix))
  pos_of_max <- which(abs(rownames) == max(abs(rownames)))
  pos_of_mod <- which(rownames %% 10 == 0)
  keep <- c(pos_of_max, pos_of_mod)
  rownames[-keep] <- ""
  return(rownames)
}

pdf(file = paste0("~/public_html/plots/homology/indel_heatmaps/", "HOMOLOGY_NEW_medium_amplicon_prediction.SINGLEGRADIENT.pdf"), width = unit(4.5, "in"), height = unit(4.5, "in"))
lgd <- Legend(
  at = medium.homology.breaks,
  # col_fun = col_fun_prop,
  legend_gp = gpar(fill = col_fun_prop(medium.homology.breaks)),
  # breaks = c(0:11),
  # at = c(0,  1e-6 - 1e-7,5e-6,1e-5,5e-5,1e-4,5e-4,1e-3,5e-3,1e-2,5e-2,0.1,0.5),
  # labels = c("0",  "1e-6","5e-6","1e-5","5e-5","1e-4","5e-4","1e-3","5e-3","1e-2","5e-2","0.1","0.5"),
  title = paste0("Homology"),
  title_position = "lefttop-rot",
  row_gap = unit(0, "cm"),
  grid_width = unit(5, "mm"),
  grid_height = unit(5, "mm"),
  column_gap = unit(0, "cm"),
  #  legend_height = unit(1, "in"),
  title_gp = gpar(fontsize = 12, fontface = 4),
  labels_gp = gpar(fontsize = 10)
)
Heatmap(
  medium.homology,
  row_order = rev(rownames(medium.homology)),
  row_title = "Distance from cut site",
  row_title_side = "right",
  row_title_gp = gpar(fontsize = 8),
  col = col_fun_prop,
  rect_gp = gpar(col = "white", lwd = 0.001),
  # top_annotation = deletion_histogram,
  # right_annotation = row_histogram,
  column_title_side = "top",
  column_title = "Deletion size",
  column_title_gp = gpar(fontsize = 8),
  column_order = colnames(medium.homology),
  show_column_names = T,
  show_row_names = T,
  column_labels = col_labels(medium.homology),
  row_labels = row_labels(medium.homology),
  column_names_side = "top",
  na_col = "white",
  width = ncol(medium.homology) * unit(0.45, "mm"),
  height = nrow(medium.homology) * unit(0.45, "mm"),
  show_heatmap_legend = F,
  row_names_gp = gpar(fontsize = 5, fontfamily = "mono", fontface = "bold"),
  column_names_gp = gpar(fontsize = 5, fontfamily = "mono", fontface = "bold")
) -> hm
draw(hm)
draw(lgd, x = unit(0.6, "in"), y = unit(2, "in"))
dev.off()

pdf(file = paste0("~/public_html/plots/homology/indel_heatmaps/", "HOMEOLOGY_NEW_medium_amplicon_prediction.REQUESTEDGRADIENT.pdf"), width = unit(4.5, "in"), height = unit(4.5, "in"))
lgd <- Legend(
  at = medium.homeology.breaks,
  # col_fun = col_fun_prop,
  legend_gp = gpar(fill = col_fun_prop(medium.homeology.breaks)),
  # breaks = c(0:11),
  # at = c(0,  1e-6 - 1e-7,5e-6,1e-5,5e-5,1e-4,5e-4,1e-3,5e-3,1e-2,5e-2,0.1,0.5),
  # labels = c("0",  "1e-6","5e-6","1e-5","5e-5","1e-4","5e-4","1e-3","5e-3","1e-2","5e-2","0.1","0.5"),
  title = paste0("Homeology"),
  title_position = "lefttop-rot",
  row_gap = unit(0, "cm"),
  grid_width = unit(5, "mm"),
  grid_height = unit(5, "mm"),
  column_gap = unit(0, "cm"),
  #  legend_height = unit(1, "in"),
  title_gp = gpar(fontsize = 12, fontface = 4),
  labels_gp = gpar(fontsize = 10)
)
Heatmap(
  medium.homeology,
  row_order = rev(rownames(medium.homeology)),
  row_title = "Distance from cut site",
  row_title_side = "right",
  row_title_gp = gpar(fontsize = 8),
  col = col_fun_prop,
  rect_gp = gpar(col = "white", lwd = 0.001),
  # top_annotation = deletion_histogram,
  # right_annotation = row_histogram,
  column_title_side = "top",
  column_title = "Deletion size",
  column_title_gp = gpar(fontsize = 8),
  column_order = colnames(medium.homeology),
  show_column_names = T,
  show_row_names = T,
  column_labels = col_labels(medium.homeology),
  row_labels = row_labels(medium.homeology),
  column_names_side = "top",
  na_col = "white",
  width = ncol(medium.homeology) * unit(0.45, "mm"),
  height = nrow(medium.homeology) * unit(0.45, "mm"),
  show_heatmap_legend = F,
  row_names_gp = gpar(fontsize = 5, fontfamily = "mono", fontface = "bold"),
  column_names_gp = gpar(fontsize = 5, fontfamily = "mono", fontface = "bold")
) -> hm
draw(hm)
draw(lgd, x = unit(0.6, "in"), y = unit(1.75, "in"))
dev.off()

pdf(file = paste0("~/public_html/plots/homology/indel_heatmaps/", "HOMOLOGY_NEW_medium_amplicon_prediction.REQUESTEDGRADIENT.pdf"), width = unit(4.5, "in"), height = unit(4.5, "in"))
lgd <- Legend(
  at = medium.homology.breaks,
  # col_fun = col_fun_prop,
  legend_gp = gpar(fill = col_fun_prop(medium.homology.breaks)),
  # breaks = c(0:11),
  # at = c(0,  1e-6 - 1e-7,5e-6,1e-5,5e-5,1e-4,5e-4,1e-3,5e-3,1e-2,5e-2,0.1,0.5),
  # labels = c("0",  "1e-6","5e-6","1e-5","5e-5","1e-4","5e-4","1e-3","5e-3","1e-2","5e-2","0.1","0.5"),
  title = paste0("Homology"),
  title_position = "lefttop-rot",
  row_gap = unit(0, "cm"),
  grid_width = unit(5, "mm"),
  grid_height = unit(5, "mm"),
  column_gap = unit(0, "cm"),
  #  legend_height = unit(1, "in"),
  title_gp = gpar(fontsize = 12, fontface = 4),
  labels_gp = gpar(fontsize = 10)
)
Heatmap(
  medium.homology,
  row_order = rev(rownames(medium.homology)),
  row_title = "Distance from cut site",
  row_title_side = "right",
  row_title_gp = gpar(fontsize = 8),
  col = col_fun_prop,
  rect_gp = gpar(col = "white", lwd = 0.001),
  # top_annotation = deletion_histogram,
  # right_annotation = row_histogram,
  column_title_side = "top",
  column_title = "Deletion size",
  column_title_gp = gpar(fontsize = 8),
  column_order = colnames(medium.homology),
  show_column_names = T,
  show_row_names = T,
  column_labels = col_labels(medium.homology),
  row_labels = row_labels(medium.homology),
  column_names_side = "top",
  na_col = "white",
  width = ncol(medium.homology) * unit(0.45, "mm"),
  height = nrow(medium.homology) * unit(0.45, "mm"),
  show_heatmap_legend = F,
  row_names_gp = gpar(fontsize = 5, fontfamily = "mono", fontface = "bold"),
  column_names_gp = gpar(fontsize = 5, fontfamily = "mono", fontface = "bold")
) -> hm
draw(hm)
draw(lgd, x = unit(0.6, "in"), y = unit(2, "in"))
dev.off()

pdf(file = paste0("~/public_html/plots/homology/indel_heatmaps/", "HOMEOLOGY_NEW_medium_amplicon_prediction.SINGLEGRADIENT.50percent.pdf"), width = unit(4.5, "in"), height = unit(4.5, "in"))
lgd <- Legend(
  at = medium.homeology.breaks,
  # col_fun = col_fun_prop,
  legend_gp = gpar(fill = col_fun_prop(medium.homeology.breaks)),
  # breaks = c(0:11),
  # at = c(0,  1e-6 - 1e-7,5e-6,1e-5,5e-5,1e-4,5e-4,1e-3,5e-3,1e-2,5e-2,0.1,0.5),
  # labels = c("0",  "1e-6","5e-6","1e-5","5e-5","1e-4","5e-4","1e-3","5e-3","1e-2","5e-2","0.1","0.5"),
  title = paste0("Homeology"),
  title_position = "lefttop-rot",
  row_gap = unit(0, "cm"),
  grid_width = unit(5, "mm"),
  grid_height = unit(5, "mm"),
  column_gap = unit(0, "cm"),
  #  legend_height = unit(1, "in"),
  title_gp = gpar(fontsize = 12, fontface = 4),
  labels_gp = gpar(fontsize = 10)
)
Heatmap(
  medium.homeology,
  row_order = rev(rownames(medium.homeology)),
  row_title = "Distance from cut site",
  row_title_side = "right",
  row_title_gp = gpar(fontsize = 8),
  col = col_fun_prop,
  rect_gp = gpar(col = "white", lwd = 0.001),
  # top_annotation = deletion_histogram,
  # right_annotation = row_histogram,
  column_title_side = "top",
  column_title = "Deletion size",
  column_title_gp = gpar(fontsize = 8),
  column_order = colnames(medium.homeology),
  show_column_names = T,
  show_row_names = T,
  column_labels = col_labels(medium.homeology),
  row_labels = row_labels(medium.homeology),
  column_names_side = "top",
  na_col = "white",
  width = ncol(medium.homeology) * unit(0.45, "mm"),
  height = nrow(medium.homeology) * unit(0.45, "mm"),
  show_heatmap_legend = F,
  row_names_gp = gpar(fontsize = 5, fontfamily = "mono", fontface = "bold"),
  column_names_gp = gpar(fontsize = 5, fontfamily = "mono", fontface = "bold")
) -> hm
draw(hm)
draw(lgd, x = unit(0.6, "in"), y = unit(1.75, "in"))
dev.off()

pdf(file = paste0("~/public_html/plots/homology/indel_heatmaps/", "HOMOLOGY_NEW_medium_amplicon_prediction.TWOGRADIENT.pdf"), width = unit(4.5, "in"), height = unit(4.5, "in"))
lgd <- Legend(
  at = medium.homology.breaks,
  # col_fun = col_fun_prop,
  legend_gp = gpar(fill = col_fun_prop(medium.homology.breaks)),
  # breaks = c(0:11),
  # at = c(0,  1e-6 - 1e-7,5e-6,1e-5,5e-5,1e-4,5e-4,1e-3,5e-3,1e-2,5e-2,0.1,0.5),
  # labels = c("0",  "1e-6","5e-6","1e-5","5e-5","1e-4","5e-4","1e-3","5e-3","1e-2","5e-2","0.1","0.5"),
  title = paste0("Homology"),
  title_position = "lefttop-rot",
  row_gap = unit(0, "cm"),
  grid_width = unit(5, "mm"),
  grid_height = unit(5, "mm"),
  column_gap = unit(0, "cm"),
  #  legend_height = unit(1, "in"),
  title_gp = gpar(fontsize = 12, fontface = 4),
  labels_gp = gpar(fontsize = 10)
)
Heatmap(
  medium.homology,
  row_order = rev(rownames(medium.homology)),
  row_title = "Distance from cut site",
  row_title_side = "right",
  row_title_gp = gpar(fontsize = 8),
  col = col_fun_prop,
  rect_gp = gpar(col = "white", lwd = 0.001),
  # top_annotation = deletion_histogram,
  # right_annotation = row_histogram,
  column_title_side = "top",
  column_title = "Deletion size",
  column_title_gp = gpar(fontsize = 8),
  column_order = colnames(medium.homology),
  show_column_names = T,
  show_row_names = T,
  column_labels = col_labels(medium.homology),
  row_labels = row_labels(medium.homology),
  column_names_side = "top",
  na_col = "white",
  width = ncol(medium.homology) * unit(0.45, "mm"),
  height = nrow(medium.homology) * unit(0.45, "mm"),
  show_heatmap_legend = F,
  row_names_gp = gpar(fontsize = 5, fontfamily = "mono", fontface = "bold"),
  column_names_gp = gpar(fontsize = 5, fontfamily = "mono", fontface = "bold")
) -> hm
draw(hm)
draw(lgd, x = unit(0.6, "in"), y = unit(2, "in"))
dev.off()

pdf(file = paste0("~/public_html/plots/homology/indel_heatmaps/", "HOMEOLOGY_NEW_medium_amplicon_prediction.TWOGRADIENT.pdf"), width = unit(4.5, "in"), height = unit(4.5, "in"))
lgd <- Legend(
  at = medium.homeology.breaks,
  # col_fun = col_fun_prop,
  legend_gp = gpar(fill = col_fun_prop(medium.homeology.breaks)),
  # breaks = c(0:11),
  # at = c(0,  1e-6 - 1e-7,5e-6,1e-5,5e-5,1e-4,5e-4,1e-3,5e-3,1e-2,5e-2,0.1,0.5),
  # labels = c("0",  "1e-6","5e-6","1e-5","5e-5","1e-4","5e-4","1e-3","5e-3","1e-2","5e-2","0.1","0.5"),
  title = paste0("Homeology"),
  title_position = "lefttop-rot",
  row_gap = unit(0, "cm"),
  grid_width = unit(5, "mm"),
  grid_height = unit(5, "mm"),
  column_gap = unit(0, "cm"),
  #  legend_height = unit(1, "in"),
  title_gp = gpar(fontsize = 12, fontface = 4),
  labels_gp = gpar(fontsize = 10)
)
Heatmap(
  medium.homeology,
  row_order = rev(rownames(medium.homeology)),
  row_title = "Distance from cut site",
  row_title_side = "right",
  row_title_gp = gpar(fontsize = 8),
  col = col_fun_prop,
  rect_gp = gpar(col = "white", lwd = 0.001),
  # top_annotation = deletion_histogram,
  # right_annotation = row_histogram,
  column_title_side = "top",
  column_title = "Deletion size",
  column_title_gp = gpar(fontsize = 8),
  column_order = colnames(medium.homeology),
  show_column_names = T,
  show_row_names = T,
  column_labels = col_labels(medium.homeology),
  row_labels = row_labels(medium.homeology),
  column_names_side = "top",
  na_col = "white",
  width = ncol(medium.homeology) * unit(0.45, "mm"),
  height = nrow(medium.homeology) * unit(0.45, "mm"),
  show_heatmap_legend = F,
  row_names_gp = gpar(fontsize = 5, fontfamily = "mono", fontface = "bold"),
  column_names_gp = gpar(fontsize = 5, fontfamily = "mono", fontface = "bold")
) -> hm
draw(hm)
draw(lgd, x = unit(0.6, "in"), y = unit(1.75, "in"))
dev.off()


pdf(file = paste0("~/public_html/plots/homology/indel_heatmaps/", "HOMOLOGY_NEW_large_amplicon_prediction.SINGLEGRADIENT.pdf"), width = unit(5.5, "in"), height = unit(5.5, "in"))
lgd <- Legend(
  at = large.homology.breaks,
  # col_fun = col_fun_prop,
  legend_gp = gpar(fill = col_fun_prop(large.homology.breaks)),
  # breaks = c(0:11),
  # at = c(0,  1e-6 - 1e-7,5e-6,1e-5,5e-5,1e-4,5e-4,1e-3,5e-3,1e-2,5e-2,0.1,0.5),
  # labels = c("0",  "1e-6","5e-6","1e-5","5e-5","1e-4","5e-4","1e-3","5e-3","1e-2","5e-2","0.1","0.5"),
  title = paste0("Homology"),
  title_position = "lefttop-rot",
  row_gap = unit(0, "cm"),
  grid_width = unit(5, "mm"),
  grid_height = unit(5, "mm"),
  column_gap = unit(0, "cm"),
  #  legend_height = unit(1, "in"),
  title_gp = gpar(fontsize = 12, fontface = 4),
  labels_gp = gpar(fontsize = 10)
)
Heatmap(
  large.homology,
  row_order = rev(rownames(large.homology)),
  row_title = "Distance from cut site",
  row_title_side = "right",
  row_title_gp = gpar(fontsize = 8),
  col = col_fun_prop,
  rect_gp = gpar(col = "white", lwd = 0.001),
  # top_annotation = deletion_histogram,
  # right_annotation = row_histogram,
  column_title_side = "top",
  column_title = "Deletion size",
  column_title_gp = gpar(fontsize = 8),
  column_order = colnames(large.homology),
  show_column_names = T,
  show_row_names = T,
  column_labels = col_labels(large.homology),
  row_labels = row_labels(large.homology),
  column_names_side = "top",
  na_col = "white",
  width = ncol(large.homology) * unit(0.45, "mm"),
  height = nrow(large.homology) * unit(0.45, "mm"),
  show_heatmap_legend = F,
  row_names_gp = gpar(fontsize = 5, fontfamily = "mono", fontface = "bold"),
  column_names_gp = gpar(fontsize = 5, fontfamily = "mono", fontface = "bold")
) -> hm
draw(hm)
draw(lgd, x = unit(0.6, "in"), y = unit(2.5, "in"))
dev.off()

pdf(file = paste0("~/public_html/plots/homology/indel_heatmaps/", "HOMEOLOGY_NEW_large_amplicon_prediction.SINGLEGRADIENT.50percent.pdf"), width = unit(5.5, "in"), height = unit(5.5, "in"))
lgd <- Legend(
  at = large.homeology.breaks,
  # col_fun = col_fun_prop,
  legend_gp = gpar(fill = col_fun_prop(large.homeology.breaks)),
  # breaks = c(0:11),
  # at = c(0,  1e-6 - 1e-7,5e-6,1e-5,5e-5,1e-4,5e-4,1e-3,5e-3,1e-2,5e-2,0.1,0.5),
  # labels = c("0",  "1e-6","5e-6","1e-5","5e-5","1e-4","5e-4","1e-3","5e-3","1e-2","5e-2","0.1","0.5"),
  title = paste0("Homeology"),
  title_position = "lefttop-rot",
  row_gap = unit(0, "cm"),
  grid_width = unit(5, "mm"),
  grid_height = unit(5, "mm"),
  column_gap = unit(0, "cm"),
  #  legend_height = unit(1, "in"),
  title_gp = gpar(fontsize = 12, fontface = 4),
  labels_gp = gpar(fontsize = 10)
)
Heatmap(
  large.homeology,
  row_order = rev(rownames(large.homeology)),
  row_title = "Distance from cut site",
  row_title_side = "right",
  row_title_gp = gpar(fontsize = 8),
  col = col_fun_prop,
  rect_gp = gpar(col = "white", lwd = 0.001),
  # top_annotation = deletion_histogram,
  # right_annotation = row_histogram,
  column_title_side = "top",
  column_title = "Deletion size",
  column_title_gp = gpar(fontsize = 8),
  column_order = colnames(large.homeology),
  show_column_names = T,
  show_row_names = T,
  column_labels = col_labels(large.homeology),
  row_labels = row_labels(large.homeology),
  column_names_side = "top",
  na_col = "white",
  width = ncol(large.homeology) * unit(0.45, "mm"),
  height = nrow(large.homeology) * unit(0.45, "mm"),
  show_heatmap_legend = F,
  row_names_gp = gpar(fontsize = 5, fontfamily = "mono", fontface = "bold"),
  column_names_gp = gpar(fontsize = 5, fontfamily = "mono", fontface = "bold")
) -> hm
draw(hm)
draw(lgd, x = unit(0.6, "in"), y = unit(2.5, "in"))
dev.off()


pdf(file = paste0("~/public_html/plots/homology/indel_heatmaps/", "HOMOLOGY_NEW_large_amplicon_prediction.REQUESTEDGRADIENT.pdf"), width = unit(5.5, "in"), height = unit(5.5, "in"))
lgd <- Legend(
  at = large.homology.breaks,
  # col_fun = col_fun_prop,
  legend_gp = gpar(fill = col_fun_prop(large.homology.breaks)),
  # breaks = c(0:11),
  # at = c(0,  1e-6 - 1e-7,5e-6,1e-5,5e-5,1e-4,5e-4,1e-3,5e-3,1e-2,5e-2,0.1,0.5),
  # labels = c("0",  "1e-6","5e-6","1e-5","5e-5","1e-4","5e-4","1e-3","5e-3","1e-2","5e-2","0.1","0.5"),
  title = paste0("Homology"),
  title_position = "lefttop-rot",
  row_gap = unit(0, "cm"),
  grid_width = unit(5, "mm"),
  grid_height = unit(5, "mm"),
  column_gap = unit(0, "cm"),
  #  legend_height = unit(1, "in"),
  title_gp = gpar(fontsize = 12, fontface = 4),
  labels_gp = gpar(fontsize = 10)
)
Heatmap(
  large.homology,
  row_order = rev(rownames(large.homology)),
  row_title = "Distance from cut site",
  row_title_side = "right",
  row_title_gp = gpar(fontsize = 8),
  col = col_fun_prop,
  rect_gp = gpar(col = "white", lwd = 0.001),
  # top_annotation = deletion_histogram,
  # right_annotation = row_histogram,
  column_title_side = "top",
  column_title = "Deletion size",
  column_title_gp = gpar(fontsize = 8),
  column_order = colnames(large.homology),
  show_column_names = T,
  show_row_names = T,
  column_labels = col_labels(large.homology),
  row_labels = row_labels(large.homology),
  column_names_side = "top",
  na_col = "white",
  width = ncol(large.homology) * unit(0.45, "mm"),
  height = nrow(large.homology) * unit(0.45, "mm"),
  show_heatmap_legend = F,
  row_names_gp = gpar(fontsize = 5, fontfamily = "mono", fontface = "bold"),
  column_names_gp = gpar(fontsize = 5, fontfamily = "mono", fontface = "bold")
) -> hm
draw(hm)
draw(lgd, x = unit(0.6, "in"), y = unit(2.5, "in"))
dev.off()

pdf(file = paste0("~/public_html/plots/homology/indel_heatmaps/", "HOMEOLOGY_NEW_large_amplicon_prediction.REQUESTEDGRADIENT.pdf"), width = unit(5.5, "in"), height = unit(5.5, "in"))
lgd <- Legend(
  at = large.homeology.breaks,
  # col_fun = col_fun_prop,
  legend_gp = gpar(fill = col_fun_prop(large.homeology.breaks)),
  # breaks = c(0:11),
  # at = c(0,  1e-6 - 1e-7,5e-6,1e-5,5e-5,1e-4,5e-4,1e-3,5e-3,1e-2,5e-2,0.1,0.5),
  # labels = c("0",  "1e-6","5e-6","1e-5","5e-5","1e-4","5e-4","1e-3","5e-3","1e-2","5e-2","0.1","0.5"),
  title = paste0("Homeology"),
  title_position = "lefttop-rot",
  row_gap = unit(0, "cm"),
  grid_width = unit(5, "mm"),
  grid_height = unit(5, "mm"),
  column_gap = unit(0, "cm"),
  #  legend_height = unit(1, "in"),
  title_gp = gpar(fontsize = 12, fontface = 4),
  labels_gp = gpar(fontsize = 10)
)
Heatmap(
  large.homeology,
  row_order = rev(rownames(large.homeology)),
  row_title = "Distance from cut site",
  row_title_side = "right",
  row_title_gp = gpar(fontsize = 8),
  col = col_fun_prop,
  rect_gp = gpar(col = "white", lwd = 0.001),
  # top_annotation = deletion_histogram,
  # right_annotation = row_histogram,
  column_title_side = "top",
  column_title = "Deletion size",
  column_title_gp = gpar(fontsize = 8),
  column_order = colnames(large.homeology),
  show_column_names = T,
  show_row_names = T,
  column_labels = col_labels(large.homeology),
  row_labels = row_labels(large.homeology),
  column_names_side = "top",
  na_col = "white",
  width = ncol(large.homeology) * unit(0.45, "mm"),
  height = nrow(large.homeology) * unit(0.45, "mm"),
  show_heatmap_legend = F,
  row_names_gp = gpar(fontsize = 5, fontfamily = "mono", fontface = "bold"),
  column_names_gp = gpar(fontsize = 5, fontfamily = "mono", fontface = "bold")
) -> hm
draw(hm)
draw(lgd, x = unit(0.6, "in"), y = unit(2.5, "in"))
dev.off()

pdf(file = paste0("~/public_html/plots/homology/indel_heatmaps/", "HOMOLOGY_NEW_large_amplicon_prediction.TWOGRADIENT.pdf"), width = unit(5.5, "in"), height = unit(5.5, "in"))
lgd <- Legend(
  at = large.homology.breaks,
  # col_fun = col_fun_prop,
  legend_gp = gpar(fill = col_fun_prop(large.homology.breaks)),
  # breaks = c(0:11),
  # at = c(0,  1e-6 - 1e-7,5e-6,1e-5,5e-5,1e-4,5e-4,1e-3,5e-3,1e-2,5e-2,0.1,0.5),
  # labels = c("0",  "1e-6","5e-6","1e-5","5e-5","1e-4","5e-4","1e-3","5e-3","1e-2","5e-2","0.1","0.5"),
  title = paste0("Homology"),
  title_position = "lefttop-rot",
  row_gap = unit(0, "cm"),
  grid_width = unit(5, "mm"),
  grid_height = unit(5, "mm"),
  column_gap = unit(0, "cm"),
  #  legend_height = unit(1, "in"),
  title_gp = gpar(fontsize = 12, fontface = 4),
  labels_gp = gpar(fontsize = 10)
)
Heatmap(
  large.homology,
  row_order = rev(rownames(large.homology)),
  row_title = "Distance from cut site",
  row_title_side = "right",
  row_title_gp = gpar(fontsize = 8),
  col = col_fun_prop,
  rect_gp = gpar(col = "white", lwd = 0.001),
  # top_annotation = deletion_histogram,
  # right_annotation = row_histogram,
  column_title_side = "top",
  column_title = "Deletion size",
  column_title_gp = gpar(fontsize = 8),
  column_order = colnames(large.homology),
  show_column_names = T,
  show_row_names = T,
  column_labels = col_labels(large.homology),
  row_labels = row_labels(large.homology),
  column_names_side = "top",
  na_col = "white",
  width = ncol(large.homology) * unit(0.45, "mm"),
  height = nrow(large.homology) * unit(0.45, "mm"),
  show_heatmap_legend = F,
  row_names_gp = gpar(fontsize = 5, fontfamily = "mono", fontface = "bold"),
  column_names_gp = gpar(fontsize = 5, fontfamily = "mono", fontface = "bold")
) -> hm
draw(hm)
draw(lgd, x = unit(0.6, "in"), y = unit(2.5, "in"))
dev.off()

pdf(file = paste0("~/public_html/plots/homology/indel_heatmaps/", "HOMEOLOGY_NEW_large_amplicon_prediction.TWOGRADIENT.pdf"), width = unit(5.5, "in"), height = unit(5.5, "in"))
lgd <- Legend(
  at = large.homeology.breaks,
  # col_fun = col_fun_prop,
  legend_gp = gpar(fill = col_fun_prop(large.homeology.breaks)),
  # breaks = c(0:11),
  # at = c(0,  1e-6 - 1e-7,5e-6,1e-5,5e-5,1e-4,5e-4,1e-3,5e-3,1e-2,5e-2,0.1,0.5),
  # labels = c("0",  "1e-6","5e-6","1e-5","5e-5","1e-4","5e-4","1e-3","5e-3","1e-2","5e-2","0.1","0.5"),
  title = paste0("Homeology"),
  title_position = "lefttop-rot",
  row_gap = unit(0, "cm"),
  grid_width = unit(5, "mm"),
  grid_height = unit(5, "mm"),
  column_gap = unit(0, "cm"),
  #  legend_height = unit(1, "in"),
  title_gp = gpar(fontsize = 12, fontface = 4),
  labels_gp = gpar(fontsize = 10)
)
Heatmap(
  large.homeology,
  row_order = rev(rownames(large.homeology)),
  row_title = "Distance from cut site",
  row_title_side = "right",
  row_title_gp = gpar(fontsize = 8),
  col = col_fun_prop,
  rect_gp = gpar(col = "white", lwd = 0.001),
  # top_annotation = deletion_histogram,
  # right_annotation = row_histogram,
  column_title_side = "top",
  column_title = "Deletion size",
  column_title_gp = gpar(fontsize = 8),
  column_order = colnames(large.homeology),
  show_column_names = T,
  show_row_names = T,
  column_labels = col_labels(large.homeology),
  row_labels = row_labels(large.homeology),
  column_names_side = "top",
  na_col = "white",
  width = ncol(large.homeology) * unit(0.45, "mm"),
  height = nrow(large.homeology) * unit(0.45, "mm"),
  show_heatmap_legend = F,
  row_names_gp = gpar(fontsize = 5, fontfamily = "mono", fontface = "bold"),
  column_names_gp = gpar(fontsize = 5, fontfamily = "mono", fontface = "bold")
) -> hm
draw(hm)
draw(lgd, x = unit(0.6, "in"), y = unit(2.5, "in"))
dev.off()
