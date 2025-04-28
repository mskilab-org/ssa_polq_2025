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
library(ggtext)
library(rstatix)
library(ggforce)
library(cowplot)
library(stringi)
library(scico)
library(scales)
library(khroma)
library(ggsci)
library(circlize)
library(ggnewscale)

setDTthreads(1)


#' Generate Tornado Plot for Variant Data
#'
#' This function processes a dataset containing variant (deletion, insertion, mutation) information,
#' groups and summarizes the data based on various criteria, and generates a tornado plot visualizing
#' indel density along with variant characteristics. It supports returning specific variant plots (e.g., only
#' insertions, only deletions, or only mutations) or combining them into one overall plot.
#'
#' @param data A data.table or data.frame containing variant data with fields such as "cell_line", "genotype",
#'             "deletion", "insertion_sizes", "substitution_positions", and others required for summarization.
#' @param plot.ins Logical. If TRUE, insertion variants will be processed and plotted.
#' @param plot.mut Logical. If TRUE, mutation variants will be processed and plotted.
#' @param x.spacing Numeric. Determines the spacing of breaks on the x-axis (sequence position).
#' @param adjust.size Numeric. A value used to adjust positional parameters for plotting (e.g., centering insertions).
#' @param gts A character vector specifying which genotypes to include in the analysis.
#' @param all Logical. If TRUE, all data is used to generate the main plot; if FALSE, processing may be specific to a variant type.
#' @param hom.min Integer. Minimum homology value to filter deletion events when returning deletion-only plots.
#' @param return.ins Logical. If TRUE, only the insertion plot is returned.
#' @param return.del Logical. If TRUE, only the deletion plot is returned.
#' @param return.mut Logical. If TRUE, only the mutation plot is returned.
#' @param homology.field Character. Name of the field in the data that contains homology information.
#' @param del.grouping Optional. Grouping parameter used to split deletion events based on their size.
#'
#' @return A list containing the generated ggplot object(s) and the processed data.table used for plotting.
#'
#' @examples
#' \dontrun{
#'   # Assuming 'my_data' is a data.table with the required fields and 'gts' is a vector of genotype names:
#'   result <- tornado.graph.fn(my_data, plot.ins = TRUE, plot.mut = FALSE, x.spacing = 20,
#'                              adjust.size = 69, gts = c("genotype1", "genotype2"))
#'   # To extract the plot:
#'   plot_obj <- result[[1]]
#'   # To extract the processed data:
#'   processed_data <- result[[2]]
#' }
tornado.graph.fn = function(data,
                            plot.ins = T,
                            plot.mut = F,
                            x.spacing = 20,
                            adjust.size = 69,
                            gts = NA_character_,
                            all = F,
                            hom.min = NA_integer_,
                            return.ins = F,
                            return.del = F,
                            return.mut = F,
                            homology.field = "homology",
                            del.grouping = NULL){
  if(adjust.size == 69) {
    max.del.size = 100
    x.max = 110
  }
  else {
    max.del.size = 201
    x.max = 250}
  group_by1 <- function(x, cols = NULL) group_by(x, across(all_of(c("cell_line", "genotype", cols))))
  group_by2 <- function(x, hom = F) {
    cols <-c("genotype", "first_del", "size", "total_dels", "total_ins", "total_indels",
             "total_muts", "total_modified", "total_reads", "ins", "del", "mut")
    if (hom == T ){ cols <- c(cols, "homology", "homeology")}
    group_by(x, across(all_of(cols)))
  }
  if(!"deletion_density" %in% names(data)){
    data.gt = data[genotype %in% gts][del == T]
    data.gt$genotype <- factor(data.gt$genotype, levels = gts)
    data.gt$homology <- data.gt %>% dd(homology.field)
    data.gt <- data.gt %>%
      group_by1 %>%
      summarize(
        first_del,
        size = deletion,
        reads, ins,del,mut,sample,
        homology, homeology,
        across(starts_with("total_"), ~ sum(unique(.)), .names = "{.col}")
      )%>%
      ungroup() %>%
      group_by2(hom = T) %>%
      summarize(reads = sum(reads),
                samples = paste0(unique(sample), collapse = ", ")) %>%
      as.data.table()
    data.gt[, class := "del"][, fill.color := "#E0115F"]
    data.ins <- data[genotype %in% gts][ins == T] %>% #[mut == F][del == F] %>%
      group_by1 %>%
      summarize(
        first_del = adjust.size - insertion_sizes / 2,
        size = insertion_sizes,
        reads, ins,del,mut,sample,
        across(starts_with("total_"), ~ sum(unique(.)), .names = "{.col}")
      ) %>%
      ungroup() %>%
      group_by2(hom = F) %>%
      summarize(reads = sum(reads),
                samples = paste0(unique(sample), collapse = ", ")) %>%
      as.data.table()
    data.ins[size == 1, homology := -2]
    data.ins[size > 1,  homology := -1]
    data.ins[, class := "ins"][, fill.color := "black"]
    data.mut <- data[genotype %in% gts][mut == T][ins == F][del == F] %>%
      rowwise() %>%
      mutate(start_pos = substitution_positions[[1]],
             end_pos = rev(substitution_positions)[[1]]) %>%
      as.data.table() %>%
      group_by1 %>%
      summarize(
        first_del = start_pos,
        size = (end_pos + 1 - start_pos),
        reads, ins,del,mut,sample,
        across(starts_with("total_"), ~ sum(unique(.)), .names = "{.col}")
      ) %>%
      ungroup() %>%
      group_by2(hom = F) %>%
      summarize(reads = sum(reads),
                samples = paste0(unique(sample), collapse = ", ")) %>%
      as.data.table()
    data.mut[, homology := -3]
    data.mut[, class := "mut"][, fill.color := "red"]
    data.all <- data.gt
    if(plot.ins == T)
      data.all <- data.ins %>% rbind(data.all, fill = T)
    if(plot.mut == T) {
      data.all <- data.mut %>% rbind(data.all, fill = T)}
    data.all[, deletion_density := reads / case_when(
                                             plot.ins == T & plot.mut == T ~  total_modified,
                                             plot.ins == T ~ total_indels,
                                             plot.mut == T ~ total_dels + total_muts,
                                             T ~ total_dels)]
    #browser()
    data.all %>% arrange(desc(class), size, desc(homology), first_del) %>%
      group_by(genotype) %>%
      mutate(
        ymin = cumsum(c(0, head(deletion_density, -1))),
        ymax = cumsum(deletion_density),
        grouping = factor(case_when(
          size <= 4 & size >= 1 & class == "del" ~ "1-4bp",
          size <= 25 & size >= 5 & class == "del" ~ "5-25bp",
          size > 25 & class == "del" ~ ">25bp",
          T ~ NA
        ), levels = c("1-4bp", "5-25bp", ">25bp")),
        class = factor(class, levels = c("mut","ins", "del"))) %>%
      arrange(desc(class), grouping) %>%
      as.data.table()-> data.all
    y.min <- 0; y.max <- 1
  } else {
    #browser()
    data.all <- data
    if(return.ins) {
      data.all <- data.all[class == "ins"]
    } else if(return.mut) {
      data.all <- data.all[class == "mut"]
    } else if (return.del & is.null(del.grouping)) {
      data.all <- data.all[class == "del"]
      unique.del.groups <- as.numeric(unique(sort(data.all$grouping)))
      del.tranche.groupings = list(unique.del.groups[length(unique.del.groups)])
      hom.min.groupings = list(NA_integer_, 5)
      del.tranche.groupings = rep(del.tranche.groupings, length(hom.min.groupings))
      mapply(function(del.grouping, hom.min.grouping){
        tornado.graph.fn(data, plot.ins = plot.ins,
                         plot.mut = plot.mut,
                         x.spacing = x.spacing, adjust.size = adjust.size,
                         gts = gts, return.del = T, del.grouping = del.grouping,
                         hom.min = hom.min.grouping) -> group.call
        return(group.call[[1]])
      }, del.tranche.groupings, hom.min.groupings,SIMPLIFY = FALSE) -> all.del.plots
      return(list(do.call(plot_grid, c(all.del.plots,
                                       list(nrow = 1,
                                            rel_widths = c(3,1,3,3.8)))), data))
    } else if (return.del & !is.null(del.grouping)) {
      data.all <- data.all[class == "del"][as.numeric(grouping) %in% del.grouping]
      if(!is.na(hom.min)) {
        #browser()
        data.all <- data.all[homology >= hom.min]
        data.all %>%
         group_by(genotype) %>%
         arrange(desc(size)) %>%
          mutate(
            ymax = 1 - cumsum(c(0, head(deletion_density, -1))),
            ymin = 1 - cumsum(deletion_density),
          ) %>%
          arrange(desc(class), grouping) %>%
          as.data.table()-> data.all
      }
    }
    y.min <- floor(min(data.all$ymin, na.rm = T) * 1000) / 1000
    y.max <- ceiling(max(data.all$ymax, na.rm = T) * 1000) / 1000
  }
  fill.values <- c(rep("white", 1), col_fun_prop(1:max(data.all$homology, na.rm = T)))
  names(fill.values) <- as.factor(c(0:max(data.all$homology, na.rm = T)))
  fill.labels <- paste0(as.character(0:max(data.all$homology, na.rm = T)), "bp del")
  lapply(gts, function(gt){
    data.all.gt <- data.all[genotype == gt]
    print(gt)
    x.min <- -x.max
    n.breaks.mut <- 5
    if(return.mut){n.breaks.mut <- 3}
    if(return.ins){n.breaks.mut <- 3}
    if(return.del){n.breaks.mut <- 4}
    if(!is.na(hom.min)) {
      n.breaks.mut <- 2
    }
    if ((plot.ins | return.ins) & return.mut == F & return.del == F) {
      fill.values = c("#FF8C00","#FFDAB9", fill.values);
      fill.labels = c("1bp ins", ">1bp ins", fill.labels)
      names(fill.values) <- as.factor(as.numeric(c(-2, -1, names(fill.values)[!names(fill.values) == ""])))
    }
    if ((plot.mut | return.mut) & return.ins == F & return.del == F) {
      fill.values = c("#b04040", fill.values); fill.labels = c("snv", fill.labels)
      names(fill.values) <- as.factor(as.numeric(c(-3, names(fill.values)[!names(fill.values) == ""])))
    }
    if(return.mut | return.ins | return.del) {
      variant.classes <- data.all.gt %>% as.data.table %>%
        arrange(class, grouping) %>%
        group_by(class, grouping) %>%
        summarize(reads = sum(reads), deletion_density = sum(deletion_density),
                  ymin = min(ymin, na.rm = T),
                  ymax = max(ymax, na.rm = T)) %>%
        as.data.table()
    } else{
      #browser()
      variant.classes <- data.all.gt %>% as.data.table %>%
        group_by(class, grouping) %>%
        summarize(reads = sum(reads), deletion_density = sum(deletion_density)) %>%
      ungroup() %>%
        arrange(class, grouping) %>%
        mutate(
          ymin = cumsum(c(0, head(deletion_density, -1))),
          ymax = cumsum(deletion_density)
        ) %>%
        as.data.table()
    }
    if(all == F) {
      theme2 <- theme(axis.text.x = element_text(hjust = 0.6, size = 3, margin = margin(t = 1.2), angle = 0),
            plot.title = element_text(hjust = 0.5, size = 6, margin = margin(b = 0.4), face = "bold", color = "white"),
            plot.subtitle = element_text(hjust = 0.5, size = 4, margin = margin(t = 0, b = 0.2)),
            plot.caption = element_text(hjust = 0.5, size = 6, vjust = 1),
            axis.title.y = element_blank(),#text(size = 4, margin = margin(b = 2, r = 1.75)),
            axis.text.y = element_text(size = 3, vjust = 0.5, margin = margin(r = 0.4)),
            axis.line = element_line(size = 0.1),
            axis.ticks = element_line(size = 0.1),     # Major ticks
            axis.ticks.length = unit(0.02, "cm"),
            axis.title.x = element_blank(),
            legend.position = "none",
            plot.margin = unit(c(0.2, 1.2, 0.2, 0.2), "mm")) 
    } else{
      theme2 <- theme(legend.position = "left")
    }
    #### DISTANCES ###
    x.axis.distance = 2 * x.max
    pad = x.axis.distance / 40
    size.marker = 10 * pad
    tranche.marker = 5 * pad
    ggplot(data.all.gt) +
      geom_rect(aes(
        fill = as.factor(homology),
        ymin = x.min - pad,
        ymax = x.max,
        xmin = ymin,
        xmax = ymax
      )) +
      scale_y_continuous(name = "Sequence Position (bp)",
                         limits = c(x.min - pad, x.max + tranche.marker), #size.marker),
                         breaks = seq(floor(x.min - pad),
                                      x.max + pad)[seq(floor(x.min - pad), x.max + pad) %% x.spacing == 0],
                         labels = seq(floor(x.min - pad),
                                      x.max + pad)[seq(floor(x.min - pad), x.max + pad) %% x.spacing == 0],
                         expand = c(0, 0)) +
      scale_x_continuous(name = "Indel Density",
                         limits = c(y.min, y.max),
                         n.breaks = n.breaks.mut,
                         labels = function(x) sprintf("%.2f", x),
                         expand = c(0, 0)) +
      scale_fill_manual(values = fill.values,
                        labels = fill.labels,
                        name = paste0(homology.field, " / event")) +
      new_scale_fill() +
      geom_rect(aes(
        ymin = first_del - adjust.size,
        ymax = first_del + as.numeric(size) - adjust.size,
        xmin = ymin,
        xmax = ymax#,
      ), fill = "#808080") +
      new_scale_fill() +
      geom_rect(aes(
        fill = grouping,
        ymin = x.max, #+ tranche.marker,
        ymax = x.max + tranche.marker,
        xmin = ymin,
        xmax = ymax
      ), data = variant.classes, color = "black", size = 0.08) +
      scale_fill_manual(values = col.groupings,
                        na.value = "black",
                        name = "Deletion size") +
      theme_pubr() +
      labs(title =  gt,
           subtitle = gsub("e+0","",paste0("n=", scientific(sum(data.all.gt$reads))))
           ) +
      theme(axis.text.x = element_text(hjust = 0.6, size = 3, margin = margin(t = 1.5)),
            plot.title = element_text(hjust = 0.5, size = 6, margin = margin(b = 0.4), face = "bold"),
            plot.subtitle = element_text(hjust = 0.5, size = 4, margin = margin(t = 0, b = 0.2)),
            plot.caption = element_text(hjust = 0.5, size = 6, vjust = 1),
            axis.title.y = element_text(size = 4, margin = margin(b = 2, r = 1.75)),
            axis.text.y = element_text(size = 3.25, vjust = 0.5, margin = margin(b = 1,
                                                                                 r = 0.7)),
            axis.line = element_line(size = 0.1),
            axis.ticks = element_line(size = 0.2),     # Major ticks
            axis.ticks.length = unit(0.04, "cm"),
            axis.title.x = element_blank(),
            legend.key.size = unit(0.15, 'cm'), #change legend key size
            legend.key.height = unit(0.15, 'cm'), #change legend key height
            legend.key.width = unit(0.15, 'cm'), #change legend key width
            legend.spacing.y = unit(0.02, 'cm'),
            legend.title = element_text(size=3, margin = margin(b = 2)),
            legend.text = element_text(size=3, margin = margin(l = 0.4)),
            plot.margin = unit(c(0.2, 1.1, 0.2, 0.2), "mm"),
            legend.position = "right") +
      theme2-> p
    return(p)
  }) -> all.plots
  #browser()
  get_legend(all.plots[[1]]) -> legend
  all.plots <- lapply(1:length(all.plots),
                      function(x) return(all.plots[[x]] + theme(legend.position = "none")))
  if(!(return.ins | return.del | return.mut)){
    do.call(plot_grid, c(all.plots, ncol = 1))-> plots
    plot_grid(legend,
              plots, nrow = 1, rel_widths = c(1, 4)) -> plots
  } else{
    all.plots <- lapply(1:length(all.plots),
                        function(x) return(all.plots[[x]]))
    do.call(plot_grid,
            c(all.plots,
              ncol = 1, align = "h"))-> plots
  }
  return(list(plots, data.all))
}

#' Generate Combined Tornado Plot for a Specific Cell Line
#'
#' This function creates a combined tornado plot for a given cell line by processing variant data.
#' It filters the dataset based on the cell line and optionally specific genotypes, generates intermediate
#' plots for deletions, insertions, and mutations by calling tornado.graph.fn, and arranges them into an overall
#' plot. The final plot is saved to a PDF file.
#'
#' @param data A data.table or data.frame containing variant data.
#' @param cell A character string indicating the cell line identifier used to filter the data.
#' @param gts Optional. A character vector of genotype names to be used. If NULL, genotypes are inferred from the data.
#' @param suffix Character. Indicates plot size parameters; e.g., "short" for smaller plots and "long" (or any alternative)
#'               for larger plots. Determines output file naming and certain plotting aesthetics.
#' @param plot.ins Logical. If TRUE, insertion variants are included in the plot.
#' @param plot.mut Logical. If TRUE, mutation variants are included in the plot.
#' @param field Character. The field name (column) to be used for homology information, such as "homology" or "homeology".
#' @param run.twice Logical. If TRUE, the function calls itself a second time with an alternate field (e.g., "homeology")
#'                to generate a second set of plots.
#' @param rel.widths Numeric vector. Specifies the relative widths for arranging the composite plot layout.
#' @param output_dir Character. Directory path where the output PDF file will be saved.
#'
#' @return This function does not return a value. It saves the generated plot as a PDF file in the specified output directory.
#'
#' @examples
#' \dontrun{
#'   # With a dataset 'variant_data' and cell line "CellA":
#'   tornado_plot(variant_data, cell = "CellA", gts = c("WT", "Mutant"),
#'                suffix = "short", plot.ins = TRUE, plot.mut = FALSE,
#'                field = "homology", run.twice = TRUE,
#'                rel.widths = c(1.8, 0.4, 0.4, 1.3),
#'                output_dir = "~/public_html/plots/homology/tornado/")
#' }
#' @export
tornado_plot = function(data, cell, gts = NULL,
                        suffix = "short",
                        plot.ins = T, plot.mut = F,
                        field = "homology",
                        run.twice = T,
                        rel.widths = c(1.8,0.4,0.4,1.3),
                        output_dir = "~/public_html/plots/homology/tornado/"){
  if(suffix == "short"){
    suffix2 <- "245bp"; x.spacing <- 25; adjust.size = 69; max.del.size = 100
  }else{
    suffix2 <- "503bp"; x.spacing <- 40; adjust.size = 251; max.del.size = 201
  }
  if(is.null(gts))
    gts <- data[(cell_line == cell)]$genotype %>% unique; data <- data[(cell_line == cell)]
  ins.plots <- NULL
  tornado.graph.fn(data,
                   plot.ins = plot.ins,
                   plot.mut = plot.mut,
                   x.spacing, adjust.size,
                   gts, all = T,
                   homology.field = field) -> all.plots
  tornado.graph.fn(all.plots[[2]],
                   plot.ins = plot.ins,
                   plot.mut = plot.mut,
                   x.spacing, adjust.size,
                   gts, homology.field = field,
                   return.del = T) -> del.plots
  ins.plots <- NULL
  if(plot.ins)
    tornado.graph.fn(all.plots[[2]],
                     plot.ins = plot.ins,
                     plot.mut = plot.mut,
                     x.spacing, adjust.size,
                     gts, homology.field = field,
                     return.ins = plot.ins) -> ins.plots
  mut.plots <- NULL
  if(plot.mut)
    tornado.graph.fn(all.plots[[2]],
                     plot.ins = plot.ins,
                     plot.mut = plot.mut,
                     x.spacing, adjust.size,
                     gts, homology.field = field,
                     return.mut = plot.mut) -> mut.plots
  plots.coerced = list(all.plots[[1]],
                       mut.plots[[1]],
                       ins.plots[[1]],
                       del.plots[[1]])
  do.call(plot_grid,
          c(plots.coerced[!sapply(plots.coerced, is.null)],
            list(rel_widths = rel.widths[!sapply(plots.coerced, is.null)],
                 nrow = 1, align = "v")))-> plots
  suffix.all <- paste0(gts, collapse = "-")
  if(plot.ins == T)
    suffix.all <- paste0("INS_", suffix.all)
  if(plot.mut == T)
    suffix.all <- paste0("MUT_", suffix.all)
  dir.to.create <- paste0(output_dir, cell, "_", suffix, "_", field)
  dir.create(dir.to.create)
  ggsave2(
    paste0(output_dir,
           cell, "_", suffix, "_", field, "/",
           cell, "_", suffix, "_", suffix.all, ".pdf" ),
    plot = plots,
    width = 5.2,
    height = 0.75 * length(gts),
    dpi = 600,
    units = "in"
  )
  if(run.twice)
    tornado_plot(data, cell, gts,
                 suffix = suffix,
                 plot.ins, plot.mut,
                 field = "homeology",
                 run.twice = F)
}

##############################################################################

#' jrafailov Monday, April 28, 2025 11:47:24
#' process all datasets as follows:

genotype_sets <- list(
  list(data = short.data, cell = "DLD1", gts = c("WT", "PolQ-PolDel", "PolQ-HelDel", "PolQ-PolDHelD")),
  list(data = short.data, cell = "U2OS", gts = c("DMSO", "RP6685-10uM", "RP9913-10uM", "ART558-10uM", "ART558-25uM", "NVB-50uM", "NVB-100uM")),
  list(data = short.data, cell = "Hek293T", gts = c("siNT", "siBRCA2", "siRAD52", "siPolQ")),
  list(data = short.data, cell = "Hek293T", gts = c("siNT", "siRAD52", "siPolQ", "siBRCA2", "siBRCA2_PolQ", "siBRCA2_R52")),
  list(data = short.data, cell = "Hek293T", gts = c("siNT", "siBRCA2", "siPolH", "siBRCA2_PolH", "siPolL", "siBRCA2_PolL", "siPolQ", "siBRCA2_PolQ")),
  list(data = short.data, cell = "Hek293T", gts = c("siNT", "siBRCA2", "siRAD52", "siBRCA2_R52", "siPolQ", "siBRCA2_PolQ", "siPolQ_R52"))
)

mclapply(genotype_sets, function(set) {
  mcmapply(function(plot.ins, plot.mut) {
``    tornado_plot(set$data, set$cell, gts = set$gts, plot.ins = plot.ins, plot.mut = plot.mut)
  }, plot.ins = c(T, F, T), plot.mut = c(F, F, T), mc.cores = 3)
}, mc.cores = 6)

genotype_sets <- list(
  list(data = long.data, cell = "DLD1", gts = c("WT", "PolQ-PolDel", "PolQ-HelDel", "PolQ-PolDHelD")),
  list(data = long.data, cell = "U2OS", gts = c("DMSO", "RP6685-10uM", "RP9913-10uM", "ART558-25uM", "NVB-100uM")),
  list(data = long.data, cell = "Hek293T", gts = c("siNT", "siBRCA2", "siBRCA2_RAD52", "siBRCA2_PolQ")),
  list(data = long.data, cell = "Hek293T", gts = c("siNT", "siBRCA2", "siBRCA2_RAD52", "siBRCA2_PolQ", "siBRCA2_PolH", "siBRCA2_PolL"))
)

mclapply(genotype_sets, function(set) {
  mcmapply(function(plot.ins, plot.mut) {
    tornado_plot(set$data, set$cell, gts = set$gts, plot.ins = plot.ins, plot.mut = plot.mut, suffix = "long")
  }, plot.ins = c(T, F, T), plot.mut = c(F, F, T), mc.cores = 3)
}, mc.cores = 4)



