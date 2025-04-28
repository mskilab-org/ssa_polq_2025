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

#' jrafailov Tuesday, Oct 15, 2024 10:37:03 AM
#' hek293T 245bp samples pull in from scratch

q_window = "cut_site"

project <- short[2]

projects = c("Project_13108", # HEK293T 245bp
"PE250_Hek", # HEK293T 245bp
"PE150_extra", # HEK293T 245bp
"PE150_DLD1", # DLD1 245bp
"Project_14818", # DLD1 inhibitor 245bp
"Project_14789", #Hek293T 503bp
"Project_14783", #mixed Hek293T/DLD1 503bp
"Project_14629", #DLD1 503bp
"Project_15074") #DLD1 inhibitor 503bp)

project <- "Project_13108"

short <- c("Project_13108", # HEK293T 245bp
           "PE250_Hek", # HEK293T 245bp
           "PE150_extra", # HEK293T 245bp
           "PE150_DLD1", # DLD1 245bp
           "Project_14818")

long <- c("Project_14789", #Hek293T 503bp
          "Project_14783", #mixed Hek293T/DLD1 503bp
          "Project_14629", #DLD1 503bp
          "Project_15074")

genewiz <- c("PE250_Hek", # HEK293T 245bp GENEWIZ!!!!!
             "PE150_extra") # HEK293T 245bp GENEWIZ !!!)
igo <- c("Project_13108", # HEK293T 245bp IGO
         "PE150_DLD1", # DLD1 245bp IGO
         "Project_14818",
         "Project_14789", #Hek293T 503bp
         "Project_14783", #mixed Hek293T/DLD1 503bp
         "Project_14629", #DLD1 503bp
         "Project_15074")

mclapply(long., function(project){
  print(q_window)
  print(project)
  if(q_window == "cut_site"){
    dt <- readRDS(paste0("/gpfs/commons/home/sbrylka/Projects/homeology/",project,"/dt.rds"))
  }
  if(q_window == "whole_amplicon"){
    dt <- readRDS(paste0("~/projects/polq/",project,"/dt2.rds"))
  }
  if(project %in% short) {pos <- 69; end_pos = 245}
  if(project %in% long) {pos <- 251; end_pos = 503}
  if(project %in% genewiz)
    grp <- "GENEWIZ"
  if(project %in% igo)
    grp <- "IGO"
  stack_plot_data <- mclapply(1:nrow(dt),
                              function(x){
                                path = paste0(dt$crisp.out[x],"/Alleles_frequency_table.txt")
                                print(file.exists(path))
                                a <- read.delim(path) %>% as.data.table()
                                total_reads = as.numeric(sum(a$X.Reads))
                                w <- a %>%
                                  filter(n_deleted > 0 | n_inserted > 0 | n_mutated > 0) %>%
                                  as.data.table()
                                w$deletion_coordinates <- gsub("\\[|\\]|\\(|\\)|\\s", "",
                                                               w$deletion_coordinates) %>%
                                  str_split(",") %>% lapply(as.numeric)
                                w$insertion_sizes <- gsub("\\[|\\]", "", w$insertion_sizes) %>% lapply(as.numeric) %>% unlist
                                w$insertion_coordinates <- gsub("\\[|\\]|\\(|\\)|\\s", "", w$insertion_coordinates) %>% str_split(",") %>% lapply(as.numeric)
                                w$substitution_positions <- gsub("\\[|\\]|\\s", "", w$substitution_positions) %>% str_split(",") %>% lapply(as.numeric)
                                w$substitution_values <-   gsub("\\[|\\]", "", w$substitution_values) %>% str_split(" ") %>% lapply(function(x)return( gsub("\\'", "", x)))
                                w$first_del <- sapply(w$deletion_coordinates, function(x)return(x[1]))
                                w$last_del  <- sapply(w$deletion_coordinates, function(x)return(x[2]))
                                w$del <- !(is.na(w$first_del) & is.na(w$last_del))
                                w$ins <- !is.na(w$insertion_sizes)
                                w$mut <- !is.na(w$substitution_positions)
                                w <- w[!first_del == 0 | is.na(first_del)]
                                w <- w[!last_del == (end_pos - 1) | is.na(last_del)]
                                w <- w[first_del <= pos & last_del >= pos | is.na(first_del) | is.na(last_del)]
                                w[, deletion_size := ifelse(!is.na(first_del), last_del - first_del, NA)]
                                #total_dels = as.numeric(sum(w[!is.na(deletion_coordinates)]$X.Reads))
                                total_modified = as.numeric(sum(w$X.Reads))
                                total_ins = as.numeric(sum(w[ins == T]$X.Reads))
                                total_dels = as.numeric(sum(w[del == T]$X.Reads))
                                total_muts = as.numeric(sum(w[mut == T & ins == F & del == F]$X.Reads))
                                  d <- w %>%
                                    mutate(sample = dt$sample[x], reads = X.Reads, deletion = n_deleted) %>%
                                    select(sample, reads, deletion_coordinates, deletion,
                                           insertion_sizes, insertion_coordinates,substitution_positions,
                                           substitution_positions, first_del, last_del, del, ins, mut) %>%
                                    mutate(total_dels = total_dels,
                                           total_ins = total_ins,
                                           total_indels = total_dels + total_ins,
                                           total_modified = total_modified,
                                           total_muts = total_muts,
                                           total_reads = total_reads) %>%
                                    as.data.table()
                                return(d)
                              }, mc.cores = 16)
  stack_plot_data <- do.call(rbind, stack_plot_data)
  if(project == "Project_13108"){
    w <- stack_plot_data %>%
    rowwise() %>%
    mutate(parsed = str_split(sample,pattern = "_")) %>%
    mutate(experiment = parsed[1],
           cell_line = parsed[2],
           genotype = parsed[3]) %>%
  select(-parsed) %>%
    data.table()
    w$group <- grp
    summary = w[, .(sample, experiment, cell_line, genotype, total_reads, total_modified, total_indels, total_dels, total_muts, total_ins)] %>%
    unique %>%
    group_by(experiment, genotype, cell_line) %>%
    summarize(
      n = n(),
      samples = paste0(unique(sample), collapse = ", "),
      total_reads = sum(total_reads),
      total_modified = sum(total_modified),
      total_indels = sum(total_indels),
      total_dels = sum(total_dels),
      total_ins = sum(total_ins),
      total_muts = sum(total_muts)
    ) %>%
    as.data.table()
    summary$group <- grp
  saveRDS(w, paste0("~/projects/polq/jr_stash/all_data2/", project, "_data.rds"))
  saveRDS(summary, paste0("~/projects/polq/jr_stash/all_data2/", project, "_summary.rds"))
	}
if(project == "PE250_Hek" | project == "PE150_extra"){
  w <- stack_plot_data %>%
    rowwise() %>%
    mutate(parsed = str_split(sample,pattern = "-")) %>%
    mutate(experiment = parsed[1],
           cell_line = parsed[2],
           genotype = paste0(parsed[3],"_", parsed[4])) %>%
    select(-parsed) %>%
    data.table()
  w$group <- grp
  summary = w[, .(sample, experiment, cell_line, genotype, total_reads, total_modified, total_indels, total_dels, total_muts, total_ins)] %>%
    unique %>%
    group_by(experiment, genotype, cell_line) %>%
    summarize(
      n = n(),
      samples = paste0(unique(sample), collapse = ", "),
      total_reads = sum(total_reads),
      total_modified = sum(total_modified),
      total_indels = sum(total_indels),
      total_dels = sum(total_dels),
      total_ins = sum(total_ins),
      total_muts = sum(total_muts)
    ) %>%
    as.data.table()
  summary$group <- grp
  saveRDS(w, paste0("~/projects/polq/jr_stash/all_data2/", project, "_data.rds"))
  saveRDS(summary, paste0("~/projects/polq/jr_stash/all_data2/", project, "_summary.rds"))}
  if(project == "PE150_DLD1"){
  w <- stack_plot_data %>%
    rowwise() %>%
    mutate(parsed = str_split(sample,pattern = "_")) %>%
    mutate(experiment = parsed[1],
           cell_line = "DLD1",
           genotype = parsed[2]) %>%
    mutate(genotype = gsub("DLD1-","",genotype)) %>%
    mutate(genotype = ifelse(genotype == "PolQ-HellDel", "PolQ-HelDel", genotype)) %>%
    select(-parsed) %>%
    data.table()
  w$group <- grp
  summary = w[, .(sample, experiment, cell_line, genotype, total_reads, total_modified, total_indels, total_dels, total_muts, total_ins)] %>%
    unique %>%
    group_by(experiment, genotype, cell_line) %>%
    summarize(
      n = n(),
      samples = paste0(unique(sample), collapse = ", "),
      total_reads = sum(total_reads),
      total_modified = sum(total_modified),
      total_indels = sum(total_indels),
      total_dels = sum(total_dels),
      total_ins = sum(total_ins),
      total_muts = sum(total_muts)
    ) %>%
    as.data.table()
  summary$group <- grp
  saveRDS(w, paste0("~/projects/polq/jr_stash/all_data2/", project, "_data.rds"))
  saveRDS(summary, paste0("~/projects/polq/jr_stash/all_data2/", project, "_summary.rds"))
  }
  if(project == "Project_14818" | project == "Project_15074"){
  w <- stack_plot_data %>%
    rowwise() %>%
    mutate(parsed = str_split(sample,pattern = "_")) %>%
    mutate(experiment = parsed[1],
           cell_line = "U2OS",
           treatment = parsed[3],
           volume = parsed[4]) %>%
    mutate(volume = ifelse(volume == "C", NA,volume)) %>%
    select(-parsed) %>%
    data.table()
  w$group <- grp
  summary = w[, .(sample, experiment, cell_line, treatment, volume, total_reads, total_modified, total_indels, total_dels, total_muts, total_ins)] %>%
    unique %>%
    group_by(experiment, cell_line, treatment, volume) %>%
    summarize(
      n = n(),
      samples = paste0(unique(sample), collapse = ", "),
      total_reads = sum(total_reads),
      total_modified = sum(total_modified),
      total_indels = sum(total_indels),
      total_dels = sum(total_dels),
      total_ins = sum(total_ins),
      total_muts = sum(total_muts)
    ) %>%
    as.data.table()
  summary$group <- grp
  saveRDS(w, paste0("~/projects/polq/jr_stash/all_data2/", project, "_data.rds"))
  saveRDS(summary, paste0("~/projects/polq/jr_stash/all_data2/", project, "_summary.rds"))
  }
  if(project == "Project_14783" | project == "Project_14789" | project == "Project_14629"){
  w <- stack_plot_data %>%
    rowwise() %>%
    mutate(parsed = str_split(sample,pattern = "_")) %>%
    mutate(experiment = parsed[2],
           cell_line = parsed[3],
           genotype = paste0(parsed[4],"_",parsed[5])) %>%
    mutate(c = cell_line) %>%
    mutate(cell_line = ifelse(genotype == "C_IGO", "DLD1", cell_line),
           genotype = ifelse(genotype == "C_IGO", gsub("DLD1-","", c), genotype)) %>%
    select(-c(parsed,c)) %>%
    data.table()
  summary = w[, .(sample, experiment, cell_line, genotype, total_reads, total_modified, total_indels, total_dels, total_muts, total_ins)] %>%
    unique %>%
    group_by(experiment, genotype, cell_line) %>%
    summarize(
      n = n(),
      samples = paste0(unique(sample), collapse = ", "),
      total_reads = sum(total_reads),
      total_modified = sum(total_modified),
      total_indels = sum(total_indels),
      total_dels = sum(total_dels),
      total_ins = sum(total_ins),
      total_muts = sum(total_muts)
    ) %>%
    as.data.table()
  summary$group <- grp
  saveRDS(w, paste0("~/projects/polq/jr_stash/all_data2/", project, "_data.rds"))
  saveRDS(summary, paste0("~/projects/polq/jr_stash/all_data2/", project, "_summary.rds"))}
}, mc.cores = 5)



#' jrafailov Wednesday, Oct 16, 2024 04:24:24 PM
#' QC metrics

project = c("Project_13108", # HEK293T 245bp
            "PE250_Hek", # HEK293T 245bp
            "PE150_extra", # HEK293T 245bp
            "PE150_DLD1", # DLD1 245bp
            "Project_14818", # DLD1 inhibitor 245bp
            "Project_14783", #mixed Hek293T/DLD1 503bp
            "Project_14789", #Hek293T 503bp
            "Project_14629", #DLD1 503bp
            "Project_15074") #DLD1 inhibitor 503bp)

short <- c("Project_13108", # HEK293T 245bp
           "PE250_Hek", # HEK293T 245bp
           "PE150_extra", # HEK293T 245bp
           "PE150_DLD1", # DLD1 245bp
           "Project_14818")

long <- c("Project_14789", #Hek293T 503bp
          "Project_14783", #mixed Hek293T/DLD1 503bp
          "Project_14629", #DLD1 503bp
          "Project_15074")



lapply(short, function(x){
  readRDS(paste0(paste0("~/projects/polq/jr_stash/all_data/", x, "_summary.rds"))) -> a
  a$project <- x
  return(a)
}) %>% rbindlist(fill = T) %>%
  mutate(
    genotype = gsub("_C", "", genotype)
  ) %>%
  mutate(
    genotype = case_when(
      is.na(genotype) ~ paste0(treatment, ifelse(is.na(volume), "", paste0("-", volume))),
   T ~ genotype
    )
  )-> short.summary

lapply(short, function(x){
  readRDS(paste0(paste0("~/projects/polq/jr_stash/all_data2/", x, "_summary.rds"))) -> a
  a$project <- x
  return(a)
}) %>% rbindlist(fill = T) %>%
  mutate(
    genotype = gsub("_C", "", genotype)
  ) %>%
  mutate(
    genotype = case_when(
      is.na(genotype) ~ paste0(treatment, ifelse(is.na(volume), "", paste0("-", volume))),
      T ~ genotype
    )
  )-> short.summary2


short.n.summary = short.summary2 %>%
  transmute(genotype, experiment, cell_line, sequencing_n = n, samples = paste0("[", samples, "]")) %>%
  group_by(cell_line, genotype) %>%
  summarize(
    biological_n = n(),
    technical_n = sum(sequencing_n),
    experiment = paste0(unique(experiment), collapse = ", "),
    samples = paste0(unique(samples), collapse = ", " )
  ) %>%
  unique %>%
  as.data.table()

readr::write_tsv(short.n.summary, "~/Projects/tmp/245bp.samplesummary.tsv")

short.n = short.summary[, .(experiment, genotype, cell_line)] %>%
  group_by(genotype, cell_line) %>%
  summarize(n = n(),
            n.exp = length(unique(experiment)))

short.qc.reads <- short.summary2 %>%
  group_by(genotype, cell_line, experiment, group) %>%
  summarize(
    total_reads = sum(total_reads),
    total_modified = sum(total_modified),
    total_indels = sum(total_indels),
    total_ins = sum(total_ins),
    total_dels = sum(total_dels),
    total_muts = sum(total_muts)
  ) %>%
  ungroup() %>%
  pivot_longer(
    starts_with("total"),
    names_prefix = "total_",
    values_to = "n",
    names_to = "field"
  ) %>%
  as.data.table()
short.qc.reads$field <- factor(short.qc.reads$field,
                               levels = c("dels", "ins", "muts", "indels", "modified", "reads"))

short.qc.frac = short.summary2 %>%
  group_by(genotype, cell_line, experiment, group) %>%
  summarize(
    total_reads = sum(total_reads),
    total_modified = sum(total_modified),
    total_indels = sum(total_indels),
    total_ins = sum(total_ins),
    total_dels = sum(total_dels),
    total_muts = sum(total_muts)
  ) %>%
  ungroup() %>%  mutate(
    perc_modified = total_modified / total_reads,
    perc_indels = total_indels / total_reads,
    perc_dels = total_dels / total_reads,
    perc_ins = total_ins / total_reads,
    perc_muts = total_muts / total_reads
    ) %>%
  select(-starts_with("total")) %>%
  pivot_longer(
    starts_with("perc"),
    names_prefix = "perc_",
    values_to = "n",
    names_to = "field"
  ) %>%
  as.data.table()
short.qc.frac$field <- factor(short.qc.frac$field,
                               levels = c("dels", "ins", "muts", "indels", "modified", "reads"))


lapply(unique(short.qc.reads$cell_line), function(x){
  ggplot(short.qc.reads[cell_line == x] ,
         aes(x = experiment, y = n)) +
    geom_bar(aes(fill = field), stat = "identity" , color = "black") +
    facet_grid(field ~ genotype + group,
#               ncol = length(unique(short.qc.reads[cell_line == x]$genotype)),
               scales = "free",
               switch = "x",
               space = "free_x") +
  scale_fill_scico_d(palette = "batlow") +
#  scale_color_scico_d(palette = "batlow") +
    #scale_color_manual(values = c(rep(color0, 2), rep("black", 5))) +
    theme_bw() +
  theme(
    panel.grid.major.x = element_blank(),
    axis.text.y = element_text(size = 6),
    axis.text.x = element_text(size = 5, angle = 90, vjust = 0.5, hjust = 0.5),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_line(size = .1, color = "black"),
    panel.grid.minor.y = element_blank(),
    strip.background = element_blank(),
    axis.ticks.x = element_blank(),
    strip.text.x = element_text(size = 6, angle = 30, vjust = 0.5, hjust = 0.5),
    strip.text.y = element_text(size = 10),
    strip.clip = "off",
    legend.position = "none",
    plot.margin = unit(c(0, 0, 0, 0), "cm")
    ) +
  ylab("reads") +
  labs(title = "raw number of reads", subtitle = x) +
  scale_y_continuous(labels = label_scientific())+
    scale_x_discrete(drop = T) +
  xlab("")-> all
  length.n = length(unique(short.qc.reads[cell_line == x]$experiment)) * length(unique(short.qc.reads[cell_line == x]$genotype))
  ppdf(plot(all), paste0("~/public_html/plots/homology/QC_new/245bp-",
                       x,  "-READS.pdf"),width = length.n * 0.137255, height = 5.5)
})

lapply(unique(short.qc.frac$cell_line), function(x){
  ggplot(short.qc.frac[cell_line == x] ,
         aes(x = experiment, y = n)) +
    geom_bar(aes(fill = field), stat = "identity" , color = "black") +
    facet_grid(field ~ genotype + group,
                                        #               ncol = length(unique(short.qc.reads[cell_line == x]$genotype)),
               switch = "x",
               scales = "free",
               space = "free_x") +
  scale_fill_scico_d(palette = "batlow") +
#  scale_color_scico_d(palette = "batlow") +
    #scale_color_manual(values = c(rep(color0, 2), rep("black", 5))) +
    theme_bw() +
    theme(
      panel.grid.major.x = element_blank(),
      axis.text.y = element_text(size = 6),
      axis.text.x = element_text(size = 5, angle = 90, vjust = 0.5, hjust = 0.5),
      panel.grid.minor.x = element_blank(),
      panel.grid.major.y = element_line(size = .1, color = "black"),
      panel.grid.minor.y = element_blank(),
      strip.background = element_blank(),
      axis.ticks.x = element_blank(),
      strip.text.x = element_text(size = 6, angle = 30, vjust = 0.5, hjust = 0.5),
      strip.text.y = element_text(size = 10),
      strip.clip = "off",
      legend.position = "none",
      plot.margin = unit(c(0, 0, 0, 0), "cm")
    ) +
    ylab("% of total reads") +
    labs(title = "fraction of total reads", subtitle = x) +
    scale_y_continuous(labels = label_scientific())+
    scale_x_discrete(drop = T) +
  xlab("")-> all
  length.n = length(unique(short.qc.reads[cell_line == x]$experiment)) * length(unique(short.qc.reads[cell_line == x]$genotype))
  ppdf(plot(all), paste0("~/public_html/plots/homology/QC_new/245bp-",
                       x,  "-FRAC.pdf"),width = length.n * 0.137255, height = 5.5)
})

lapply(unique(short.qc.frac$cell_line), function(x){
  data = short.qc.frac[cell_line == x]
  data.per.replicate = data %>%
    group_by(field) %>%
    t_test(n ~ genotype) %>%
    add_y_position(fun = "mean_sd", step.increase = 0.05, scales = "free_y")
  ggplot(data, aes(x = genotype, y = n)) +
    geom_bar(stat = "summary", fun = "mean", position = position_dodge(), fill = "white", color = "black") +
    geom_errorbar(stat = "summary",fun.data = "mean_se",position = position_dodge(0.9),width = 0.5) +
    geom_jitter(aes(color = group), height = 0, width = 0.4, size = 0.4) +
    labs(x = "genotype", y = "% of all reads", fill = "platform") +
    theme_minimal() +
    facet_wrap(~ field,
               ncol = 1,
               scales = "free_y") +
    stat_pvalue_manual(data.per.replicate, hide.ns = TRUE) +
    scale_fill_scico(palette = "lipari") +
    theme(
      axis.text.x = element_text(size = 6, angle = 45, vjust = 0.5, hjust = 0.5),
      legend.position = "bottom",
      legend.text = element_text(size = 7),    # Adjust the text size
      legend.title = element_text(size = 7)
    ) +
    ylab("Fraction of Reads") +
    labs(title = "statistical evaluations of proportions", subtitle = x)-> all
  length.n = unique(data$genotype) %>% length
  ppdf(print(all), paste0("~/public_html/plots/homology/QC_new/245bp-",
                         x,  "-PVALUE.pdf"),width = 0.75 * length.n, height = 7)
})


lapply(long, function(x){
  readRDS(paste0(paste0("~/projects/polq/jr_stash/all_data/", x, "_summary.rds"))) -> a
  return(a)
}) %>% rbindlist(fill = T) %>%
  mutate(
    genotype = case_when(
      is.na(genotype) ~ paste0(treatment, ifelse(is.na(volume), "", paste0("-", volume))),
      T ~ genotype
    )
  ) -> long.summary

lapply(long, function(x){
  readRDS(paste0(paste0("~/projects/polq/jr_stash/all_data2/", x, "_summary.rds"))) -> a
  return(a)
}) %>% rbindlist(fill = T) %>%
  mutate(
    genotype = case_when(
      is.na(genotype) ~ paste0(treatment, ifelse(is.na(volume), "", paste0("-", volume))),
      T ~ genotype
    )
  ) -> long.summary2


long.n <- long.summary2[, .(genotype, cell_line)] %>%
  group_by(genotype, cell_line) %>%
  summarize(n = n()) %>%
  arrange(cell_line)

long.n.summary = long.summary2 %>%
  transmute(genotype, experiment, cell_line, sequencing_n = n, samples = paste0("[", samples, "]")) %>%
  group_by(cell_line, genotype) %>%
  summarize(
    biological_n = n(),
    technical_n = sum(sequencing_n),
    experiment = paste0(unique(experiment), collapse = ", "),
    samples = paste0(unique(samples), collapse = ", " )
  ) %>%
  unique %>%
  as.data.table()

readr::write_tsv(long.n.summary, "~/Projects/tmp/503bp.samplesummary.tsv")

long.qc.reads <- long.summary2 %>%
  group_by(genotype, cell_line, experiment, group) %>%
  summarize(
    total_reads = sum(total_reads),
    total_modified = sum(total_modified),
    total_indels = sum(total_indels),
    total_ins = sum(total_ins),
    total_dels = sum(total_dels),
    total_muts = sum(total_muts)
  ) %>%
  ungroup() %>%
  pivot_longer(
    starts_with("total"),
    names_prefix = "total_",
    values_to = "n",
    names_to = "field"
  ) %>%
  as.data.table()
long.qc.reads$field <- factor(long.qc.reads$field,
                               levels = c("dels", "ins", "muts", "indels", "modified", "reads"))

long.qc.frac = long.summary2 %>%
  group_by(genotype, cell_line, experiment, group) %>%
  summarize(
    total_reads = sum(total_reads),
    total_modified = sum(total_modified),
    total_indels = sum(total_indels),
    total_ins = sum(total_ins),
    total_dels = sum(total_dels),
    total_muts = sum(total_muts)
  ) %>%
  ungroup() %>%  mutate(
    perc_modified = total_modified / total_reads,
    perc_indels = total_indels / total_reads,
    perc_dels = total_dels / total_reads,
    perc_ins = total_ins / total_reads,
    perc_muts = total_muts / total_reads
    ) %>%
  select(-starts_with("total")) %>%
  pivot_longer(
    starts_with("perc"),
    names_prefix = "perc_",
    values_to = "n",
    names_to = "field"
  ) %>%
  as.data.table()
long.qc.frac$field <- factor(long.qc.frac$field,
                               levels = c("dels", "ins", "muts", "indels", "modified", "reads"))

lapply(unique(long.qc.reads$cell_line), function(x){
  ggplot(long.qc.reads[cell_line == x] ,
         aes(x = experiment, y = n)) +
    geom_bar(aes(fill = field, color = field), stat = "identity") +
    facet_grid(field ~ genotype + group,
#               ncol = length(unique(short.qc.reads[cell_line == x]$genotype)),
               scales = "free",
               switch = "x",
               space = "free_x") +
  scale_fill_scico_d(palette = "batlow") +
  scale_color_scico_d(palette = "batlow") +
    #scale_color_manual(values = c(rep(color0, 2), rep("black", 5))) +
    theme_bw() +
  theme(
    panel.grid.major.x = element_blank(),
    axis.text.y = element_text(size = 6),
    axis.text.x = element_text(size = 5, angle = 45, vjust = 0.5, hjust = 0.5),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_line(size = .1, color = "black"),
    panel.grid.minor.y = element_blank(),
    strip.background = element_blank(),
    axis.ticks.x = element_blank(),
    strip.text.x = element_text(size = 6, angle = 0, vjust = 0.5, hjust = 0.5),
    strip.text.y = element_text(size = 10),
    strip.clip = "off",
    legend.position = "none",
    plot.margin = unit(c(0, 0, 0, 0), "cm")
    ) +
    ylab("") +
  scale_x_discrete(drop = T) +
  labs(title = "raw number of reads", subtitle = x) +
  xlab("")-> all
  length.n = length(unique(long.qc.reads[cell_line == x]$experiment)) * length(unique(long.qc.reads[cell_line == x]$genotype))
  ppdf(plot(all), paste0("~/public_html/plots/homology/QC_new/503bp-",
                       x,  "-READS.pdf"),width = length.n / 4, height = 6)
})

lapply(unique(long.qc.frac$cell_line), function(x){
  ggplot(long.qc.frac[cell_line == x] ,
         aes(x = experiment, y = n)) +
    geom_bar(aes(fill = field, color = field), stat = "identity") +
    facet_grid(field ~ genotype + group,
               switch = "x",
               scales = "free",
               space = "free_x") +
  scale_fill_scico_d(palette = "batlow") +
  scale_color_scico_d(palette = "batlow") +
    #scale_color_manual(values = c(rep(color0, 2), rep("black", 5))) +
    theme_bw() +
    theme(
      panel.grid.major.x = element_blank(),
      axis.text.y = element_text(size = 6),
      axis.text.x = element_text(size = 5, angle = 45, vjust = 0.5, hjust = 0.5),
      panel.grid.minor.x = element_blank(),
      panel.grid.major.y = element_line(size = .1, color = "black"),
      panel.grid.minor.y = element_blank(),
      strip.background = element_blank(),
      axis.ticks.x = element_blank(),
      strip.text.x = element_text(size = 5, angle = 0, vjust = 0.5, hjust = 0.5),
      strip.text.y = element_text(size = 10),
      strip.clip = "off",
      legend.position = "none",
      plot.margin = unit(c(0, 0, 0, 0), "cm")
    ) +
    ylab("") +
    scale_x_discrete(drop = T) +
    labs(title = "fraction of total reads", subtitle = x) +
  xlab("")-> all
  length.n = length(unique(long.qc.reads[cell_line == x]$experiment)) * length(unique(long.qc.reads[cell_line == x]$genotype))
  ppdf(plot(all), paste0("~/public_html/plots/homology/QC_new/503bp-",
                       x,  "-FRAC.pdf"),width = length.n / 4, height = 6)
})

lapply(unique(long.qc.frac$cell_line), function(x){
  data = long.qc.frac[cell_line == x]
  data.per.replicate = data %>%
    group_by(field) %>%
    t_test(n ~ genotype) %>%
    add_y_position(fun = "mean_sd", step.increase = 0.05, scales = "free_y")
  ggplot(data, aes(x = genotype, y = n)) +
    geom_bar(stat = "summary", fun = "mean", position = position_dodge(), fill = "white", color = "black") +
    geom_errorbar(stat = "summary",fun.data = "mean_se",position = position_dodge(0.9),width = 0.5) +
    geom_jitter(aes(color = group), height = 0, width = 0.4, size = 0.9) +
    labs(x = "genotype", y = "% of all reads", fill = "platform") +
    theme_minimal() +
    facet_wrap(~ field,
               ncol = 1,
               scales = "free_y") +
    stat_pvalue_manual(data.per.replicate, hide.ns = TRUE) +
    scale_fill_scico(palette = "lipari") +
    theme(
      axis.text.x = element_text(size = 8, angle = 45, vjust = 0.5, hjust = 0.5)
    ) +
    ylab("Fraction of Reads") +
    labs(title = "statistical evaluations of proportions", subtitle = x)-> all
  ppdf(print(all), paste0("~/public_html/plots/homology/QC_new/503bp-",
                         x,  "-PVALUE.pdf"),width = 8.5, height = 11)
})


#' jrafailov Thursday, Oct 17, 2024 01:35:20 AM
#' merge the heatmaps with the real homology values and collapse the samples into correct biological replicate #s

project = c("Project_13108", # HEK293T 245bp
            "PE250_Hek", # HEK293T 245bp
            "PE150_extra", # HEK293T 245bp
            "PE150_DLD1", # DLD1 245bp
            "Project_14818", # DLD1 inhibitor 245bp
            "Project_14789", #Hek293T 503bp
            "Project_14783", #mixed Hek293T/DLD1 503bp
            "Project_14789", #Hek293T 503bp
            "Project_14629", #DLD1 503bp
            "Project_15074") #DLD1 inhibitor 503bp)

short <- c("Project_13108", # HEK293T 245bp IGO
           "PE250_Hek", # HEK293T 245bp GENEWIZ!!!!!
           "PE150_extra", # HEK293T 245bp GENEWIZ !!!
           "PE150_DLD1", # DLD1 245bp IGO
           "Project_14818") #U2OS 245bp IGO

long <- c("Project_14789", #Hek293T 503bp
          "Project_14783", #mixed Hek293T/DLD1 503bp
          "Project_14629", #DLD1 503bp
          "Project_15074")

mclapply(short, function(x){
  print(x)
  readRDS(paste0(paste0("~/projects/polq/jr_stash/all_data2/", x, "_data.rds"))) -> a
  a$project <- x
  return(a)
}, mc.cores = 4) %>% rbindlist(fill = T) %>%
  mutate(
    genotype = gsub("_C", "", genotype)
  ) %>%
  mutate(
    genotype = case_when(
      is.na(genotype) ~ paste0(treatment, ifelse(is.na(volume), "", paste0("-", volume))),
      T ~ genotype
    )
  )-> short.data


small.values <- readRDS("~/projects/polq/jr_stash/final_prediction_matrix.SMALL.HBONDS.rds") %>% as.data.table() %>%
  mutate(position = position - 1) %>%
  transmute(position, deletion_size, homology, homeology = homeology0.8) %>%
  as.data.table()

short.data <- short.data %>%
  merge(small.values, by.x = c("first_del", "deletion"), by.y = c("position", "deletion_size")) %>%
  as.data.table()

col_fun_prop = readRDS("~/projects/polq/jr_stash/col_fun_prop/all.rds")

mclapply(long, function(x){
  print(x)
  readRDS(paste0(paste0("~/projects/polq/jr_stash/all_data2/", x, "_data.rds"))) -> a
  a$project <- x
  return(a)
}, mc.cores = 4) %>% rbindlist(fill = T) %>%
  mutate(
    genotype = gsub("_C", "", genotype)
  ) %>%
  mutate(
    genotype = case_when(
      is.na(genotype) ~ paste0(treatment, ifelse(is.na(volume), "", paste0("-", volume))),
      T ~ genotype
    )
  )-> long.data

long.values <- readRDS("~/projects/polq/jr_stash/final_prediction_matrix.LARGE.HBONDS.rds") %>% as.data.table() %>%
  mutate(position = position - 1) %>%
  transmute(position, deletion_size, homology, homeology = homeology0.8) %>%
  as.data.table()

long.data <- long.data %>%
  merge(long.values, by.x = c("first_del", "deletion"), by.y = c("position", "deletion_size")) %>%
  as.data.table()
long.data$group <- "IGO"

write_tsv(short.data, "~/public_html/plots/homology/short_245_data_ALL.tsv")
write_tsv(long.data, "~/public_html/plots/homology/long_503_data_ALL.tsv")

write_xlsx(long.data[cell_line == "Hek293T"], "~/public_html/plots/homology/long_503_data_Hek293T.xlsx")

write_csv(long.data[cell_line == "Hek293T"], "~/public_html/plots/homology/long_503_data_Hek293T.csv")
write_csv(long.data[cell_line == "DLD1"], "~/public_html/plots/homology/long_503_data_DLD1.csv")
write_csv(long.data[cell_line == "U2OS"], "~/public_html/plots/homology/long_503_data_U2OS.csv")


y <- function(data, x, type2, gt = NULL, class, suffix = ""){
  if(is.null(gt))
    gt <- data[(cell_line == x)]$genotype %>% unique
  data <- data[(cell_line == x) & (genotype %in% gt)] %>%
    select(c(cell_line, genotype, experiment, first_del, deletion, reads,
             total_dels, total_ins, total_indels, total_modified, total_muts, total_reads,
             homology, homeology)) %>%
    group_by(cell_line, genotype) %>%
    summarize(
      first_del,
      size = deletion,
      reads,
      homology, homeology,
      total_dels = sum(unique(total_dels)),
      total_ins = sum(unique(total_ins)),
      total_muts = sum(unique(total_muts)),
      total_indels = sum(unique(total_indels)),
      total_modified = sum(unique(total_modified)),
      total_reads = sum(unique(total_reads))
    ) %>%
    ungroup() %>%
    group_by(genotype, size, homology, homeology, total_dels, total_ins, total_indels, total_muts, total_modified, total_reads)%>%
    summarize(reads = sum(reads)) %>%
    as.data.table()
  data[, deletion_density := reads / total_dels]
  data[, log_total_dels_of_size := log(sum(reads), 10), by = c("genotype", "size")]
  data[, proportion_of_dels := reads / (sum(reads)), by = c("genotype", "size")]
  data[, fraction_of_log_bar := proportion_of_dels * log_total_dels_of_size]
  data <- data %>%
    pivot_longer(
      c(homology, homeology),
      names_to = "type", values_to = "bp"
    ) %>% filter(type == type2) %>%
    mutate(genotype = factor(genotype, levels = gt)) %>%
    as.data.table()
  #browser()
  breaks.n = c(0:ceiling(max(data[, sum(fraction_of_log_bar), by = c("genotype", "size", "type")]$V1)))
  labels.n = paste0("10<sup>", breaks.n, "</sup>")
  breaks.x = c(1, seq(from = 0, to = plyr::round_any(max(data$size), 10, ceiling), 10)[-c(1)])
  #browser()
  ggplot(data,
       aes(x = size, y = fraction_of_log_bar)) +
  geom_bar(aes(fill = factor(bp), color = factor(bp)), stat = "identity", width = 0.8, linewidth= 0.03)+
  facet_wrap( ~ genotype,
             nrow = length(unique(data$genotype)),
             strip.position = "right") +
  scale_fill_manual(values = col_fun_prop(0:max(data$bp))) +
  scale_color_manual(values = c(col_fun_prop(0:max(data$bp)))) +
  theme_bw() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(size = 8),
    axis.text.x = element_text(size = 6),
    axis.text.y = element_markdown(size = 6),
    axis.title.x = element_text(size = 9),
    axis.title.y = element_text(size = 9)
    #plot.margin = unit(c(0, 0, 0, 0), "cm")
  ) +
    scale_y_continuous(limits = c(log(1, 10), max(breaks.n)),
                     breaks = c(log(1, 10),seq(1, max(breaks.n))),
                     labels = labels.n) +
    ylab("Reads") +
    #scale_y_continuous(trans = "log10") +
    scale_x_continuous(breaks = breaks.x, expand = c(0, 0),
                       limits = c(0, max(breaks.x) + 1)) +
    xlab("Deletion Size") +
  labs(subtitle = paste0(x," 245-bp ", suffix),
       fill = "bp",
       color = "bp")-> really.small
  dir.create(paste0("~/public_html/plots/homology/READCOUNT/",x, "_", class))
  ppdf(print(really.small),
     paste0("~/public_html/plots/homology/READCOUNT/",x, "_", class, "/",
            x, "_", type2, "_",suffix, "_", class, "_", paste0(gt, collapse = "-"), "_LOG_CUSTOM.pdf"),
     width = 6,
     height = 1.5 * length(gt))
})

bad.projects.short <- c("PE250_Hek", # HEK293T 245bp GENEWIZ!!!!!
"PE150_extra") # HEK293T 245bp GENEWIZ !!!

y(short.data,
  x = "Hek293T",
  type2 = "homology",
  c("siNT", "siBRCA2", "siPolQ", "siRAD52"), class = "small", suffix = "ALL_COHORTS")
y(short.data, x = "Hek293T",
  type2 = "homeology",
  c("siNT", "siBRCA2", "siPolQ", "siRAD52"), class = "small", suffix = "ALL_COHORTS")

y(short.data, x = "Hek293T", type2 = "homology", c("siNT", "siPolH" ,"siPolL" ,"siBRCA2_PolH" ,"siBRCA2_PolL" ,"siBRCA2_PolQ" ,"siBRCA2_R52"  ,"siPolQ_R52"), class = "small", suffix = "ALL_COHORTS_OTHERSAMPLES")

y(short.data, x = "Hek293T", type2 = "homeology", c("siNT", "siPolH" ,"siPolL" ,"siBRCA2_PolH" ,"siBRCA2_PolL" ,"siBRCA2_PolQ" ,"siBRCA2_R52"  ,"siPolQ_R52"), class = "small", suffix = "ALL_COHORTS_OTHERSAMPLES")

y(short.data[!project %in% bad.projects.short], x = "Hek293T", type2 = "homology", c("siNT", "siPolH" ,"siPolL" ,"siBRCA2_PolH" ,"siBRCA2_PolL" ,"siBRCA2_PolQ" ,"siBRCA2_R52"  ,"siPolQ_R52"), class = "small", suffix = "IGO_OTHERSAMPLES")
y(short.data[!project %in% bad.projects.short], x = "Hek293T", type2 = "homeology", c("siNT", "siPolH" ,"siPolL" ,"siBRCA2_PolH" ,"siBRCA2_PolL" ,"siBRCA2_PolQ" ,"siBRCA2_R52"  ,"siPolQ_R52"), class = "small", suffix = "IGO_OTHERSAMPLES")

y(short.data[project %in% bad.projects.short], x = "Hek293T", type2 = "homology", c("siNT", "siPolH" ,"siPolL" ,"siBRCA2_PolH" ,"siBRCA2_PolL" ,"siBRCA2_PolQ" ,"siBRCA2_R52"  ,"siPolQ_R52"), class = "small", suffix = "GENEWIZ_OTHERSAMPLES")
y(short.data[project %in% bad.projects.short], x = "Hek293T", type2 = "homeology", c("siNT", "siPolH" ,"siPolL" ,"siBRCA2_PolH" ,"siBRCA2_PolL" ,"siBRCA2_PolQ" ,"siBRCA2_R52"  ,"siPolQ_R52"), class = "small", suffix = "GENEWIZ_OTHERSAMPLES")

y(short.data[!project %in% bad.projects.short],
  x = "Hek293T", type2 = "homology", c("siNT", "siBRCA2", "siPolQ", "siRAD52"), class = "small", suffix = "IGO")
y(short.data[!project %in% bad.projects.short],
  x = "Hek293T", type2 = "homeology", c("siNT", "siBRCA2", "siPolQ", "siRAD52"), class = "small", suffix = "IGO")

y(short.data[project %in% bad.projects.short],
  x = "Hek293T", type2 = "homology", c("siNT", "siBRCA2", "siPolQ", "siRAD52"), class = "small", suffix = "GENEWIZ")
y(short.data[project %in% bad.projects.short],
  x = "Hek293T", type2 = "homeology", c("siNT", "siBRCA2", "siPolQ", "siRAD52"), class = "small", suffix = "GENEWIZ")

y(short.data, x = "DLD1", type2 = "homology", class = "small", suffix = "ALL_COHORTS")
y(short.data, x = "DLD1", type2 = "homeology", class = "small", suffix = "ALL_COHORTS_80THRESH")

y(short.data, x = "U2OS", type2 = "homology", class = "small", suffix = "ALL_COHORTS")
y(short.data, x = "U2OS", type2 = "homeology", class = "small", suffix = "ALL_COHORTS_80THRESH")

y(short.data, x = "DLD1", type2 = "homeology",  class = "small")
y(short.data, x = "DLD1", type2 = "homology",  class = "small")

y(short.data, x = "U2OS", type2 = "homology",  class = "small")
y(short.data, x = "U2OS", type2 = "homeology",  class = "small")

y(long.data, x = "DLD1", type2 = "homology",  class = "long", suffix = "IGO")
y(long.data, x = "DLD1", type2 = "homeology",  class = "long", suffix = "IGO_80THRESH")

y(long.data, x = "Hek293T", type2 = "homology", c("siNT","siBRCA2", "siBRCA2_RAD52", "siBRCA2_PolQ"), class = "long", suffix = "IGO")
y(long.data, x = "Hek293T", type2 = "homeology",  c("siNT","siBRCA2", "siBRCA2_RAD52","siBRCA2_PolQ"), class = "long", suffix = "IGO_80THRESH")

y(long.data, x = "U2OS", type2 = "homology",  class = "long", suffix = "IGO")
y(long.data, x = "U2OS", type2 = "homeology",  class = "long", suffix = "IGO_80THRESH")

z <- function(data, x, type2, gt = NULL, suffix){
  if(is.null(gt))
    gt <- data[(cell_line == x)]$genotype %>% unique
  data <- data[(cell_line == x) & (genotype %in% gt)] %>%
    #filter(first_del <= 69, last_del >= 69) %>%
    select(c(cell_line, genotype, experiment, first_del, deletion, reads,
             total_dels, total_ins, total_muts, total_indels, total_modified, total_reads, homology, homeology)) %>%
    group_by(cell_line, genotype) %>%
    summarize(
      first_del,
      size = deletion,
      reads,
      homology, homeology,
      total_dels = sum(unique(total_dels)),
      total_ins = sum(unique(total_ins)),
      total_muts = sum(unique(total_muts)),
      total_indels = sum(unique(total_indels)),
      total_modified = sum(unique(total_modified)),
      total_reads = sum(unique(total_reads))
    ) %>%
    ungroup() %>%
    group_by(genotype, size, homology, homeology, total_indels, total_dels, total_modified) %>%
    summarize(reads = sum(reads)) %>%
    as.data.table()
  #data[, total_indels := sum(reads), by = c("genotype")]
  data[, deletion_density := reads / total_modified]
  #data[, log_total_dels_of_size := log(sum(reads) + 1, 10), by = c("genotype", "size")]
  #data[, proportion_of_dels := reads / (sum(reads) + 1), by = c("genotype", "size")]
  #data[, fraction_of_log_bar := proportion_of_dels * log_total_dels_of_size]
  data <- data %>%
    pivot_longer(
      c(homology, homeology),
      names_to = "type", values_to = "bp"
    ) %>%
    mutate(genotype = factor(genotype, levels = gt)) %>% as.data.table()
  #browser()
  #breaks.n = c(0:ceiling(max(data[, sum(fraction_of_log_bar), by = c("genotype", "size", "type")]$V1)))
  #labels.n = paste0("10<sup>", breaks.n, "</sup>")
                                        #breaks.x = c(1, seq(from = 0, to = plyr::round_any(max(data$size), 10), 10)[-c(1)])
  theme_density <- theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(size = 8),
    axis.text.y = element_text(size = 4),
    plot.margin = unit(c(0, 0, 0, 0), "cm")
  )
  ggplot(data[type == type2 & size == 1],
       aes(x = size, y = deletion_density)) +
  geom_bar(aes(fill = factor(bp), color = factor(bp)), stat = "identity", linewidth= 0.1) +
  facet_wrap( ~ genotype,
             nrow = length(gt)) +
  scale_fill_manual(values = col_fun_prop(0:10)) +
  scale_color_manual(values = c(rep(color0, 2), rep("black", 5))) +
  theme_bw() +
  theme_density +
    ylab("") +
  scale_x_continuous(breaks = 1, expand = c(0, 0)) +
  xlab("")-> really.small
ggplot(data[type == type2 & size <= 25 & size > 1],
       aes(x = size, y = deletion_density)) +
  geom_bar(aes(fill = factor(bp), color = factor(bp)), stat = "identity",  linewidth= 0.1) +
  facet_wrap(~ genotype,
             nrow = length(gt),
             scales = "free") +
  scale_fill_manual(values = col_fun_prop(0:11)) +
  scale_color_manual(values = c(rep(color0, 2), rep("black", 9))) +
  theme_bw() +
  theme_density +
  ylab("") +
  xlab("") +
  scale_y_continuous(labels = label_scientific()) +
  scale_x_continuous(#limits = c(2, 25),
                     breaks = c(2, 5, 10,15, 20, 25), expand = c(0, 1))-> small
ggplot(data[type == type2 & size > 25 & size <= 59],
       aes(x = size, y = deletion_density)) +
  geom_bar(aes(fill = factor(bp), color = factor(bp)), stat = "identity", linewidth= 0.1) +
  facet_wrap(~genotype,
             nrow =length(gt),
             scales = "free") + #, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = col_fun_prop(0:11)) +
  scale_color_manual(values = c(rep(color0, 2), rep("black", 9))) +
  theme_bw() +
  theme_density +
  ylab("") +
  xlab("") +
  scale_y_continuous(labels = label_scientific()) +
  scale_x_continuous(breaks = c(26, 30,35,40,45,50,55,59,60,65,70,75,80,85,90,95,100),
                     expand = c(0, 0))-> small.2
ggplot(data[type == type2 & size > 59],
       aes(x = size, y = deletion_density)) +
  geom_bar(aes(fill = factor(bp), color = factor(bp)), stat = "identity", linewidth= 0.1) +
  facet_wrap(~genotype,
             nrow = length(gt),
             scales = "free") +
  scale_fill_manual(values = col_fun_prop(0:11)) +
  scale_color_manual(values = c(rep(color0, 2), rep("black", 9))) +
  theme_bw() +
  theme_density +
  ylab("") +
  xlab("") +
  scale_y_continuous(labels = label_scientific()) +
  scale_x_continuous(breaks = c(26, 30,35,40,45,50,55,60,65,70,75,80,85,90,95,100), expand = c(0, 0))-> small.3
y <- cowplot::plot_grid(#all +theme(legend.position = "none"),
  really.small + theme(legend.position = "none"),
  small + theme(legend.position = "none"),
  small.2 + theme(legend.position = "none"),
  small.3 + theme(legend.position = "none"),
  align = "h",
  rel_widths = c(0.5,2,3,3),
  nrow = 1,
  greedy = T,
  scale = rep(0.9, 3))
  dir.create(paste0("~/public_html/plots/homology/DENSITY/",x, "_short"))
  ppdf(print(y),
     paste0("~/public_html/plots/homology/DENSITY/",x, "_short/",
            x, "_SHORT_", suffix, "_", type2, "_", paste0(gt, collapse = "-"), ".pdf"),
     width = 10,
     height = 1.5 * length(gt))
  return(data)
}

z.long <- function(data, x, type2, gt = NULL,suffix){
  if(is.null(gt))
    gt <- data[(cell_line == x)]$genotype %>% unique
  data <- data[(cell_line == x) & (genotype %in% gt)] %>%
    #filter(first_del <= 69, last_del >= 69) %>%
    select(c(cell_line, genotype, experiment, first_del, deletion, reads,
             total_dels, total_ins, total_muts, total_indels, total_modified, total_reads, homology, homeology)) %>%
    group_by(cell_line, genotype) %>%
    summarize(
      first_del,
      size = deletion,
      reads,
      homology, homeology,
      total_dels = sum(unique(total_dels)),
      total_ins = sum(unique(total_ins)),
      total_muts = sum(unique(total_muts)),
      total_indels = sum(unique(total_indels)),
      total_modified = sum(unique(total_modified)),
      total_reads = sum(unique(total_reads))
    ) %>%
    ungroup() %>%
    group_by(genotype, size, homology, homeology, total_indels, total_dels, total_modified) %>%
    summarize(reads = sum(reads)) %>%
    as.data.table()
  #data[, total_indels := sum(reads), by = c("genotype")]
  data[, deletion_density := reads / total_modified]
  #data[, log_total_dels_of_size := log(sum(reads) + 1, 10), by = c("genotype", "size")]
  #data[, proportion_of_dels := reads / (sum(reads) + 1), by = c("genotype", "size")]
  #data[, fraction_of_log_bar := proportion_of_dels * log_total_dels_of_size]
  data <- data %>%
    pivot_longer(
      c(homology, homeology),
      names_to = "type", values_to = "bp"
    ) %>%
    mutate(genotype = factor(genotype, levels = gt)) %>% as.data.table()
  theme_density <- theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(size = 8),
    axis.text.y = element_text(size = 5),
    axis.text.x = element_text(size = 5),
    plot.margin = unit(c(0, 0, 0, 0), "cm")
  )
    ggplot(data[type == type2 & size == 1],
       aes(x = size, y = deletion_density)) +
  geom_bar(aes(fill = factor(bp), color = factor(bp)), stat = "identity", linewidth= 0.1) +
  facet_wrap( ~ genotype,
             nrow = length(gt)) +
  scale_fill_manual(values = col_fun_prop(0:11)) +
  scale_color_manual(values = c(rep(color0, 2), rep("black", 10))) +
  theme_bw() +
  theme_density +
    ylab("") +
   scale_x_continuous(breaks = 1, expand = c(0, 0)) +
  xlab("")-> really.small
ggplot(data[type == type2 & size <= 35 & size > 1],
       aes(x = size, y = deletion_density)) +
  geom_bar(aes(fill = factor(bp), color = factor(bp)), stat = "identity",  linewidth= 0.1) +
  facet_wrap(~ genotype,
             nrow = length(gt),
             scales = "free") +
  scale_fill_manual(values = col_fun_prop(0:11)) +
  scale_color_manual(values = c(rep(color0, 2), rep("black", 10))) +
  theme_bw() +
  theme_density +
  ylab("") +
  xlab("") +
  scale_y_continuous(labels = label_scientific()) +
  scale_x_continuous(limits = c(1, 36),
    breaks = c(2, 5, 10,15, 20, 25, 30, 35),
    expand = c(0, 1))-> small
ggplot(data[type == type2 & size > 35 & size <= 119],
       aes(x = size, y = deletion_density)) +
  geom_bar(aes(fill = factor(bp), color = factor(bp)), stat = "identity", linewidth= 0.1) +
  facet_wrap(~genotype,
             nrow =length(gt),
             scales = "free") + #, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = col_fun_prop(0:11)) +
  scale_color_manual(values = c(rep(color0, 2), rep("black", 10))) +
  theme_bw() +
  theme_density +
  ylab("") +
  xlab("") +
  scale_y_continuous(labels = label_scientific()) +
  scale_x_continuous(limits = c(34, 120),
    breaks = c(26,30,36,40,50,60,70,80,90,100,110,119,130,140,144,150,160,170,180,190,200),
                     expand = c(0, 1))-> small.2
ggplot(data[type == type2 & size > 120],
       aes(x = size, y = deletion_density)) +
  geom_bar(aes(fill = factor(bp), color = factor(bp)), stat = "identity", linewidth= 0.1) +
  facet_wrap(~genotype,
             nrow = length(gt),
             scales = "free") +
  scale_fill_manual(values = col_fun_prop(0:11)) +
  scale_color_manual(values = c(rep(color0, 2), rep("black", 10))) +
  theme_bw() +
  theme_density +
  ylab("") +
  xlab("") +
  scale_y_continuous(labels = label_scientific()) +
  scale_x_continuous(limits = c(119, max(data$size, na.rm =T) + 1),
    breaks = c(26,30,40,50,60,70,79,80,90,100,110,120,130,140,145,150,160,170,180,190,200),
                     expand = c(0, 1))-> small.3
y <- cowplot::plot_grid(#all +theme(legend.position = "none"),
  really.small + theme(legend.position = "none"),
  small + theme(legend.position = "none"),
  small.2 + theme(legend.position = "none"),
  small.3 + theme(legend.position = "none"),
  align = "h",
  rel_widths = c(0.5, 2,3, 3),
  nrow = 1,
  scale = rep(0.9, 4))
  dir.create(paste0("~/public_html/plots/homology/DENSITY/",x, "_long"))
  ppdf(print(y),
     paste0("~/public_html/plots/homology/DENSITY/",x, "_long/",
            x, "_LONG_", suffix, "_",  type2, "_", paste0(gt, collapse = "-"), ".pdf"),
     width = 10,
     height = 1.5 * length(gt))
  return(data)
}

color0 <- "#e3e3e3"

#' even number of samples in each genotype: sample.data
sample.data <- short.data[experiment %in% c("C242a", "C243a" ,"C244a", "C226a", "C227a" ,"C228a")]

#HEKT SMALL
hek.small.all <- z(short.data, "Hek293T",
  gt = c("siNT", "siPolQ", "siRAD52", "siBRCA2"),
  type2 = "homology",
  suffix = "ALL")
hek.small.all <- z(short.data, "Hek293T",
                   gt = c("siNT", "siPolQ", "siRAD52", "siBRCA2"),
                   type2 = "homeology",
                   suffix = "ALL")
hek.small.igo <- z(short.data[group == "IGO"], "Hek293T",
                   gt = c("siNT", "siPolQ", "siRAD52", "siBRCA2"),
                   type2 = "homology",
                   suffix = "IGO")
hek.small.igo <- z(short.data[group == "IGO"], "Hek293T",
                   gt = c("siNT", "siPolQ", "siRAD52", "siBRCA2"),
                   type2 = "homeology",
                   suffix = "IGO")
hek.small.genewiz <- z(short.data[group == "GENEWIZ"], "Hek293T",
                   gt = c("siNT", "siPolQ", "siRAD52", "siBRCA2"),
                   type2 = "homology",
                   suffix = "GENEWIZ")
hek.small.genewiz <- z(short.data[group == "GENEWIZ"], "Hek293T",
                   gt = c("siNT", "siPolQ", "siRAD52", "siBRCA2"),
                   type2 = "homeology",
                   suffix = "GENEWIZ")

#HEKT SMALL SAMPLE (even number of genotypes)
hek.small.all <- z(sample.data, "Hek293T",
  gt = c("siNT", "siPolQ", "siRAD52", "siBRCA2"),
  type2 = "homology",
  suffix = "ALL_common_experiments")
hek.small.all <- z(sample.data, "Hek293T",
                   gt = c("siNT", "siPolQ", "siRAD52", "siBRCA2"),
                   type2 = "homeology",
                   suffix = "ALL_common_experiments")
hek.small.genewiz <- z(sample.data[group == "GENEWIZ"], "Hek293T",
                   gt = c("siNT", "siPolQ", "siRAD52", "siBRCA2"),
                   type2 = "homology",
                   suffix = "GENEWIZ_common_experiments")
hek.small.genewiz <- z(sample.data[group == "GENEWIZ"], "Hek293T",
                   gt = c("siNT", "siPolQ", "siRAD52", "siBRCA2"),
                   type2 = "homeology",
                   suffix = "GENEWIZ_common_experiments")


# all genewiz samples hek293T
hek.small.genewiz <- z(short.data[project %in% bad.projects], "Hek293T",
  gt = c("siNT", "siPolQ", "siRAD52"),
  type2 = "homology",
  suffix = "genewiz_all_samples")
z(short.data[project %in% bad.projects.short], "Hek293T",
  gt = c("siNT", "siPolQ", "siRAD52"),
  type2 = "homeology",
  suffix = "genewiz_all_samples")

hek.small.genewiz <- z(short.data[project %in% bad.projects], "Hek293T",
                       gt = c("siNT", "siBRCA2", "siBRCA2-PolQ", "siBRCA2_R52"),
                       type2 = "homology",
                       suffix = "genewiz_all_samples")
z(short.data[project %in% bad.projects.short], "Hek293T",
  gt = c("siNT", "siPolQ", "siRAD52"),
  type2 = "homeology",
  suffix = "genewiz_all_samples")



# DLD samples all small
dld.small.counts <- z(short.data, "DLD1",
                          "homeology", suffix = "ALL")
dld.small.counts <- z(short.data, "DLD1",
                          "homology", suffix = "ALL")

# U2OS samples all small
u2os.small.counts <- z(short.data, "U2OS",
                       "homology", suffix = "ALL")
u2os.small.counts <- z(short.data, "U2OS",
                       "homeology", suffix = "ALL")

# long samples
z.long(long.data, "Hek293T", "homology", suffix = "ALL")
z.long(long.data, "Hek293T", "homeology", suffix = "ALL")

z.long(long.data, "Hek293T", "homology", c("siNT", "siBRCA2", "siBRCA2_RAD52", "siBRCA2_PolQ"),suffix = "ALL")
z.long(long.data, "Hek293T", "homeology", c("siNT", "siBRCA2", "siBRCA2_RAD52", "siBRCA2_PolQ"), suffix = "ALL")

z.long(long.data, "DLD1", "homeology", suffix = "ALL")
z.long(long.data, "DLD1", "homology", suffix = "ALL")
z.long(long.data, "U2OS", "homology", suffix = "ALL")
z.long(long.data, "U2OS", "homeology", suffix = "ALL")


#' jrafailov Friday, Dec 06, 2024 11:21:32 AM
#'

"~/projects/polq/final_data/polq.rds" %>% readRDS -> a
"~/projects/polq/final_data/brca2.rds" %>% readRDS -> b
"~/projects/polq/final_data/rad52.rds" %>% readRDS -> c
"~/projects/polq/final_data/wt.rds" %>% readRDS -> d


data <- long.data
x <- "Hek293T"
gt = NULL
type2 = "homeology"


min = 1
max = 25
min.homo = 2
stacked = F
homo_homeo = "homo"


heatmap.fn <- function(data, x,gt = NULL, min,max,order = NULL,min.homo = 2,stacked = F,homo_homeo = "homeo", size = "long"){
  if(is.null(gt))
    gt <- data[(cell_line == x)]$genotype %>% unique
  type.n <-case_when(
    homo_homeo == "homeo" ~ "homeology",
    T ~ "homology")
  data.total <- data[(cell_line == x) & (genotype %in% gt)] %>%
    select(c(cell_line, genotype, experiment, first_del, deletion, reads,
             total_dels, total_ins, total_indels, total_modified, total_muts, total_reads,
             homology, homeology, group)) %>%
    group_by(cell_line, genotype, group) %>%
    summarize(
      first_del,
      size = deletion,
      reads,
      homology, homeology,
      total_dels = sum(unique(total_dels)),
      total_ins = sum(unique(total_ins)),
      total_muts = sum(unique(total_muts)),
      total_indels = sum(unique(total_indels)),
      total_modified = sum(unique(total_modified)),
      total_reads = sum(unique(total_reads))
    ) %>%
    ungroup() %>%
    group_by(genotype, group, size, homology, homeology, total_dels, total_modified) %>%
    summarize(reads = sum(reads)) %>%
    pivot_longer(
      c(homology, homeology),
      names_to = "type", values_to = "bp"
    ) %>%
    filter(type == type.n) %>%
    mutate(genotype = factor(genotype, levels = gt)) %>%
    as.data.table()
  data.total[, deletion_density := reads / total_modified]
  data.per.replicate <- data[(cell_line == x) & (genotype %in% gt)] %>%
    select(c(cell_line, genotype, experiment, first_del, deletion, reads,
             total_dels, total_ins, total_indels, total_modified, total_muts, total_reads,
             homology, homeology, group)) %>%
    group_by(cell_line, genotype, experiment, group) %>%
    summarize(
      first_del,
      size = deletion,
      reads,
      homology, homeology,
      total_dels = sum(unique(total_dels)),
      total_ins = sum(unique(total_ins)),
      total_muts = sum(unique(total_muts)),
      total_indels = sum(unique(total_indels)),
      total_modified = sum(unique(total_modified)),
      total_reads = sum(unique(total_reads))
    ) %>%
    ungroup() %>%
    group_by(genotype, experiment, group, size, homology, homeology, total_dels, total_modified) %>%
    summarize(reads = sum(reads)) %>%
    pivot_longer(
      c(homology, homeology),
      names_to = "type", values_to = "bp"
    ) %>%
    filter(type == type.n) %>%
    mutate(genotype = factor(genotype, levels = gt)) %>%
    as.data.table()
  data.per.replicate[, deletion_density := reads / total_modified]
  #browser()
  bp.to.plot.all <- data.total %>%
    filter(type == type.n,
           size >= min,
           size <= max,
           bp >= min.homo) %>%
    group_by(genotype, bp, group) %>%
    summarize(deletion_density = sum(deletion_density),
              total_dels = unique(total_dels)) %>%
    as.data.table()
  bp.to.plot.rep <- data.per.replicate %>%
    filter(type == type.n,
           size >= min,
           size <= max,
           bp >= min.homo) %>%
    group_by(genotype, bp, experiment, group) %>%
    summarize(deletion_density = sum(deletion_density),
              total_dels = unique(total_dels)) %>%
    mutate(score = deletion_density) %>%
    as.data.table()
  scores.per.replicate <- bp.to.plot.rep %>%
    mutate(score = bp * deletion_density) %>%
    group_by(genotype, experiment, group) %>%
    summarize(score = sum(score),
              deletion_density = sum(deletion_density),
              total_dels = unique(total_dels)) %>%
    ungroup %>%
    group_by(genotype) %>%
    mutate(
      sd.score = sd(score),
      mean.score = mean(score),
      se.score = sd.score / sqrt(n()),
      sd.deletion_density = sd(deletion_density),
      mean.deletion_density = mean(deletion_density),
      se.deletion_density = sd.deletion_density / sqrt(n())
    ) %>%
    as.data.table()
  scores.per.replicate.group <- bp.to.plot.rep %>%
    mutate(score = bp * deletion_density) %>%
    group_by(genotype, experiment, group) %>%
    summarize(score = sum(score),
              deletion_density = sum(deletion_density),
              total_dels = unique(total_dels)) %>%
    as.data.table()
  control.mean.scores <-  bp.to.plot.rep[as.numeric(genotype) == 1]%>%
    mutate(score = bp * deletion_density) %>%
    group_by(genotype, experiment, group) %>%
    summarize(score = mean(sum(score)),
              deletion_density = mean(sum(deletion_density)),
              total_dels = sum(unique(total_dels))) %>%
    as.data.table()
#browser()
  merge.norm <- control.mean.scores[,-c("genotype")]
  names(merge.norm) <- paste0("control_", names(merge.norm))
  norm.scores.per.replicate <- bp.to.plot.rep %>%
    mutate(score = bp * deletion_density) %>%
    group_by(genotype, experiment, group) %>%
    summarize(score = sum(score),
              deletion_density = sum(deletion_density),
              total_dels = unique(total_dels))%>%
  merge(merge.norm, by.x = "experiment", by.y = "control_experiment")  %>%
    mutate(norm_score = score / control_score,
           norm_deletion_density = deletion_density / control_deletion_density) %>%
    ungroup %>%
    group_by(genotype) %>%
    mutate(
      sd.norm_score = sd(norm_score),
      mean.norm_score = mean(norm_score),
      se.norm_score = sd.norm_score / sqrt(n()),
      sd.norm_deletion_density = sd(norm_deletion_density),
      mean.norm_deletion_density = mean(norm_deletion_density),
      se.norm_deletion_density = sd.norm_deletion_density / sqrt(n())
    ) %>%
    as.data.table()
  #browser()
  stats.score <- scores.per.replicate %>%
    as.data.table() %>%
    t_test(score ~ genotype ) %>%
    #filter(!p.adj.signif == "ns") %>%
    add_y_position(fun = "max", step.increase = 0.12) %>%
    filter(!p.adj.signif == "ns")
  #browser()
  tryCatch(
    expr = {
      stats.density <- scores.per.replicate %>%
        as.data.table() %>%
        t_test(deletion_density ~ genotype ) %>%
                                        #filter(!p.adj.signif == "ns") %>%
        add_y_position(fun = "max", step.increase = 0.12)      %>%
        filter(!p.adj.signif == "ns")
    },
    error = function(e){
      print("hey")
      stats.score %>%
        transmute(group1, group2, p.adj.signif = "ns", y.position) %>%
        filter(!p.adj.signif == "ns") -> stats.density
    }
  )
  if(!exists("stats.density")){
    stats.score %>%
      transmute(group1, group2, p.adj.signif = "ns", y.position) %>%
      filter(!p.adj.signif == "ns")-> stats.density
  }
  stats.norm.score <- norm.scores.per.replicate %>%
    as.data.table() %>%
    t_test(norm_score ~ genotype) %>%
    add_y_position(fun = "max") %>%
    filter(!p.adj.signif == "ns")
  tryCatch(
    expr = {
      stats.norm.density <- norm.scores.per.replicate %>%
        as.data.table() %>%
        t_test(norm_deletion_density ~ genotype ) %>%
        #filter(!p.adj.signif == "ns") %>%
        add_y_position(fun = "max", step.increase = 0.12) %>%
        filter(!p.adj.signif == "ns")
    },
    error = function(e){
      print("hey2")
      stats.norm.score %>%
        transmute(group1, group2, p.adj.signif = "ns", y.position) %>%
        filter(!p.adj.signif == "ns")-> stats.norm.density
    }
  )
  if(!exists("stats.norm.density")){
    stats.norm.score %>%
      transmute(group1, group2, p.adj.signif = "ns", y.position) %>%
      filter(!p.adj.signif == "ns")-> stats.norm.density
  }
  browser()
  ggplot(bp.to.plot.all, aes(x = genotype, y = deletion_density)) +
    geom_bar(aes(fill = factor(bp), color = factor(bp)), stat = "identity", lwd = 0.1, width = 0.4)+
    #stat_pvalue_manual(stats, label = "p.adj.signif") +
    scale_fill_manual(values = col_fun_prop(sort(unique(bp.to.plot.all$bp)))) +
    scale_color_manual(values = rep("black", length(unique(bp.to.plot.all$bp)))) +
    theme_pubr() +
    theme(panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.major.y = element_line(size = .1, color = "gray"),
          panel.grid.minor.y = element_blank(),
          strip.background = element_blank(),
          strip.text = element_text(size = 8),
          axis.text.x = element_text(size = 8, angle = 45, vjust = 0.5, hjust = 0.5),
          axis.ticks = element_blank(),
          legend.position = "right",
          plot.title = element_text(size = 11, face = "bold"),
          plot.subtitle = element_text(face = "italic"),
          #plot.margin = unit(c(0, 0, 0, 0), "cm")
          ) +
    facet_grid(cols = vars(group))+
    labs(title = paste0(x, "-", size,": ",  max, " >= deletion size >= ", min, ", at least ", min.homo, " bp ", type.n),
         subtitle = "cumulative deletion density per genotype",
         color = "bp", fill = "bp") +
    ylab("% of alterations") +
    xlab("")-> all
  ggplot(bp.to.plot.rep, aes(x = experiment, y = deletion_density)) +
    geom_bar(aes(fill = factor(bp), color = factor(bp)), stat = "identity", lwd = 0.1) +
    #stat_pvalue_manual(stats.density, label = "p.adj.signif") +
    scale_fill_manual(values = col_fun_prop(sort(unique(bp.to.plot.rep$bp)))) +
    scale_color_manual(values = rep("black", length(unique(bp.to.plot.rep$bp)))) +
    facet_grid( ~ genotype + group,
              scales = "free_x",
              #strip.position = "bottom",
              space = "free_x") +
               #space = "free_x") +
    theme_bw() +
    theme(panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.major.y = element_line(size = .1, color = "gray"),
          panel.grid.minor.y = element_blank(),
          strip.background = element_blank(), #element_rect(color = "black", fill = NA, size = 0.03),
          strip.text = element_text(size = 6, angle = 0),
          axis.text.x = element_text(size = 6, angle = 90),
          axis.ticks = element_blank(),
          strip.clip = "off",
          legend.position = "none",
          plot.subtitle = element_text(face = "italic")
                                        #plot.margin = unit(c(0, 0, 0, 0), "cm")
          ) +
    labs(subtitle = "Alteration proportions per replicate") +
    ylab("% of deletions") +
    xlab("")-> all.rep
  theme.right <- theme(panel.grid.major.x = element_blank(),
                       panel.grid.minor.x = element_blank(),
                       panel.grid.major.y = element_blank(),#element_line(size = .1, color = "black"),
                       panel.grid.minor.y = element_blank(),
                       axis.ticks.x = element_blank(),
                       axis.text.x = element_text(angle = 30, size = 5, vjust = 0.5, hjust = 0.5),
                       strip.background = element_blank(),
                       strip.text = element_text(size = 8),
                       legend.position = "none"
                                        #plot.margin = unit(c(0, 0, 0, 0), "cm")
                       )
  ggplot(scores.per.replicate, aes(x = genotype, y = score)) +
    geom_bar(stat = "summary", color = "black", fill = "skyblue", width = 0.7, lwd = 0.1) +
    geom_jitter(aes(color = group), size = 1, width = 0.3, height = 0) +
    geom_errorbar(aes(ymin = mean.score - sd.score, ymax = mean.score + sd.score),
                  width = 0.4, color = "black", lwd = 0.2) +
    stat_pvalue_manual(stats.score, label = "p.adj.signif", label.size = 2) +
    theme_pubr() +
    theme.right +
    #labs(title = paste0("Deletion sizes ", min, "-", max, ", inclusive")) +
    ylab("Weighted") +
    xlab("")-> all.score
  ggplot(norm.scores.per.replicate, aes(x = genotype, y = norm_score)) +
    geom_bar(stat = "summary", color = "black", fill = "#D3D3FF", width = 0.7, lwd = 0.1) +
    geom_jitter(aes(color = group), size = 1, width = 0.3, height = 0) +
    geom_errorbar(aes(ymin = mean.norm_score - sd.norm_score,
                      ymax = mean.norm_score + sd.norm_score),
                  width = 0.4, color = "black", lwd = 0.2) +
    stat_pvalue_manual(stats.norm.score, label = "p.adj.signif", label.size = 2) +
    theme_pubr() +
    theme.right +
    #labs(title = paste0("Deletion sizes ", min, "-", max, ", inclusive")) +
    ylab("Normalized/weighted") +
    xlab("")-> all.score.norm
  ggplot(scores.per.replicate, aes(x = genotype, y = deletion_density)) +
    geom_bar(stat = "summary", color = "black", fill = "skyblue", width = 0.7, lwd = 0.1) +
    geom_jitter(aes(x = genotype, y = deletion_density, color = group),
                data = scores.per.replicate,size = 1, width = 0.3, height = 0, inherit.aes = F) +
    geom_errorbar(aes(ymin = mean.deletion_density - sd.deletion_density,
                      ymax = mean.deletion_density + sd.deletion_density),
                  width = 0.4, color = "black", lwd = 0.2) +
    stat_pvalue_manual(stats.density, label = "p.adj.signif", label.size = 2) +
    theme_pubr() +
    theme.right +
    #labs(title = paste0("Deletion sizes ", min, "-", max, ", inclusive")) +
    ylab("deletion proportions") +
    #scale_y_continuous(expand = c(0,NA), limits = c(0,NA))+
    xlab("")-> all.density
  ggplot(norm.scores.per.replicate, aes(x = genotype, y = norm_deletion_density)) +
    geom_bar(stat = "summary", color = "black", fill = "#D3D3FF", width = 0.7, lwd = 0.1) +
    geom_jitter(aes(x = genotype, y = norm_deletion_density, color = group),
                data = norm.scores.per.replicate,size = 1, width = 0.3, height = 0, inherit.aes = F) +
    geom_errorbar(aes(ymin = mean.norm_deletion_density - sd.norm_deletion_density,
                      ymax = mean.norm_deletion_density + sd.norm_deletion_density),
                  width = 0.4, color = "black", lwd = 0.2) +
    stat_pvalue_manual(stats.norm.density, label = "p.adj.signif", label.size = 2) +
    theme_pubr() +
    theme.right +
    #labs(title = paste0("Deletion sizes ", min, "-", max, ", inclusive")) +
    ylab("Normalized") +
    #scale_y_continuous(expand = c(0,NA), limits = c(0,NA))+
    xlab("")-> all.density.norm
  #browser()
  first_col = plot_grid(all, all.rep, ncol = 1)
  second_col = plot_grid(all.density, all.density.norm, all.score, all.score.norm, ncol = 1, align = "v")
  #browser()
  xxx <- plot_grid(#all +theme(legend.position = "none"),
    first_col,
    second_col,
    #ncol = 1
    #ncol = 1,
    #align = "h",
    rel_widths = c(1.2,0.8),
    #rel_heights = c(1,1.2),
    scale = c(1, 1),
    nrow = 1)
    #browser()
  dir.to.save <- paste0("~/public_html/plots/homology/HEATMAPS/", x, "_", size, "/")
  if(!dir.exists(dir.to.save))
    dir.create(dir.to.save)
  ppdf(plot(xxx), paste0(dir.to.save,
                         x, "_", size,"_",  min, '-', max, "size_", min.homo,"bp",homo_homeo, "___",paste0(gt, collapse = "-"), ".pdf"),
       width = 8.5, height = 11)
}

data <- short.data
x <- "Hek293T"
gt = c("siNT", "siPolQ", "siBRCA2", "siRAD52")
min = 3
max = 25
min.homo = 2
stacked = F
homo_homeo = "homo"
size <- "size"

heatmap.fn(short.data,
           "Hek293T",
           c("siNT", "siPolQ", "siBRCA2", "siRAD52"), 3, 25, min.homo = 3, stacked = F, homo_homeo = "homo")




debug(heatmap.fn)

undebug(heatmap.fn)

short.data.hold <-  short.data

short.data <- short.data.hold

heatmap.fn(short.data, "Hek293T", c("siNT", "siPolQ", "siBRCA2", "siRAD52"), 3, 25, min.homo = 3, stacked = F, homo_homeo = "homo")
heatmap.fn(short.data, "Hek293T", c("siNT", "siPolQ", "siBRCA2", "siRAD52"), 3, 25, min.homo = 3, stacked = F, homo_homeo = "homeo")
heatmap.fn(short.data, "Hek293T", c("siNT", "siPolQ", "siBRCA2", "siRAD52"), 25, 100, min.homo = 3, stacked = F, homo_homeo = "homo")
heatmap.fn(short.data, "Hek293T", c("siNT", "siPolQ", "siBRCA2", "siRAD52"), 25, 100, min.homo = 3, stacked = F, homo_homeo = "homeo")
heatmap.fn(short.data, "Hek293T", c("siNT", "siPolQ", "siBRCA2", "siRAD52"), 50, 100, min.homo = 3, stacked = F, homo_homeo = "homo")
heatmap.fn(short.data, "Hek293T", c("siNT", "siPolQ", "siBRCA2", "siRAD52"), 50, 100, min.homo = 3, stacked = F, homo_homeo = "homeo")
heatmap.fn(short.data, "Hek293T", c("siNT", "siPolQ", "siBRCA2", "siRAD52"), 1, 25, min.homo = 2, stacked = F, homo_homeo = "homo")
heatmap.fn(short.data, "Hek293T", c("siNT", "siPolQ", "siBRCA2", "siRAD52"), 1, 25, min.homo = 2, stacked = F, homo_homeo = "homeo")
heatmap.fn(short.data, "Hek293T", c("siNT", "siPolQ", "siBRCA2", "siRAD52"), 1, 25, min.homo = 0, stacked = F, homo_homeo = "homo")
heatmap.fn(short.data, "Hek293T", c("siNT", "siPolQ", "siBRCA2", "siRAD52"), 1, 25, min.homo = 0, stacked = F, homo_homeo = "homeo")
heatmap.fn(short.data, "Hek293T", c("siNT", "siPolQ", "siBRCA2", "siRAD52"), 1, 100, min.homo = 0, stacked = F, homo_homeo = "homo")
heatmap.fn(short.data, "Hek293T", c("siNT", "siPolQ", "siBRCA2", "siRAD52"), 1, 100, min.homo = 0, stacked = F, homo_homeo = "homeo")
heatmap.fn(short.data, "Hek293T", c("siNT", "siPolQ", "siBRCA2", "siRAD52"), 1, 100, min.homo = 2, stacked = F, homo_homeo = "homo")
heatmap.fn(short.data, "Hek293T", c("siNT", "siPolQ", "siBRCA2", "siRAD52"), 1, 100, min.homo = 2, stacked = F, homo_homeo = "homeo")
heatmap.fn(short.data, "Hek293T", c("siNT", "siPolQ", "siBRCA2", "siRAD52"), 26, 100, min.homo = 0, stacked = F, homo_homeo = "homo")
heatmap.fn(short.data, "Hek293T", c("siNT", "siPolQ", "siBRCA2", "siRAD52"), 26, 100, min.homo = 0, stacked = F, homo_homeo = "homeo")
heatmap.fn(short.data, "Hek293T", c("siNT", "siPolQ", "siBRCA2", "siRAD52"), 26, 100, min.homo = 2, stacked = F, homo_homeo = "homo")
heatmap.fn(short.data, "Hek293T", c("siNT", "siPolQ", "siBRCA2", "siRAD52"), 26, 100, min.homo = 2, stacked = F, homo_homeo = "homeo")
heatmap.fn(short.data, "Hek293T", c("siNT", "siPolQ", "siBRCA2", "siRAD52"), 50, 100, min.homo = 0, stacked = F, homo_homeo = "homo")
heatmap.fn(short.data, "Hek293T", c("siNT", "siPolQ", "siBRCA2", "siRAD52"), 50, 100, min.homo = 0, stacked = F, homo_homeo = "homeo")
heatmap.fn(short.data, "Hek293T", c("siNT", "siPolQ", "siBRCA2", "siRAD52"), 50, 100, min.homo = 2, stacked = F, homo_homeo = "homo")
heatmap.fn(short.data, "Hek293T", c("siNT", "siPolQ", "siBRCA2", "siRAD52"), 50, 100, min.homo = 2, stacked = F, homo_homeo = "homeo")
heatmap.fn(short.data, "Hek293T", c("siNT", "siPolQ", "siPolL", "siPolH"), 1, 25, min.homo = 2, stacked = F, homo_homeo = "homo")
heatmap.fn(short.data, "Hek293T", c("siNT", "siPolQ", "siPolL", "siPolH"), 1, 25, min.homo = 2, stacked = F, homo_homeo = "homeo")
heatmap.fn(short.data, "Hek293T", c("siNT", "siPolQ", "siPolL", "siPolH"), 1, 25, min.homo = 0, stacked = F, homo_homeo = "homo")
heatmap.fn(short.data, "Hek293T", c("siNT", "siPolQ", "siPolL", "siPolH"), 1, 25, min.homo = 0, stacked = F, homo_homeo = "homeo")
heatmap.fn(short.data, "Hek293T", c("siNT", "siPolQ", "siPolL", "siPolH"), 1, 100, min.homo = 0, stacked = F, homo_homeo = "homo")
heatmap.fn(short.data, "Hek293T", c("siNT", "siPolQ", "siPolL", "siPolH"), 1, 100, min.homo = 0, stacked = F, homo_homeo = "homeo")
heatmap.fn(short.data, "Hek293T", c("siNT", "siPolQ", "siPolL", "siPolH"), 1, 100, min.homo = 2, stacked = F, homo_homeo = "homeo")
heatmap.fn(short.data, "Hek293T", c("siNT", "siPolQ", "siPolL", "siPolH"), 1, 100, min.homo = 2, stacked = F, homo_homeo = "homeo")
heatmap.fn(short.data, "Hek293T", c("siNT", "siPolQ", "siPolL", "siPolH"), 26, 100, min.homo = 0, stacked = F, homo_homeo = "homo")
heatmap.fn(short.data, "Hek293T", c("siNT", "siPolQ", "siPolL", "siPolH"), 26, 100, min.homo = 0, stacked = F, homo_homeo = "homeo")
heatmap.fn(short.data, "Hek293T", c("siNT", "siPolQ", "siPolL", "siPolH"), 26, 100, min.homo = 2, stacked = F, homo_homeo = "homo")
heatmap.fn(short.data, "Hek293T", c("siNT", "siPolQ", "siPolL", "siPolH"), 26, 100, min.homo = 2, stacked = F, homo_homeo = "homeo")
heatmap.fn(short.data, "Hek293T", c("siNT", "siPolQ", "siPolL", "siPolH"), 50, 100, min.homo = 0, stacked = F, homo_homeo = "homo")
heatmap.fn(short.data, "Hek293T", c("siNT", "siPolQ", "siPolL", "siPolH"), 50, 100, min.homo = 0, stacked = F, homo_homeo = "homeo")
heatmap.fn(short.data, "Hek293T", c("siNT", "siPolQ", "siPolL", "siPolH"), 50, 100, min.homo = 2, stacked = F, homo_homeo = "homo")
heatmap.fn(short.data, "Hek293T", c("siNT", "siPolQ", "siPolL", "siPolH"), 50, 100, min.homo = 2, stacked = F, homo_homeo = "homeo")
heatmap.fn(short.data, "Hek293T", c("siNT", "siPolQ", "siPolL", "siPolH"), 3, 25, min.homo = 3, stacked = F, homo_homeo = "homo")
heatmap.fn(short.data, "Hek293T", c("siNT", "siPolQ", "siPolL", "siPolH"), 3, 25, min.homo = 3, stacked = F, homo_homeo = "homeo")
heatmap.fn(short.data, "Hek293T", c("siNT", "siPolQ", "siPolL", "siPolH"), 25, 100, min.homo = 3, stacked = F, homo_homeo = "homo")
heatmap.fn(short.data, "Hek293T", c("siNT", "siPolQ", "siPolL", "siPolH"), 25, 100, min.homo = 3, stacked = F, homo_homeo = "homeo")
heatmap.fn(short.data, "Hek293T", c("siNT", "siPolQ", "siPolL", "siPolH"), 50, 100, min.homo = 3, stacked = F, homo_homeo = "homo")
heatmap.fn(short.data, "Hek293T", c("siNT", "siPolQ", "siPolL", "siPolH"), 50, 100, min.homo = 3, stacked = F, homo_homeo = "homeo")

od <- c("siNT", "siRAD52","siPolQ", "siBRCA2", "siBRCA2_R52","siBRCA2_PolQ")

debug(heatmap.fn)

heatmap.fn(short.data, "Hek293T", od, 5, 25, min.homo = 3, stacked = F, homo_homeo = "homo")
heatmap.fn(short.data, "Hek293T", od, 25, 100, min.homo = 3, stacked = F, homo_homeo = "homo")

heatmap.fn(short.data, "Hek293T", od, 5, 25, min.homo = 3, stacked = F, homo_homeo = "homeo")
heatmap.fn(short.data, "Hek293T", od, 25, 100, min.homo = 3, stacked = F, homo_homeo = "homeo")

heatmap.fn(short.data[group == "GENEWIZ"], "Hek293T", od, 5, 25, min.homo = 3, stacked = F, homo_homeo = "homo")
heatmap.fn(short.data[group == "GENEWIZ"], "Hek293T", od, 25, 100, min.homo = 3, stacked = F, homo_homeo = "homo")



heatmap.fn(short.data, "Hek293T", od, 1, 25, min.homo = 2, stacked = F, homo_homeo = "homo")
heatmap.fn(short.data, "Hek293T", od, 1, 25, min.homo = 2, stacked = F, homo_homeo = "homo")

heatmap.fn(short.data[first_del == 69], "Hek293T", od, 12, 12, min.homo = 5, stacked = F, homo_homeo = "homo")
heatmap.fn(short.data[first_del == 69], "Hek293T", od, 12, 12, min.homo = 5, stacked = F, homo_homeo = "homo")




heatmap.fn(short.data, "Hek293T", c("siNT", "siBRCA2", "siBRCA2_PolH", "siBRCA2_PolL", "siBRCA2_PolQ"), 1, 25, min.homo = 2, stacked = F, homo_homeo = "homo")
heatmap.fn(short.data, "Hek293T", c("siNT", "siBRCA2", "siBRCA2_PolH", "siBRCA2_PolL", "siBRCA2_PolQ"), 1, 25, min.homo = 2, stacked = F, homo_homeo = "homeo")
heatmap.fn(short.data, "Hek293T", c("siNT", "siBRCA2", "siBRCA2_PolH", "siBRCA2_PolL", "siBRCA2_PolQ"), 1, 25, min.homo = 0, stacked = F, homo_homeo = "homo")
heatmap.fn(short.data, "Hek293T", c("siNT", "siBRCA2", "siBRCA2_PolH", "siBRCA2_PolL", "siBRCA2_PolQ"), 1, 25, min.homo = 0, stacked = F, homo_homeo = "homeo")
heatmap.fn(short.data, "Hek293T", c("siNT", "siBRCA2", "siBRCA2_PolH", "siBRCA2_PolL", "siBRCA2_PolQ"), 1, 100, min.homo = 0, stacked = F, homo_homeo = "homo")
heatmap.fn(short.data, "Hek293T", c("siNT", "siBRCA2", "siBRCA2_PolH", "siBRCA2_PolL", "siBRCA2_PolQ"), 1, 100, min.homo = 0, stacked = F, homo_homeo = "homeo")
heatmap.fn(short.data, "Hek293T", c("siNT", "siBRCA2", "siBRCA2_PolH", "siBRCA2_PolL", "siBRCA2_PolQ"), 1, 100, min.homo = 2, stacked = F, homo_homeo = "homo")
heatmap.fn(short.data, "Hek293T", c("siNT", "siBRCA2", "siBRCA2_PolH", "siBRCA2_PolL", "siBRCA2_PolQ"), 1, 100, min.homo = 2, stacked = F, homo_homeo = "homeo")
heatmap.fn(short.data, "Hek293T", c("siNT", "siBRCA2", "siBRCA2_PolH", "siBRCA2_PolL", "siBRCA2_PolQ"), 26, 100, min.homo = 0, stacked = F, homo_homeo = "homo")
heatmap.fn(short.data, "Hek293T", c("siNT", "siBRCA2", "siBRCA2_PolH", "siBRCA2_PolL", "siBRCA2_PolQ"), 26, 100, min.homo = 0, stacked = F, homo_homeo = "homeo")
heatmap.fn(short.data, "Hek293T", c("siNT", "siBRCA2", "siBRCA2_PolH", "siBRCA2_PolL", "siBRCA2_PolQ"), 26, 100, min.homo = 2, stacked = F, homo_homeo = "homo")
heatmap.fn(short.data, "Hek293T", c("siNT", "siBRCA2", "siBRCA2_PolH", "siBRCA2_PolL", "siBRCA2_PolQ"), 26, 100, min.homo = 2, stacked = F, homo_homeo = "homeo")
heatmap.fn(short.data, "Hek293T", c("siNT", "siBRCA2", "siBRCA2_PolH", "siBRCA2_PolL", "siBRCA2_PolQ"), 50, 100, min.homo = 0, stacked = F, homo_homeo = "homo")
heatmap.fn(short.data, "Hek293T", c("siNT", "siBRCA2", "siBRCA2_PolH", "siBRCA2_PolL", "siBRCA2_PolQ"), 50, 100, min.homo = 0, stacked = F, homo_homeo = "homeo")
heatmap.fn(short.data, "Hek293T", c("siNT", "siBRCA2", "siBRCA2_PolH", "siBRCA2_PolL", "siBRCA2_PolQ"), 50, 100, min.homo = 2, stacked = F, homo_homeo = "homo")
heatmap.fn(short.data, "Hek293T", c("siNT", "siBRCA2", "siBRCA2_PolH", "siBRCA2_PolL", "siBRCA2_PolQ"), 50, 100, min.homo = 2, stacked = F, homo_homeo = "homeo")
heatmap.fn(short.data, "Hek293T", c("siNT", "siBRCA2", "siBRCA2_PolH", "siBRCA2_PolL", "siBRCA2_PolQ"), 3, 25, min.homo = 3, stacked = F, homo_homeo = "homo")
heatmap.fn(short.data, "Hek293T", c("siNT", "siBRCA2", "siBRCA2_PolH", "siBRCA2_PolL", "siBRCA2_PolQ"), 3, 25, min.homo = 3, stacked = F, homo_homeo = "homeo")
heatmap.fn(short.data, "Hek293T", c("siNT", "siBRCA2", "siBRCA2_PolH", "siBRCA2_PolL", "siBRCA2_PolQ"), 25, 100, min.homo = 3, stacked = F, homo_homeo = "homo")
heatmap.fn(short.data, "Hek293T", c("siNT", "siBRCA2", "siBRCA2_PolH", "siBRCA2_PolL", "siBRCA2_PolQ"), 25, 100, min.homo = 3, stacked = F, homo_homeo = "homeo")
heatmap.fn(short.data, "Hek293T", c("siNT", "siBRCA2", "siBRCA2_PolH", "siBRCA2_PolL", "siBRCA2_PolQ"), 50, 100, min.homo = 3, stacked = F, homo_homeo = "homo")
heatmap.fn(short.data, "Hek293T", c("siNT", "siBRCA2", "siBRCA2_PolH", "siBRCA2_PolL", "siBRCA2_PolQ"), 50, 100, min.homo = 3, stacked = F, homo_homeo = "homeo")

od <- c("WT", "PolQ-PolDel", "PolQ-HelDel", "PolQ-PolDHelD")

heatmap.fn(short.data, "DLD1", od, min = 1, max = 25, min.homo = 2, stacked = F, homo_homeo = "homeo")
heatmap.fn(short.data, "DLD1", od, min = 1, max = 25, min.homo = 2, stacked = F, homo_homeo = "homo")
heatmap.fn(short.data, "DLD1", od, min = 1, max = 25, min.homo = 0, stacked = F, homo_homeo = "homeo")
heatmap.fn(short.data, "DLD1", od, min = 1, max = 25, min.homo = 0, stacked = F, homo_homeo = "homo")
heatmap.fn(short.data, "DLD1", od, min = 1, max = 100, min.homo = 2, stacked = F, homo_homeo = "homeo")
heatmap.fn(short.data, "DLD1", od, min = 1, max = 100, min.homo = 2, stacked = F, homo_homeo = "homo")
heatmap.fn(short.data, "DLD1", od, min = 1, max = 100, min.homo = 0, stacked = F, homo_homeo = "homeo")
heatmap.fn(short.data, "DLD1", od, min = 1, max = 100, min.homo = 0, stacked = F, homo_homeo = "homo")
heatmap.fn(short.data, "DLD1", od, min = 25, max = 100, min.homo = 2, stacked = F, homo_homeo = "homeo")
heatmap.fn(short.data, "DLD1", od, min = 25, max = 100, min.homo = 2, stacked = F, homo_homeo = "homo")
heatmap.fn(short.data, "DLD1", od, min = 25, max = 100, min.homo = 0, stacked = F, homo_homeo = "homeo")
heatmap.fn(short.data, "DLD1", od, min = 25, max = 100, min.homo = 0, stacked = F, homo_homeo = "homo")
heatmap.fn(short.data, "DLD1", od, min = 50, max = 100, min.homo = 2, stacked = F, homo_homeo = "homeo")
heatmap.fn(short.data, "DLD1", od, min = 50, max = 100, min.homo = 2, stacked = F, homo_homeo = "homo")
heatmap.fn(short.data, "DLD1", od, min = 50, max = 100, min.homo = 0, stacked = F, homo_homeo = "homeo")
heatmap.fn(short.data, "DLD1", od, min = 50, max = 100, min.homo = 0, stacked = F, homo_homeo = "homo")
heatmap.fn(short.data, "DLD1", od, min = 3, max = 25, min.homo = 3, stacked = F, homo_homeo = "homeo")
heatmap.fn(short.data, "DLD1", od, min = 3, max = 25, min.homo = 3, stacked = F, homo_homeo = "homo")
heatmap.fn(short.data, "DLD1", od, min = 50, max = 100, min.homo = 3, stacked = F, homo_homeo = "homeo")
heatmap.fn(short.data, "DLD1", od, min = 50, max = 100, min.homo = 3, stacked = F, homo_homeo = "homo")
heatmap.fn(short.data, "DLD1", od, min = 25, max = 100, min.homo = 3, stacked = F, homo_homeo = "homeo")
heatmap.fn(short.data, "DLD1", od, min = 25, max = 100, min.homo = 3, stacked = F, homo_homeo = "homo")

od <- c("DMSO", "NVB-50uM", "NVB-100uM", "ART558-10uM", "ART558-25uM",  "RP6685-10uM", "RP9913-10uM")



heatmap.fn(short.data, "U2OS", od, min = 1, max = 25, min.homo = 2, stacked = F, homo_homeo = "homeo")
heatmap.fn(short.data, "U2OS", od, min = 1, max = 25, min.homo = 2, stacked = F, homo_homeo = "homo")
heatmap.fn(short.data, "U2OS", od, min = 1, max = 25, min.homo = 0, stacked = F, homo_homeo = "homeo")
heatmap.fn(short.data, "U2OS", od, min = 1, max = 25, min.homo = 0, stacked = F, homo_homeo = "homo")
heatmap.fn(short.data, "U2OS", od, min = 1, max = 100, min.homo = 2, stacked = F, homo_homeo = "homeo")
heatmap.fn(short.data, "U2OS", od, min = 1, max = 100, min.homo = 2, stacked = F, homo_homeo = "homo")
heatmap.fn(short.data, "U2OS", od, min = 1, max = 100, min.homo = 0, stacked = F, homo_homeo = "homeo")
heatmap.fn(short.data, "U2OS", od, min = 1, max = 100, min.homo = 0, stacked = F, homo_homeo = "homo")
heatmap.fn(short.data, "U2OS", od, min = 25, max = 100, min.homo = 2, stacked = F, homo_homeo = "homeo")
heatmap.fn(short.data, "U2OS", od, min = 25, max = 100, min.homo = 2, stacked = F, homo_homeo = "homo")
heatmap.fn(short.data, "U2OS", od, min = 25, max = 100, min.homo = 0, stacked = F, homo_homeo = "homeo")
heatmap.fn(short.data, "U2OS", od, min = 25, max = 100, min.homo = 0, stacked = F, homo_homeo = "homo")
heatmap.fn(short.data, "U2OS", od, min = 50, max = 100, min.homo = 2, stacked = F, homo_homeo = "homeo")
heatmap.fn(short.data, "U2OS", od, min = 50, max = 100, min.homo = 2, stacked = F, homo_homeo = "homo")
heatmap.fn(short.data, "U2OS", od, min = 50, max = 100, min.homo = 0, stacked = F, homo_homeo = "homeo")
heatmap.fn(short.data, "U2OS", od, min = 50, max = 100, min.homo = 0, stacked = F, homo_homeo = "homo")
heatmap.fn(short.data, "U2OS", od, min = 50, max = 100, min.homo = 3, stacked = F, homo_homeo = "homeo")
heatmap.fn(short.data, "U2OS", od, min = 50, max = 100, min.homo = 3, stacked = F, homo_homeo = "homo")
heatmap.fn(short.data, "U2OS", od, min = 25, max = 100, min.homo = 3, stacked = F, homo_homeo = "homeo")
heatmap.fn(short.data, "U2OS", od, min = 25, max = 100, min.homo = 3, stacked = F, homo_homeo = "homo")
heatmap.fn(short.data, "U2OS", od, min = 3, max = 25, min.homo = 3, stacked = F, homo_homeo = "homeo")
heatmap.fn(short.data, "U2OS", od, min = 3, max = 25, min.homo = 3, stacked = F, homo_homeo = "homo")


c("siNT", "siBRCA2", "siBRCA2_RAD52", "siBRCA2_PolQ", "siBRCA2_PolL", "siBRCA2_PolH") -> order

heatmap.fn(long.data, "Hek293T", gt = order, min = 1, max = 25, min.homo = 2, stacked = F, homo_homeo = "homeo")
heatmap.fn(long.data, "Hek293T", gt = order, min = 1, max = 25, min.homo = 2, stacked = F, homo_homeo = "homo")
heatmap.fn(long.data, "Hek293T", gt = order, min = 1, max = 25, min.homo = 0, stacked = F, homo_homeo = "homeo")
heatmap.fn(long.data, "Hek293T", gt = order, min = 1, max = 25, min.homo = 0, stacked = F, homo_homeo = "homo")
heatmap.fn(long.data, "Hek293T", gt = order, min = 1, max = 100, min.homo = 2, stacked = F, homo_homeo = "homeo")
heatmap.fn(long.data, "Hek293T", gt = order, min = 1, max = 100, min.homo = 2, stacked = F, homo_homeo = "homo")
heatmap.fn(long.data, "Hek293T", gt = order, min = 1, max = 100, min.homo = 0, stacked = F, homo_homeo = "homeo")
heatmap.fn(long.data, "Hek293T", gt = order, min = 1, max = 100, min.homo = 0, stacked = F, homo_homeo = "homo")
heatmap.fn(long.data, "Hek293T", gt = order, min = 1, max = 250, min.homo = 2, stacked = F, homo_homeo = "homeo")
heatmap.fn(long.data, "Hek293T", gt = order, min = 1, max = 250, min.homo = 2, stacked = F, homo_homeo = "homo")
heatmap.fn(long.data, "Hek293T", gt = order, min = 1, max = 250, min.homo = 0, stacked = F, homo_homeo = "homeo")
heatmap.fn(long.data, "Hek293T", gt = order, min = 1, max = 250, min.homo = 0, stacked = F, homo_homeo = "homo")
heatmap.fn(long.data, "Hek293T", gt = order, min = 1, max = 25, min.homo = 3, stacked = F, homo_homeo = "homeo")
heatmap.fn(long.data, "Hek293T", gt = order, min = 1, max = 25, min.homo = 3, stacked = F, homo_homeo = "homo")
heatmap.fn(long.data, "Hek293T", gt = order, min = 1, max = 250, min.homo = 3, stacked = F, homo_homeo = "homeo")
heatmap.fn(long.data, "Hek293T", gt = order, min = 1, max = 250, min.homo = 3, stacked = F, homo_homeo = "homo")
heatmap.fn(long.data, "Hek293T", gt = order, min = 100, max = 250, min.homo = 3, stacked = F, homo_homeo = "homeo")
heatmap.fn(long.data, "Hek293T", gt = order, min = 100, max = 250, min.homo = 3, stacked = F, homo_homeo = "homo")

order <- c("WT", "PolQ-PolDel", "PolQ-HelDel", "PolQ-PolDHelD")

heatmap.fn(long.data, "DLD1", gt = order, min = 1, max = 25, min.homo = 2, stacked = F, homo_homeo = "homeo")
heatmap.fn(long.data, "DLD1", gt = order, min = 1, max = 25, min.homo = 2, stacked = F, homo_homeo = "homo")
heatmap.fn(long.data, "DLD1", gt = order, min = 1, max = 25, min.homo = 0, stacked = F, homo_homeo = "homeo")
heatmap.fn(long.data, "DLD1", gt = order, min = 1, max = 25, min.homo = 0, stacked = F, homo_homeo = "homo")
heatmap.fn(long.data, "DLD1", gt = order, min = 1, max = 100, min.homo = 2, stacked = F, homo_homeo = "homeo")
heatmap.fn(long.data, "DLD1", gt = order, min = 1, max = 100, min.homo = 2, stacked = F, homo_homeo = "homo")
heatmap.fn(long.data, "DLD1", gt = order, min = 1, max = 100, min.homo = 0, stacked = F, homo_homeo = "homeo")
heatmap.fn(long.data, "DLD1", gt = order, min = 1, max = 100, min.homo = 0, stacked = F, homo_homeo = "homo")
heatmap.fn(long.data, "DLD1", gt = order, min = 1, max = 250, min.homo = 2, stacked = F, homo_homeo = "homeo")
heatmap.fn(long.data, "DLD1", gt = order, min = 1, max = 250, min.homo = 2, stacked = F, homo_homeo = "homo")
heatmap.fn(long.data, "DLD1", gt = order, min = 1, max = 250, min.homo = 0, stacked = F, homo_homeo = "homeo")
heatmap.fn(long.data, "DLD1", gt = order, min = 1, max = 250, min.homo = 0, stacked = F, homo_homeo = "homo")
heatmap.fn(long.data, "DLD1", gt = order, min = 1, max = 25, min.homo = 3, stacked = F, homo_homeo = "homeo")
heatmap.fn(long.data, "DLD1", gt = order, min = 1, max = 25, min.homo = 3, stacked = F, homo_homeo = "homo")
heatmap.fn(long.data, "DLD1", gt = order, min = 1, max = 250, min.homo = 3, stacked = F, homo_homeo = "homeo")
heatmap.fn(long.data, "DLD1", gt = order, min = 1, max = 250, min.homo = 3, stacked = F, homo_homeo = "homo")
heatmap.fn(long.data, "DLD1", gt = order, min = 100, max = 250, min.homo = 3, stacked = F, homo_homeo = "homeo")
heatmap.fn(long.data, "DLD1", gt = order, min = 100, max = 250, min.homo = 3, stacked = F, homo_homeo = "homo")

unique(long.data[cell_line == "U2OS"]$genotype)

order <- c("DMSO", "NVB-100uM", "ART558-25uM",  "RP6685-10uM", "RP9913-10uM")


heatmap.fn(long.data, "U2OS", gt = order, min = 1, max = 25, min.homo = 2, stacked = F, homo_homeo = "homeo")
heatmap.fn(long.data, "U2OS", gt = order, min = 1, max = 25, min.homo = 2, stacked = F, homo_homeo = "homo")
heatmap.fn(long.data, "U2OS", gt = order, min = 1, max = 25, min.homo = 0, stacked = F, homo_homeo = "homeo")
heatmap.fn(long.data, "U2OS", gt = order, min = 1, max = 25, min.homo = 0, stacked = F, homo_homeo = "homo")
heatmap.fn(long.data, "U2OS", gt = order, min = 1, max = 100, min.homo = 2, stacked = F, homo_homeo = "homeo")
heatmap.fn(long.data, "U2OS", gt = order, min = 1, max = 100, min.homo = 2, stacked = F, homo_homeo = "homo")
heatmap.fn(long.data, "U2OS", gt = order, min = 1, max = 100, min.homo = 0, stacked = F, homo_homeo = "homeo")
heatmap.fn(long.data, "U2OS", gt = order, min = 1, max = 100, min.homo = 0, stacked = F, homo_homeo = "homo")
heatmap.fn(long.data, "U2OS", gt = order, min = 1, max = 250, min.homo = 2, stacked = F, homo_homeo = "homeo")
heatmap.fn(long.data, "U2OS", gt = order, min = 1, max = 250, min.homo = 2, stacked = F, homo_homeo = "homo")
heatmap.fn(long.data, "U2OS", gt = order, min = 1, max = 250, min.homo = 0, stacked = F, homo_homeo = "homeo")
heatmap.fn(long.data, "U2OS", gt = order, min = 1, max = 250, min.homo = 0, stacked = F, homo_homeo = "homo")
heatmap.fn(long.data, "U2OS", gt = order, min = 1, max = 25, min.homo = 3, stacked = F, homo_homeo = "homeo")
heatmap.fn(long.data, "U2OS", gt = order, min = 1, max = 25, min.homo = 3, stacked = F, homo_homeo = "homo")
heatmap.fn(long.data, "U2OS", gt = order, min = 1, max = 250, min.homo = 3, stacked = F, homo_homeo = "homeo")
heatmap.fn(long.data, "U2OS", gt = order, min = 1, max = 250, min.homo = 3, stacked = F, homo_homeo = "homo")
heatmap.fn(long.data, "U2OS", gt = order, min = 100, max = 250, min.homo = 3, stacked = F, homo_homeo = "homeo")
heatmap.fn(long.data, "U2OS", gt = order, min = 100, max = 250, min.homo = 3, stacked = F, homo_homeo = "homo")




## heatmap.fn(long.data, "U2OS",  min = 1, max = 250, min.homo = 2, stacked = F, homo_homeo = "homeo")
## heatmap.fn(long.data, "U2OS",  min = 1, max = , min.homo = 2, stacked = F, homo_homeo = "homo")
## heatmap.fn(long.data, "U2OS",  min = 1, max = , min.homo = 0, stacked = F, homo_homeo = "homeo")
## heatmap.fn(long.data, "U2OS",  min = 1, max = , min.homo = 0, stacked = F, homo_homeo = "homo")


## heatmap.fn(long.data, "U2OS",  min = 25, max = 100, min.homo = 2, stacked = F, homo_homeo = "homeo")
## heatmap.fn(long.data, "U2OS",  min = 25, max = 100, min.homo = 2, stacked = F, homo_homeo = "homo")
## heatmap.fn(long.data, "U2OS",  min = 25, max = 100, min.homo = 0, stacked = F, homo_homeo = "homeo")
## heatmap.fn(long.data, "U2OS",  min = 25, max = 100, min.homo = 0, stacked = F, homo_homeo = "homo")
## heatmap.fn(long.data, "U2OS",  min = 50, max = 100, min.homo = 2, stacked = F, homo_homeo = "homeo")
## heatmap.fn(long.data, "U2OS",  min = 50, max = 100, min.homo = 2, stacked = F, homo_homeo = "homo")
## heatmap.fn(long.data, "U2OS",   min = 50, max = 100, min.homo = 0, stacked = F, homo_homeo = "homeo")
## heatmap.fn(long.data, "U2OS",  min = 50, max = 100, min.homo = 0, stacked = F, homo_homeo = "homo")





Z111 = readRDS("~/projects/polq/jr_stash/all_data/PE250_Hek_data.rds")



#' jrafailov Tuesday, Jan 07, 2025 03:56:01 PM
#' stacked bar plot attempt


#' pull in one sample for which we have all the deletions enumerated

get_legend <- function(plot, legend = NULL) {
  gt <- ggplotGrob(plot)
  pattern <- "guide-box"
  if (!is.null(legend)) {
    pattern <- paste0(pattern, "-", legend)
  }
  indices <- grep(pattern, gt$layout$name)
  not_empty <- !vapply(
                  gt$grobs[indices],
                  inherits, what = "zeroGrob",
                  FUN.VALUE = logical(1)
                )
  indices <- indices[not_empty]
  if (length(indices) > 0) {
    return(gt$grobs[[indices[1]]])
  }
  return(NULL)
}


short <- c("Project_13108", # HEK293T 245bp IGO
           "PE250_Hek", # HEK293T 245bp GENEWIZ!!!!!
           "PE150_extra", # HEK293T 245bp GENEWIZ !!!
           "PE150_DLD1", # DLD1 245bp IGO
           "Project_14818") #U2OS 245bp IGO

long <- c("Project_14789", #Hek293T 503bp
          "Project_14783", #mixed Hek293T/DLD1 503bp
          "Project_14629", #DLD1 503bp
          "Project_15074")


mclapply(short, function(x){
  print(x)
  readRDS(paste0(paste0("~/projects/polq/jr_stash/all_data2/", x, "_data.rds"))) -> a
  a$project <- x
  return(a)
}, mc.cores = 4) %>%
  rbindlist(fill = T) %>%
  mutate(
    genotype = gsub("_C", "", genotype)
  ) %>%
  mutate(
    genotype = case_when(
      is.na(genotype) ~ paste0(treatment, ifelse(is.na(volume), "", paste0("-", volume))),
      T ~ genotype
    )
  )-> short.data

small.values <- readRDS("~/projects/polq/jr_stash/final_prediction_matrix.SMALL.HBONDS.rds") %>% as.data.table() %>%
  mutate(position = position - 1) %>%
  transmute(position, deletion_size, homology, homeology = homeology0.8) %>%
  as.data.table()

short.data <- short.data %>%
  merge(small.values, by.x = c("first_del", "deletion"), by.y = c("position", "deletion_size"), all.x = T) %>%
  as.data.table()

mclapply(long, function(x){
  print(x)
  readRDS(paste0(paste0("~/projects/polq/jr_stash/all_data2/", x, "_data.rds"))) -> a
  a$project <- x
  return(a)
}, mc.cores = 4) %>% rbindlist(fill = T) %>%
  mutate(
    genotype = gsub("_C", "", genotype)
  ) %>%
  mutate(
    genotype = case_when(
      is.na(genotype) ~ paste0(treatment, ifelse(is.na(volume), "", paste0("-", volume))),
      T ~ genotype
    )
  )-> long.data

long.values <- readRDS("~/projects/polq/jr_stash/final_prediction_matrix.LARGE.HBONDS.rds") %>% as.data.table() %>%
  mutate(position = position - 1) %>%
  transmute(position, deletion_size, homology, homeology = homeology0.8) %>%
  as.data.table()

long.data <- long.data %>%
  merge(long.values, by.x = c("first_del", "deletion"), by.y = c("position", "deletion_size"), all.x = T) %>%
  as.data.table()
  # data save

short.data.expt <- short.data %>%
  rowwise() %>%
  mutate(insertion_coordinates = unlist(insertion_coordinates) %>% paste0(collapse = ","),
         substitution_positions = unlist(substitution_positions) %>% paste0(collapse = ","),
         deletion_coordinates = unlist(deletion_coordinates) %>% paste0(collapse = ",")) %>%
  group_by(
    across(-c(starts_with("total_"), reads, sample))
  ) %>%
  summarize(reads = sum(reads),
            across(starts_with("total_"), ~ sum(unique(.)), .names = "{.col}"),
            samples = paste0(unique(sample), collapse = ", ")) %>%
  as.data.table()

long.data.expt <- long.data %>%
  rowwise() %>%
  mutate(insertion_coordinates = unlist(insertion_coordinates) %>% paste0(collapse = ","),
         substitution_positions = unlist(substitution_positions) %>% paste0(collapse = ","),
         deletion_coordinates = unlist(deletion_coordinates) %>% paste0(collapse = ",")) %>%
  group_by(
    across(-c(starts_with("total_"), reads, sample))
  ) %>%
  summarize(reads = sum(reads),
            across(starts_with("total_"), ~ sum(unique(.)), .names = "{.col}"),
            samples = paste0(unique(sample), collapse = ", ")) %>%
  as.data.table()

write_tsv(short.data.expt, "~/projects/polq/jr_stash/SHORT_DATA_experiment_collapsed.tsv")
write_tsv(long.data.expt, "~/projects/polq/jr_stash/LONG_DATA_experiment_collapsed.tsv")


  #pick a sample
  # lets do the DLD1 sample WT sample


col_fun_prop = readRDS("~/projects/polq/jr_stash/col_fun_prop/large.rds")

col_fun_prop = c(rep("white", 2),col_fun_prop(3:11))

col_fun_prop = colorRamp2(0:10,c("white", col_fun_prop(2:11)), space = "RGB")

col.groupings = wes_palette("Darjeeling1") #rev(khroma::color("high contrast")(5))

col.groupings = pal_uchicago("light")(6)
names(col.groupings) <- c("1", "2", "3-10", "11-25", "26-50", ">50")

col.groupings = pal_uchicago("light")(6)[c(1,6,3)]
names(col.groupings) <- c("1-4bp", "5-25bp", ">25bp")



col_fun_prop = readRDS("~/projects/polq/jr_stash/col_fun_prop/all.rds")
col_fun_prop = colorRamp2(c(0,1,2:11), c(rep("white",2),  col_fun_prop(2:11)), space = "RGB")


col_fun_prop = readRDS("~/projects/polq/jr_stash/col_fun_prop/ALL_round2.rds")

col_fun_prop = readRDS("~/projects/polq/jr_stash/col_fun_prop/custom.rds")
col_fun_prop = colorRamp2(c(0,1,2:11), c(rep("white",2),  col_fun_prop(2:11)), space = "RGB")


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
          #size ==  1 & class == "del" ~ "1",
          #size ==  2 & class == "del"~ "2",
          size <= 4 & size >= 1 & class == "del" ~ "1-4bp",
          size <= 25 & size >= 5 & class == "del" ~ "5-25bp",
          size > 25 & class == "del" ~ ">25bp",
          #size >= 51 & class == "del" ~ ">50",
          T ~ NA
        ), levels = c("1-4bp", "5-25bp", ">25bp")),
        ## grouping = factor(case_when(
        ##   size ==  1 & class == "del" ~ "1",
        ##   size ==  2 & class == "del"~ "2",
        ##   size <= 10 & size > 2 & class == "del" ~ "3-10",
        ##   size <= 25 & size > 10 & class == "del" ~ "11-25",
        ##   size <= 50 & size > 25 & class == "del" ~ "26-50",
        ##   size >= 51 & class == "del" ~ ">50",
        ##   T ~ NA
        ## ), levels = c("1", "2", "3-10", "11-25", "26-50", ">50")),
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
      #browser()
      unique.del.groups <- as.numeric(unique(sort(data.all$grouping)))
      del.tranche.groupings = list(#unique.del.groups[1],
                                   #unique.del.groups[length(unique.del.groups) - 1],
                                   unique.del.groups[length(unique.del.groups)])
      ## del.tranche.groupings = list(unique.del.groups[1:2],
      ##                              unique.del.groups[3:4],
      ##                              unique.del.groups[length(unique.del.groups) - 1],
      ##                              unique.del.groups[length(unique.del.groups)])
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
      # lapply(del.tranche.groupings, function(del.grouping){
      #   tornado.graph.fn(data, plot.ins = plot.ins,
      #                    plot.mut = plot.mut,
      #                    x.spacing = x.spacing, adjust.size = adjust.size,
      #                    gts = gts, return.del = T, del.grouping = del.grouping) -> group.call
      #   return(group.call[[1]])
      # }) -> all.del.plots
      # hom.min.groupings[[2]]
      #  tornado.graph.fn(data, plot.ins = plot.ins,
      #                    plot.mut = plot.mut,
      #                    x.spacing = x.spacing, adjust.size = adjust.size,
      #                    gts = gts, return.del = T, del.grouping = del.tranche.groupings[[1]],
      #                    hom.min = hom.min.groupings[[1]]) -> group.call1

      #        tornado.graph.fn(data, plot.ins = plot.ins,
      #                    plot.mut = plot.mut,
      #                    x.spacing = x.spacing, adjust.size = adjust.size,
      #                    gts = gts, return.del = T, del.grouping = del.tranche.groupings[[2]],
      #                    hom.min = hom.min.groupings[[2]]) -> group.call2
      # all.del.plots = list(group.call1[[1]], group.call2[[1]])
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
    ## if(return.mut == T){
    ##   x.max <- 20
    ## } else if (return.del == T & (1 %in% del.grouping)){
    ##   x.max <- 20
    ## } else{
    ##   x.max <- max(data.all.gt$first_del + data.all.gt$size, na.rm = T) - adjust.size
    ## }
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
      ## theme2 <- theme(legend.position = "none",
      ##                 plot.title = element_blank(),
      ##                 axis.text.x = element_text(hjust = 0.5,
      ##                                            vjust = 0.5,
      ##                                            size = 2.25, margin = margin(t = 0.4), angle = 270),
      ##                 plot.subtitle = element_text(hjust = 0.5, size = 3, margin = margin(t = 0, b = 0.2)),
      ##                 axis.title.y = element_blank(),
      ##                 axis.text.y = element_text(size = 2, vjust = 0.5, margin = margin(r = 0.4)),
      ##                 axis.line = element_line(size = 0.1),
      ##                 axis.ticks = element_line(size = 0.1),     Major ticks
      ##                 axis.ticks.length = unit(0.02, "cm"),
      ##                 axis.title.x = element_blank(),
      ##                 plot.margin = unit(c(0.2, 2, 0.2, 0.2), "mm"))
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
      #  fill = size
      ), fill = "#808080") +
      # scale_fill_scico(palette = "turku",
      #                  limits = c(1, max.del.size),
      #                  breaks = c(1, seq(25, max.del.size, 25)),
      #                  begin = 0.0,
      #                  end = 1,
      #                  direction = -1,
      #                  alpha = 1,
      #                  aesthetics = c("fill"),
      #                  name = "size") +
      # geom_rect(aes(
      #   fill = size,
      #   ymin = x.max,
      #   ymax = x.max + tranche.marker,
      #   xmin = ymin,
      #   xmax = ymax
      # ), size = 0) +
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

data <- short.data
cell <- "DLD1"
gts = NULL
suffix = "short"
plot.ins = T
plot.mut = T
zoom.min.size = 12

undebug(tornado_plot)


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









