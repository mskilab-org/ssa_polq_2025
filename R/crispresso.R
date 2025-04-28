#' sbrylka Friday, Sep 22, 2023 02:33:56 PM
#' pipeline to prepare data/graphs for homology project

library(Flow)
library(dplyr)
library(tidyr)
library(skitools)
library(ggpubr)
library(pbmcapply)

# Preparing entities for the Flow job
# choose only one path depending on which crispresso to run
# all samples have already been crispresso'd, look in "/gpfs/commons/home/sbrylka/Projects/homeology"  

#path = "/gpfs/commons/home/sbrylka/data/MSK_LDM/concat_fastq"
#path = "/gpfs/commons/home/sbrylka/data/MSK_LDM/Project_15074/all_fastq"
#path = "/gpfs/commons/home/sbrylka/data/MSK_LDM/PE250/240bp_PCR_PE250_NGS_Genewiz" # PE250 Hek
#path = "/gpfs/commons/home/sbrylka/data/MSK_LDM/Project_13108/all_fastq" # PE150 Hek293
#path = "/gpfs/commons/home/sbrylka/data/MSK_LDM/Project_14818/all_fastq" # PE150 PolQ
#path = "/gpfs/commons/home/sbrylka/data/MSK_LDM/PE150_DLD1/all_fastq" # PE150 DLD1
#path = "/gpfs/commons/home/sbrylka/data/MSK_LDM/PE150_extra/PE150_genewiz_extra" # PE150 DLD1 extra
#path = "/gpfs/commons/home/sbrylka/data/MSK_LDM/Project_14783"
#path = "/gpfs/commons/home/sbrylka/data/MSK_LDM/Project_14789"
path = "/gpfs/commons/home/sbrylka/data/MSK_LDM/Project_14629"
fastq.list <- list.files(path = path)
fastq.pairs <- data.table(
  FASTQPathFile = list.files(path, full.names = T, pattern ="fastq"),
  name = list.files(path, pattern ="fastq")) %>%
  rowwise() %>%
  mutate(
    name = str_split(name, '_'),
    sample_name = paste(name[1], name[2], name[3], name[4], name[5], sep = "_"),
    read = name[length(name)-1]) %>%
  select(-name) %>%
pivot_wider(
    id_cols = "sample_name",
    names_from = read,
    values_from = FASTQPathFile
  ) %>%
  data.table()
fqdt <- fastq.pairs %>% transmute(name = sample_name, r1 = R1, r2 = R2)


# for PE150 samples
amp <- "TTAATGTGGCTCTGGTTCTGGGTacttttatctgtcccctccaccccacagtggggccactagggacaggattggtgacagaaaagccccatccttaggcctcctccttcctagtctcctgatattgggtctaacccccacctcctgttaggcagattccttatctggtgacacacccccatttcctggagccatctctctccttgccagaacctctaaggtttGCTTACGATGGAGCCAGAGAG"

amp <- "TTACCTCTCTAGTCTGTGCTAGCTCTTCCAGCCCCCTGTCATGGCATCTTCCAGGGGTCCGAGAGCTCAGCTAGTCTTCTTCCTCCAACCCGGGCCCCTATGTCCACTTCAGGACAGCATGTTTGCTGCCTCCAGGGATCCTGTGTCCCCGAGCTGGGACCACCTTATATTCCCAGGGCCGGTTAATGTGGCTCTGGTTCTGGGTACTTTTATCTGTCCCCTCCACCCCACAGTGGGGCCACTAGGGACAGGATTGGTGACAGAAAAGCCCCATCCTTAGGCCTCCTCCTTCCTAGTCTCCTGATATTGGGTCTAACCCCCACCTCCTGTTAGGCAGATTCCTTATCTGGTGACACACCCCCATTTCCTGGAGCCATCTCTCTCCTTGCCAGAACCTCTAAGGTTTGCTTACGATGGAGCCAGAGAGGATCCTGGGAGGGAGAGCTTGGCAGGGGGTGGGAGGGAAGGGGGGGATGCGTGACCTGCCCGGTTCTCAGTGGCCA"

crisp_ent <- data.table(
  read1 = fqdt$r1,
  read2 = fqdt$r2,
  sample = fqdt$name,
  amplicon_seq = amp,
  flags = "-fg GGGGCCACTAGGGACAGGAT -fh 80 --quantification_window_size 1 --write_detailed_allele_table") %>%
  #flags = "--bam_output -fg GGGGCCACTAGGGACAGGAT -fh 80 --min_single_bp_quality 10 --min_average_read_quality 30 -amas 80 --default_min_aln_score 80 --max_paired_end_reads_overlap 503 --quantification_window_size 1 --write_detailed_allele_table") %>%
  setkey(sample)

Flow::Task("/gpfs/commons/home/sbrylka/Projects/tasks/crispresso.task")

crisp.job <- Flow::Job(
  "/gpfs/commons/home/sbrylka/Projects/tasks/crispresso.task",
  crisp_ent,
  "/gpfs/commons/home/sbrylka/Projects/homeology/Project_14629/crispresso_Project_14629_job_cut_site",
  time = '5-00',
  mem = 32, 
  cores = 8, 
  update_cores = 1
)

Flow::purge(crisp.job)

Flow::purge(crisp.job['running and some outputs present'])

Flow::purge(crisp.job['failed and some outputs present'])

Flow::srun(crisp.job)

Flow::update(crisp.job)


# Prepare Files for analysis
dt <- Flow::outputs(crisp.job['completed']) %>%
  mutate(unzip_cmd = paste0("unzip ",crisp.out,"/Alleles_frequency_table.zip -d ",crisp.out, sep = ""))

unzip <- lapply(1:nrow(dt),
                function(x){
                  system(dt$unzip_cmd[x])
                }
)

saveRDS(dt, "/gpfs/commons/home/sbrylka/Projects/homeology/Project_14629/dt.rds")
