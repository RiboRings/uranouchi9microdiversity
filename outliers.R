library(dplyr)
library(tidyr)
library(readr)

sample_loader <- function(sample_name, method = read_csv) {
  
  sample_df <- method(sample_name) %>%
    mutate(sample = gsub(".csv", "", gsub("data/", "", sample_name))) %>%
    left_join(mapping) %>%
    filter(genome == "metabat2bin_835.fna") %>%
    
    return(sample_df)
  
}

mapping <- read_tsv("data/mapping.stb", col_names = c("scaffold", "genome"))

gene_info_files <- list.files(path = "data", pattern = "gene_info*")

gene_info_list <- lapply(paste("data", gene_info_files, sep = "/"),
                         sample_loader,
                         method = read_tsv)

gene_info_df <- gene_info_list %>%
  merge_all()

gene_stats_df <- gene_info_df %>%
  group_by(gene) %>%
  summarise(pNpS = max(pNpS_variants)) %>%
  arrange(desc(pNpS))

func_annot <- read_csv2("test.gff")
write_csv2(left_join(gene_stats_df, func_annot), "data/new_func_annot.csv")
