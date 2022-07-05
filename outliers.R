library(readr)
library(dplyr)
library(reshape)

func_annot <- read_tsv("data/538_func_annot/PROKKA_07042022.tsv")
mapping <- read_tsv("data/mapping.stb", col_names = c("scaffold", "genome"))
gene_files <- list.files(pattern = "gene_info*",
                         path = "data")
gene_file_list <- lapply(paste0("data/", gene_files),
                         sample_loader,
                         method = read_tsv)

gene_df <- merge_all(gene_file_list) %>%
  left_join(mapping) %>%
  transmute(genome,
         scaffold,
         gene,
         sample = gsub("gene_info", "", gsub(".tsv", "", sample)),
         dNdS_substitutions,
         pNpS_variants,
         nucl_diversity,
         coverage,
         breadth,
         SNV_count,
         divergent_site_count,
         gene_length,
         start,
         end)

gene_df <- gene_df[!is.na(gene_df$pNpS_variants), ]
gene_stats <- gene_df %>%
  group_by(genome, gene) %>%
  summarise(MaxpNpS = max(pNpS_variants)) %>%
  slice_max(MaxpNpS, n = 10)
