snv_files <- list.files(pattern = 'snvs.+.csv',
                              path = "data")

snv_file_list <- mclapply(paste0("data/", snv_files),
                          sample_loader,
                          mc.cores = detectCores() - 1)

snv_df <- merge_all(snv_file_list) %>%
  na.omit() %>%
  mutate(sample = as.numeric(gsub("snvs", "", sample)),
         genome = gsub(".fna", "", gsub("metabat2bin_", "", genome)),
         position = as.factor(position)) %>%
  filter(genome %in% unique(top_df$genome))

write.csv2(snv_df, "data/snvs.csv")
