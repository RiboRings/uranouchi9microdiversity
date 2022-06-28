n_edge_taxa <- 10

selected_df1 <- filtered_df %>%
  arrange(desc(NucDivMean), desc(AbundMax)) %>%
  slice_head(n = n_edge_taxa)

selected_df2 <- filtered_df %>%
  arrange(NucDivMean, desc(AbundMax)) %>%
  slice_head(n = n_edge_taxa)

selected_df3 <- filtered_df %>%
  arrange(desc(AbundMax)) %>%
  slice_head(n = n_edge_taxa)

selected_df <- rbind(selected_df1,
                     selected_df2,
                     selected_df3) %>%
  distinct()

selected_genomes <- paste0(selected_df$genome, ".fna")

no_cores <- detectCores() - 1  
cl <- makeCluster(no_cores, type = "FORK")  
registerDoParallel(cl)

files <- foreach(i = 1:9) %dopar% {

    read_csv2(paste0("data/snvs", i, ".csv")) %>%
    na.omit %>%
    filter(genome %in% selected_genomes) %>%
    mutate(sample = i,
           genome = gsub(".fna", "", gsub("metabat2bin_", "", genome)),
           position = as.factor(position))
    
}

stopCluster(cl)

snvs_df <- files %>%
  merge_all()

snvs_stats <- snvs_df %>%
  group_by(genome, position) %>%
  summarise(Count = n())

snvs_df <- snvs_df %>%
  left_join(snvs_stats)

write_csv2(snvs_df, "snvs.csv")

### other stuff

file1 <- read_csv2("data/mags.csv") %>%
  select(-N50, -length)
file2 <- read_csv("genomeInformation.csv") %>%
  transmute(Bin = gsub(".fna", "", genome),
            length,
            N50)
file3 <- file1 %>%
  left_join(file2)

write_csv2(df, "data/mags.csv")

test <- micro_df
test2 <- test %>% mutate(sample = paste0("scaffolds_", sample)) %>%
  filter(genome %in% unique(df$Bin))
df <- df[order(df$Bin), ]
test2 <- test2[order(test2$genome), ]
test2 <- pivot_wider(test2, names_from = sample, values_from = true_scaffolds)
test2 <- test2 %>% column_to_rownames(var = "genome") %>%
  as.matrix()

lol <- rowMeans(test2, na.rm = TRUE) %>% as.vector() %>% as.numeric()
df$scaffolds <- lol
