library(readr)
library(dplyr)
library(doParallel)
library(reshape)
library(ggplot2)

no_cores <- detectCores() - 1  
cl <- makeCluster(no_cores, type = "FORK")  
registerDoParallel(cl)

files <- foreach(i = 1:9) %dopar% {

    read_csv2(paste0("data/snvs", i, ".csv")) %>%
    na.omit %>%
#    filter(genome %in% c("metabat2bin_8.fna",
#                         "metabat2bin_35.fna",
#                         "metabat2bin_354.fna",
#                         "metabat2bin_800.fna",
#                         "metabat2bin_449.fna")) %>%
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
  left_join(snvs_stats) %>%
  filter(Count >= 5,
         mutation_type == "N")

for (gen in unique(snvs_df$genome)) {
  
  genome_df <- snvs_df %>%
    filter(genome == gen) %>%
    select(-scaffold, -allele_count) %>%
    mutate(sample = paste0("sample_", sample)) %>%
    pivot_wider(names_from = sample,
                values_from = con_freq,
                values_fill = 0) %>%
    column_to_rownames(var = "position") %>%
    select(starts_with("sample"))
  
  names(genome_df) <- gsub("sample_", "", names(genome_df))
  genome_mat <-  genome_df[ , order(names(genome_df))] %>%
    as.matrix()
  
  limits <- c(round(min(genome_mat), 1), round(max(genome_mat), 1))
  breaks <- seq(limits[1], limits[2], by = 0.1)
  
  Heatmap(genome_mat,
          heatmap_legend_param = list(at = breaks),
          col = c("darkblue", "blue", "white", "red", "darkred"),
          name = "Base Consensus Frequency",
          cluster_rows = TRUE,
          cluster_columns = FALSE,
          row_title = "Position",
          column_title = "Sample",
          row_names_gp = gpar(fontsize = 5))
  
#  genome_df <- snvs_df %>%
#    filter(genome == gen)
#  
#  colours <- c("darkblue", "blue", "white", "red", "darkred")
#  limits <- c(round(min(genome_df$con_freq), 1), round(max(genome_df$con_freq), 1))
#  breaks <- seq(limits[1], limits[2], by = 0.1)
#  
#  ggplot(genome_df, aes(x = as.character(sample),
#                        y = position,
#                        fill = con_freq)) +
#    geom_tile() +
#    scale_fill_gradientn(name = "Consensus Base Frequency", 
#                         breaks = breaks,
#                         limits = limits,
#                         colours = colours) +
#    theme_bw() +
#    theme(axis.text.y = element_text(size = 5, angle = 45, hjust = 1),
#          legend.title = element_text(size = 7)) +
#    labs(x = "Sample",
#         y = "Genomic Site",
#         subtitle = gen)
  
    ggsave(paste0(gen, "_snvs.pdf"),
           device = "pdf",
           path = "results",
           width = 10,
           height = 20)

}
