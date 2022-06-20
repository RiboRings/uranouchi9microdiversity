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
  left_join(snvs_stats) %>%
  filter(Count >= 6,
         mutation_type %in% c("S", "N", "I"))

for (gen in unique(snvs_df$genome)) {
  
  genome_df <- snvs_df %>%
    filter(genome == gen) %>%
    transmute(position,
              con_freq,
              sample,
              mutation_type) %>%
    group_by(position) %>%
    summarise(sample,
              con_freq,
              mutation_type = ifelse(n_distinct(mutation_type) > 1, "M", mutation_type)) %>%
    group_by(position, sample) %>%
    summarise(con_freq = mean(con_freq),
              mutation_type) %>%
    distinct() %>%
    pivot_wider(names_from = sample,
                values_from = con_freq,
                values_fill = 1) %>%
    column_to_rownames(var = "position")
  
  genome_mat <-  genome_df[ , order(names(genome_df))] %>%
    select(-mutation_type) %>%
    as.matrix()
  
  set.seed(123)
  
  row_ha <- HeatmapAnnotation(Type = genome_df$mutation_type,
                              Position = as.numeric(rownames(genome_mat)),
                              which = "row")
  
  tax_df <- selected_df %>%
    filter(genome == paste0("metabat2bin_", gen)) %>%
    transmute(Tax = paste(Order, Family, Genus, Species, sep = ";"))
  
  rpkm_df <- selected_df %>%
    filter(genome == paste0("metabat2bin_", gen)) %>%
    select(starts_with("RPKM")) %>%
    melt() %>%
    mutate(variable = gsub("RPKM_", "", variable)) %>%
    filter(variable %in% colnames(genome_mat))
  
  snv_shift_vector <- 0
  
  for (i in seq(1, dim(genome_mat)[[2]] - 1)) {
    
    snv_shift_vector <- snv_shift_vector %>%
      append(round(snv_dist(genome_mat[ , i + 1], genome_mat[ , i]),
                   digits = 3))
    
  }
  
  col_ha <- HeatmapAnnotation(RPKM = anno_barplot(rpkm_df$value),
                              Distance = anno_points(snv_shift_vector),
                              which = "col")
  
  limits <- c(round(min(genome_mat, na.rm = TRUE), 1), round(max(genome_mat, na.rm = TRUE), 1))
  breaks <- seq(limits[1], limits[2], by = 0.1)
  
  pdf(paste0("results/", gen, "_snvs.pdf"),
      width = 10,
      height = 20)
  
  draw(Heatmap(genome_mat,
               heatmap_legend_param = list(at = breaks),
               col = c("darkred", "red", "orange", "yellow", "white"),
               name = "Base Consensus Frequency",
               cluster_columns = TRUE,
               column_dend_reorder = FALSE,
               cluster_rows = TRUE,
               row_title = "Position",
               column_title = tax_df$Tax,
               show_row_names = FALSE,
               right_annotation = row_ha,
               top_annotation = col_ha))
  
  dev.off()

}
