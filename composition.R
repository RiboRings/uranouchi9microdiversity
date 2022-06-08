n_top_taxa <- 50

heatmap_df <- filtered_df %>%
  slice_max(AbundMax, n = n_top_taxa)

heatmap_mat <- heatmap_df %>%
  mutate(genome = gsub("metabat2bin_", "", genome)) %>%
  select(genome, starts_with("RPKM")) %>%
  column_to_rownames(var = "genome") %>%
  as.matrix()

colnames(heatmap_mat) <- seq(1, dim(heatmap_mat)[[2]])

ha <- HeatmapAnnotation(ND = rowMeans(heatmap_df %>% select(starts_with("nucdiv"))),
                        which = "row")

Heatmap(heatmap_mat,
        heatmap_legend_param = list(at = seq(min(heatmap_mat), max(heatmap_mat), by = 10)),
        col = c("blue", "cyan", "green", "yellow", "red"),
        name = "RPKM",
        cluster_rows = TRUE,
        cluster_columns = FALSE,
        row_title = "MAG",
        column_title = "Sample",
        row_names_gp = gpar(fontsize = 5),
        left_annotation = ha)

