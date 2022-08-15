n_top_taxa <- 50

filtered_df <- filtered_df %>%
  mutate(Tax = gsub("metabat2bin_", "", paste(Order, genome, sep = ";")))

top_df <- filtered_df %>%
  slice_max(AbundMax, n = n_top_taxa)

heatmap_mat <- top_df %>%
  select(Tax, starts_with("RPKM")) %>%
  column_to_rownames(var = "Tax") %>%
  as.matrix()

colnames(heatmap_mat) <- seq(1, dim(heatmap_mat)[[2]])

ha <- HeatmapAnnotation(ND = rowMeans(top_df %>% select(starts_with("nucdiv")), na.rm = TRUE),
                        which = "row")

pdf("results/composition.pdf",
    width = 15,
    height = 10)

draw(Heatmap(heatmap_mat,
             heatmap_legend_param = list(at = seq(min(heatmap_mat), max(heatmap_mat), by = 10)),
             col = c("blue", "cyan", "green", "yellow", "red"),
             name = "RPKM",
             cluster_rows = TRUE,
             cluster_columns = FALSE,
             row_title = "MAG",
             column_title = "Sample",
             row_names_gp = gpar(fontsize = 8),
             left_annotation = ha))

dev.off()
