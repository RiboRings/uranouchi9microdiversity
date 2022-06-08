library(plotly)

filtered_df$NucDivMean <- filtered_df %>%
  select(starts_with("nucdiv")) %>%
  as.matrix() %>%
  rowMeans(na.rm = TRUE)

plot_ly(data = filtered_df,
        x = ~Completeness,
        y = ~NucDivMean,
        z = ~Contamination,
        type = "scatter3d",
        mode = "markers",
        size = ~AbundMax) %>%
  layout(scene = list(yaxis = list(title = 'Mean ND')))
