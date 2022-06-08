filtered_df$DiSiperMbpMean <- filtered_df %>%
  select(starts_with("DiSiperMbp")) %>%
  as.matrix() %>%
  rowMeans(na.rm = TRUE)

filtered_df$NonsynonimousFractionMean <- filtered_df %>%
  select(starts_with("NonsynonimousFraction")) %>%
  as.matrix() %>%
  rowMeans(na.rm = TRUE)

p <- ggplot(filtered_df, aes(x = DiSiperMbpMean,
           y = NonsynonimousFractionMean,
           colour = AbundMax)) +
  geom_point(size = 1) +
  labs(x = "Mean Divergent Sites per Mbp",
       y = "Mean Nonsynonymous Fraction",
       colour = "Max RPKM") +
  theme_classic()

ggsave("microdiversity.pdf",
       plot = p,
       device = "pdf",
       path = "results",
       width = 10,
       height = 7)
