df$MeanDiSiperMbp <- df %>%
  select(starts_with("DiSiperMbp")) %>%
  as.matrix() %>%
  rowMeans(na.rm = TRUE)

df$MeanNonsynonimousFraction <- df %>%
  select(starts_with("NonsynonimousFraction")) %>%
  as.matrix() %>%
  rowMeans(na.rm = TRUE)

disi_array <- sapply(iter(filtered_df, by = "row"),
                     max_mapper, pattern = "DiSiperMbp")
nsf_array <- sapply(iter(filtered_df, by = "row"),
                    max_mapper, pattern = "NonsynonimousFraction")

tax_df <- filtered_df %>%
  transmute(Tax = paste(Order, Family, gsub("metabat2bin_", "", genome), sep = ";"))

microdiversity_df <- data.frame("DiSiperMbpAtAbundMax" = disi_array,
                                "NonsynonimousFractionAtAbundMax" = nsf_array,
                                "AbundMax" = filtered_df$AbundMax,
                                "Tax" = tax_df$Tax) %>%
  mutate(Tax = ifelse(NonsynonimousFractionAtAbundMax > 0.5, Tax, ""))

p <- ggplot(microdiversity_df, aes(x = DiSiperMbpAtAbundMax,
            y = NonsynonimousFractionAtAbundMax,
            colour = AbundMax)) +
  geom_point(size = 1) +
  geom_text_repel(aes(label = Tax),
                  size = 3,
                  colour = "red") +
  labs(x = "Divergent Sites per Mbp at Max RPKM",
       y = "Nonsynonymous Fraction at Max RPKM",
       colour = "Max RPKM") +
  theme_classic()

ggsave("microdiversity.pdf",
       plot = p,
       device = "pdf",
       path = "results",
       width = 10,
       height = 7)
