selected_genomes <- paste0(top_df$genome, ".fna")

time_series_df1 <- top_df %>%
  select(genome,
         Tax,
         starts_with("nucdiv_")) %>%
  pivot_longer(starts_with("nucdiv"),
               values_to = "NucDiv",
               names_to = "Sample") %>%
  mutate(Sample = as.numeric(gsub("nucdiv_", "", Sample)))

time_series_df2 <- top_df %>%
  select(genome,
         starts_with("RPKM")) %>%
  pivot_longer(starts_with("RPKM"),
               values_to = "RPKM",
               names_to = "Sample") %>%
  mutate(Sample = as.numeric(gsub("RPKM_", "", Sample)))

time_series_df <- time_series_df1 %>%
  full_join(time_series_df2) %>%
  mutate(genome = gsub("metabat2bin_", "", genome))

nucdiv_range_df <- time_series_df %>%
  group_by(genome) %>%
  summarise(NucDivRange = max(NucDiv, na.rm = TRUE) - min(NucDiv, na.rm = TRUE))

time_series_df <- time_series_df %>%
  left_join(nucdiv_range_df) %>%
  arrange(desc(NucDivRange))

p <- ggplot(time_series_df, aes(x = Sample)) +
  geom_line(aes(y = RPKM), colour = "Black") +
  geom_line(aes(y = NucDiv * 1000), colour = "Dark Gray") +
  scale_x_continuous(limits = c(1, 9),
                     breaks = c(1, 3, 5, 7, 9)) +
  scale_y_continuous(name = "RPKM",
                     sec.axis = sec_axis(~./1000,
                                         name = "Nucleotide Diversity",
                                         )) +
  facet_wrap(~ NucDivRange + Tax,
             nrow = 4) +
  theme_bw() +
  theme(axis.title.y = element_text(colour = "Black"),
        axis.title.y.right = element_text(colour = "Dark Gray"),
        panel.grid = element_blank())

ggsave("multifacets.pdf",
       plot = p,
       device = "pdf",
       path = "results",
       width = 40,
       height = 10)
