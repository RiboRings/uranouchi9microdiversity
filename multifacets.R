n_edge_taxa <- 15

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

time_series_df1 <- selected_df %>%
  select(genome,
         starts_with("nucdiv_")) %>%
  pivot_longer(starts_with("nucdiv"),
               values_to = "NucDiv",
               names_to = "Sample") %>%
  mutate(Sample = as.numeric(gsub("nucdiv_", "", Sample)))

time_series_df2 <- selected_df %>%
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
  facet_wrap(~ NucDivRange + genome,
             nrow = 2) +
  theme_bw() +
  theme(axis.title.y = element_text(colour = "Black"),
        axis.title.y.right = element_text(colour = "Dark Gray"),
        panel.grid = element_blank())

ggsave("multifacets.pdf",
       plot = p,
       device = "pdf",
       path = "results",
       width = length(selected_genomes) / 2,
       height = 6)

#time_series_df <- time_series_df %>%
#  mutate(LogRPKM = log10(RPKM))
#
#ggplot(time_series_df, aes(x = Sample)) +
#  geom_line(aes(y = LogRPKM), colour = "Black") +
#  geom_line(aes(y = NucDiv * 50), colour = "Dark Gray") +
#  scale_x_continuous(limits = c(1, 9),
#                     breaks = c(1, 3, 5, 7, 9)) +
#  scale_y_continuous(name = "Log RPKM",
#                     sec.axis = sec_axis(~./50,
#                                         name = "Nucleotide Diversity",
#                     )) +
#  facet_grid(~genome) +
#  theme_bw() +
#  theme(axis.title.y = element_text(colour = "Black"),
#        axis.title.y.right = element_text(colour = "Dark Gray"))
  