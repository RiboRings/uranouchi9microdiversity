n_top_taxa <- 15

time_series_df <- filtered_df %>%
  slice_max(AbundMax, n = n_top_taxa)

time_series_df1 <- time_series_df %>%
  select(genome,
         starts_with("nucdiv_")) %>%
  pivot_longer(starts_with("nucdiv"),
               values_to = "NucDiv",
               names_to = "Sample") %>%
  mutate(Sample = as.numeric(gsub("nucdiv_", "", Sample)))

time_series_df2 <- time_series_df %>%
  select(genome,
         starts_with("RPKM")) %>%
  pivot_longer(starts_with("RPKM"),
               values_to = "RPKM",
               names_to = "Sample") %>%
  mutate(Sample = as.numeric(gsub("RPKM_", "", Sample)))

time_series_df <- time_series_df1 %>%
  full_join(time_series_df2) %>%
  mutate(genome = gsub("metabat2bin_", "", genome))

p <- ggplot(time_series_df, aes(x = Sample)) +
  geom_line(aes(y = RPKM), colour = "Black") +
  geom_line(aes(y = NucDiv * 1000), colour = "Dark Gray") +
  scale_x_continuous(limits = c(1, 9),
                     breaks = c(1, 3, 5, 7, 9)) +
  scale_y_continuous(name = "RPKM",
                     sec.axis = sec_axis(~./1000,
                                         name = "Nucleotide Diversity",
                                         )) +
  facet_grid(~genome) +
  theme_bw() +
  theme(axis.title.y = element_text(colour = "Black"),
        axis.title.y.right = element_text(colour = "Dark Gray"))

ggsave("multifacets.pdf",
       plot = p,
       device = "pdf",
       path = "results",
       width = n_top_taxa,
       height = 3)

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
