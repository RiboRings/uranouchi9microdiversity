ubiquity_mat <- micro_df %>%
  select(genome, breadth, sample) %>%
  pivot_wider(names_from = sample,
              values_from = breadth,
              values_fill = 0) %>%
  column_to_rownames(var = "genome") %>%
  as.matrix()

ubiquity <- sapply(iter(ubiquity_mat, by = "row"), function(x) sum(x > 0.50))
ubiquity_df <- data.frame("Ubiquity" = as.factor(ubiquity))

p <- ggplot(ubiquity_df, aes(x = Ubiquity)) +
  geom_bar() +
  scale_y_continuous(limits = c(0, 315),
                     breaks = seq(0, 300, by = 50)) +
  labs(y = "Number of MAGs") +
  theme_bw() +
  theme(panel.grid = element_blank())

ggsave("ubiquity.pdf",
       plot = p,
       device = "pdf",
       path = "results",
       width = 10,
       height = 7)
