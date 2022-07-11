func_annot <- read_csv2("data/func_annot.csv")

func_annot_df <- func_annot %>%
  group_by(product) %>%
  summarise(gene,
            pNpS,
            Count = n()) %>%
  arrange(desc(Count)) %>%
  filter(!is.na(pNpS),
         Count > 4)
  
p <- ggplot(func_annot_df, aes(y = pNpS, x = product)) +
  geom_boxplot() +
  labs(x = "Functional Group") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggsave("func_annot.pdf",
       plot = p,
       width = 15,
       height = 7,
       path = "results")
