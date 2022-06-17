fuhrman_df <- read_csv2("data/snvs.csv")

fuhrman_df <- fuhrman_df %>%
  select(genome,
         scaffold,
         position,
         con_base,
         A, C, T, G,
         sample)

cur_genome <- "354"

consensus_base_df <- fuhrman_df %>%
  filter(genome == cur_genome) %>%
  select(scaffold,
         position,
         con_base,
         sample) %>%
  pivot_wider(names_from = sample,
              values_from = con_base,
              values_fill = NA,
              names_prefix = "sample_") %>%
  na.omit()

consensus_base_df <- consensus_base_df[ , order(names(consensus_base_df))]

consensus_base_mat <- consensus_base_df %>%
  mutate(name = paste(scaffold, position, sep = ";")) %>%
  column_to_rownames(var = "name") %>%
  select(starts_with("sample")) %>%
  as.matrix()

set.seed(123)
Heatmap(consensus_base_mat,
        name = "Consensus Base",
        show_row_names = FALSE,
        column_labels = gsub("sample_", "", colnames(consensus_base_mat)),
        column_title = "Sample",
        row_title = "SNV")

base_count_df <- fuhrman_df %>%
  filter(genome == cur_genome) %>%
  transmute(scaffold,
            position,
            FreqA = A / (A + C + T + G),
            FreqC = C / (A + C + T + G),
            FreqT = T / (A + C + T + G),
            FreqG = G / (A + C + T + G),
            sample) %>%
  pivot_wider(names_from = sample,
              values_from = c(starts_with("Freq"))) %>%
  na.omit()

base_count_df <- base_count_df[ , order(names(base_count_df))]

freq_df_list <- list()

for (cur_sample in 1:9) {
  
  sample_consensus_base_df <- consensus_base_df %>%
    select(scaffold,
           position,
           matches(paste0("sample_", cur_sample)))
  
  colnames(sample_consensus_base_df) <- c("scaffold", "position", "con_base")
  base_count_df$con_base <- sample_consensus_base_df$con_base
  
  tmp_df <- sapply(iter(base_count_df, by = "row"),
                   freq_mapper) %>%
    t() %>%
    as.data.frame()
  
    colnames(tmp_df) <- paste0("Freq", seq(1, 9))
  tmp_df$con_base <- sample_consensus_base_df$con_base
  
  freq_df_list[[cur_sample]] <- tmp_df

}

avg_freq_df <- sapply(freq_df_list, function(x) colMeans(x %>% select(-con_base))) %>%
  t() %>%
  as.data.frame() %>%
  melt(variable_name = "Time") %>%
  transmute(Time = gsub("Freq", "", Time),
            Frequency = value)

avg_freq_df$RefSample <- as.factor(rep(seq(1, 9), 9))

p <- ggplot(avg_freq_df, aes(x = Time,
                             y = Frequency,
                             colour = RefSample,
                             group = RefSample)) +
  geom_point() +
  geom_line() +
  scale_y_continuous(limits = c(0.25, 1)) +
  labs(y = "Average Consensus Base Frequency",
       colour = "Reference Sample") +
  theme_bw()

ggsave("mean_example.pdf",
       plot = p,
       height = 7,
       width = 7,
       path = "results",
       device = "pdf")
