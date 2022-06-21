mat <- df %>%
  select(starts_with("RPKM")) %>%
  as.matrix()

df$AbundMean <- rowMeans(mat)
df$AbundMedian <- rowMedians(mat)
df$AbundMin <- rowMins(mat)
df$AbundMax <- rowMaxs(mat)
df$AbundRange <- rowMaxs(mat) - rowMins(mat)
df$AbundSd <- rowSds(mat)

abund_max_idx <- c()
for (i in 1:nrow(rpkm_mat)) abund_max_idx <- append(abund_max_idx, which(rpkm_mat[i, ] == rowMaxs(rpkm_mat, na.rm = TRUE)[[i]])[[1]])
df$AbundMaxIdx <- abund_max_idx

df$AbundRatioMean <- rowMeans(rpkm_shift_mat)
df$AbundRatioMedian <- rowMedians(rpkm_shift_mat)
df$AbundRatioMin <- rowMins(rpkm_shift_mat)
df$AbundRatioMax <- rowMaxs(rpkm_shift_mat)
df$AbundRatioRange <- rowMaxs(rpkm_shift_mat) - rowMins(rpkm_shift_mat)
df$AbundRatioSd <- rowSds(rpkm_shift_mat)

stats <- lapply(names(df)[names(df) %>% startsWith("Abund")],
                plot_distribution)

(stats[[1]] /
    stats[[2]] /
    stats[[3]]) |
  (stats[[4]] /
     stats[[5]] /
     stats[[6]])

ratio_stats <- lapply(names(df)[names(df) %>% startsWith("AbundRatio")],
                           plot_distribution)

(ratio_stats[[1]] /
    ratio_stats[[2]] /
    ratio_stats[[3]]) |
  (ratio_stats[[4]] /
     ratio_stats[[5]] /
     ratio_stats[[6]])

min_coverage_threshold <- 5

filtered_df <- df %>%
  filter(MinCov > min_coverage_threshold) %>%
  arrange(desc(AbundMax))

mean_range_plot1 <- ggplot(df, aes(x = AbundMean, y = AbundRange, colour = AbundRatioRange)) +
  geom_point() +
  theme_classic()
mean_range_plot2 <- ggplot(filtered_df, aes(x = AbundMean, y = AbundRange, colour = AbundRatioRange)) +
  geom_point() +
  theme_classic()

mean_range_plot1 / mean_range_plot2

#ggplot(filtered_df, aes(y = AbundRange)) +
#  geom_boxjitter() +
#  theme_bw()
