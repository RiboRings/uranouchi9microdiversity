stats_df <- data.frame("Completeness" = mean(df$Completeness),
                       "Contamination" = mean(df$Contamination),
                       "MeanCoverage" = mean(df$MeanCov),
                       "MinCoverage" = min(df$MinCov),
                       "MeanBreadth" = mean(rowMeans(ubiquity_mat)),
                       "DiSiperMbp" = mean(df$MeanDiSiperMbp),
                       "NonSynonymousFraction" = mean(df$MeanNonsynonimousFraction, na.rm = TRUE),
                       "NucDiv" = mean(rowMeans(nucdiv_mat)),
                       "MaxRPKM" = mean(df$AbundMax)) %>%
  t()

colnames(stats_df) <- "U9"

write.csv2(stats_df, "results/statistics.csv")
