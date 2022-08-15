library(readr)
library(dplyr)
library(tidyr)
library(purrr)
library(reshape)
library(MatrixGenerics)
library(ggplot2)
library(ggrepel)
library(patchwork)
library(ComplexHeatmap)
library(tibble)
library(doParallel)
library(stringr)

df <- read_csv2("data/mags.csv")

files <- list()
for (i in seq(1, 9)) {
  
  files[[i]] <- read_csv2(paste0("data/", i, ".csv")) %>%
    mutate(sample = i)
  
}

micro_df <- merge_all(files) %>%
  na.omit() %>%
  mutate(genome = gsub(".fna", "", genome)) %>%
  group_by(genome)

micro_stats <- micro_df %>%
  summarise(MeanCov = mean(coverage),
            MinCov = min(coverage),
            MaxCov = max(coverage))

df <- micro_stats %>%
  left_join(df, by = c("genome" = "Bin"))

rpkm_mat <- df %>%
  select(starts_with("RPKM")) %>%
  as.matrix()
rpkm_mat[rpkm_mat == 0] <- min(rpkm_mat[rpkm_mat > 0]) / 2
rpkm_shift_mat <- matrix(nrow = dim(rpkm_mat)[[1]], ncol = dim(rpkm_mat)[[2]] - 1)

nucdiv_mat <- df %>%
  select(starts_with("nucdiv"))

nucdiv_mat[is.na(nucdiv_mat)] <- min(nucdiv_mat[!is.na(nucdiv_mat)]) / 2
nucdiv_shift_mat <- matrix(nrow = dim(nucdiv_mat)[[1]], ncol = dim(nucdiv_mat)[[2]] - 1)

for (i in seq(1, dim(rpkm_mat)[[2]] - 1)) {
  
  rpkm_shift_mat[ , i] <- map2_dbl(rpkm_mat[ , i + 1], rpkm_mat[ , i], `/`)

}

colnames(rpkm_shift_mat) <- sapply(seq(1, 8), function(x) paste0("Ratio_", x))
df <- df %>%
  cbind(as.data.frame(rpkm_shift_mat))
