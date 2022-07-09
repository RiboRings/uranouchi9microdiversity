plot_distribution <- function(cur_stat) {
  
  ggplot(df, aes(x = eval(parse(text = cur_stat)))) +
    geom_histogram(bins = 40) +
    labs(x = cur_stat,
         y = "Count") +
    theme_bw()
  
}

snv_dist <- function(arg1, arg2) {
  
  return(sum(abs(arg1 - arg2)) / length(arg1))
  
}

freq_mapper <- function(cur_pos) {
  
  freq_array <- cur_pos %>%
    select(contains(paste0("Freq", cur_pos$con_base))) %>%
    as.vector() %>%
    as.numeric()
  
  return(freq_array)
  
}

max_mapper <- function(cur_pos, pattern) {
  
  max_val <- cur_pos %>%
    select(contains(paste(pattern, cur_pos$AbundMaxIdx, sep = "_"))) %>%
    as.vector() %>%
    as.numeric()
  
  return(max_val)
  
}

sample_loader <- function(sample_name, method = read_csv) {
  
  sample_df <- method(sample_name) %>%
    mutate(sample = gsub(".csv", "", gsub("data/", "", sample_name))) %>%
    left_join(mapping) %>%
    filter(genome == "metabat2bin_835.fna") %>%
    
  return(sample_df)
  
}
