plot_distribution <- function(cur_stat) {
  
  ggplot(df, aes(x = eval(parse(text = cur_stat)))) +
    geom_histogram(bins = 40) +
    labs(x = cur_stat,
         y = "Count") +
    theme_bw()
  
}

snv_dist <- function(arg1, arg2) {
  
  return(round(sum(abs(arg1 - arg2)), digits = 1))
  
}
