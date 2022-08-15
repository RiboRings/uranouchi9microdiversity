gen_list <- data.frame(genome = paste0(filtered_df$genome, ".fna"), AbundMaxIdx = filtered_df$AbundMaxIdx)
mapping <- read_tsv("data/mapping.stb", col_names = c("scaffold", "genome")) %>%
  right_join(gen_list)
gene_info_files <- list.files(path = "data", pattern = "gene_info*")

gene_info_list <- lapply(paste("data", gene_info_files, sep = "/"),
                         sample_loader,
                         method = read_tsv,
                         gen_list = gen_list$genome)

gene_info_df <- gene_info_list %>%
  merge_all()

gene_info_df <- gene_info_df %>%
  group_by(genome, gene) %>%
  summarise(pNpS = max(pNpS_variants, na.rm = TRUE)) %>%
  arrange(desc(pNpS))

man_annot <- read_csv2("data/manual_annot.txt")
sig_annot <- read_tsv("data/sig_func_annot.tsv", col_names = c("significance", "gene", "KO", "threshold", "score", "Evalue", "product"))
func_cats <- read_tsv("../KEGG_pathway_ko_uniq.txt")

func_annot <- man_annot %>%
  full_join(sig_annot) %>%
  left_join(func_cats, by = c("KO" = "ko"))

big_df <- gene_info_df %>%
  left_join(func_annot)

big_df$level2_pathway_name[big_df$KO %in% c("PF13568", "PF18998", "PF03257")] <- "Membrane transport"
big_df <- big_df %>%
  filter(!is.na(level2_pathway_name) & !is.na(pNpS) & level2_pathway_name != "Poorly characterized")
big_df$pNpS[big_df$pNpS == -Inf] <- 0

for (cur_genome in gen_list) {
  
  info_df <- filtered_df %>%
    filter(genome == gsub(".fna", "", cur_genome)) %>%
    transmute(Tax = paste(Domain, Phylum, Class, Order, Family, Genus, Species, sep = ";"),
              Contamination,
              Completeness)
  
  p1 <- ggplot(big_df, aes(x = pNpS, y = reorder(level2_pathway_name, pNpS, median, na.rm = TRUE))) +
    geom_boxplot() +
    scale_x_continuous(limits = c(0, 3),
                       breaks = seq(0, 3, by = 0.5)) +
    labs(y = "Functional Category") +
    facet_grid(~ genome) +
    theme_bw() +
    theme(axis.title.x = element_blank())
  
  top_genes_df <- func_annot %>%
    filter(!is.na(pNpS),
           product != "uncharacterized protein")
  
  top_genes_df <- top_genes_df[c(1:30, seq(nrow(top_genes_df) - 30, nrow(top_genes_df))), ] %>%
    mutate(product = ifelse(nchar(product) > 100, KO, product))

  p2 <- ggplot(top_genes_df, aes(x = pNpS,
                                 y = reorder(product, pNpS, median, na.rm = TRUE),
                                 fill = site)) +
    geom_col(aes(pNpS)) +
    labs(y = "Gene",
         fill = "Product Type") +
    scale_x_continuous(limits = c(0, 3),
                       breaks = seq(0, 3, by = 0.5)) +
    theme_classic()
  
  p <- p1 / p2 +
    plot_layout(guides = "collect") +
    plot_annotation(title = "Functional Annotation",
                    subtitle = paste0("MAG: ", info_df$Tax,
                                      " Completeness: ", info_df$Completeness,
                                      " Contamination: ", info_df$Contamination))
  
  ggsave(paste0("func_annot", gsub(".fna", "", gsub("metabat2bin_", "", cur_genome)), ".pdf"),
         device = "pdf",
         plot = p,
         width = 20,
         height = 10,
         path = "results")
  
}

func_mat <- big_df %>%
  group_by(genome, level2_pathway_name) %>%
  summarise(MeanpNpS = mean(pNpS)) %>%
  pivot_wider(names_from = level2_pathway_name,
              values_from = MeanpNpS,
              values_fill = 0) %>%
  column_to_rownames(var = "genome") %>%
  as.matrix()

func_mat[func_mat > 3] <- NA

my_tax <- filtered_df$Tax[filtered_df$genome %in% gsub(".fna", "", rownames(func_mat))]
my_tax_order <- order(as.numeric(gsub(".*\\;", "", gsub("o__", "", my_tax))))
my_mat_order <- order(as.numeric(gsub("metabat2bin_", "", gsub(".fna", "", rownames(func_mat)))))

my_tax <- my_tax[my_tax_order]
func_mat <- func_mat[my_mat_order, ]

rownames(func_mat) <- my_tax

pdf("results/func_hm.pdf",
    width = 10,
    height = 20)
draw(Heatmap(func_mat,
        name = "Mean pNpS",
        heatmap_legend_param = list(at = seq(min(func_mat, na.rm = TRUE), max(func_mat, na.rm = TRUE), by = 0.5)),
        row_title = "MAG",
        column_title = "Functional Group",
        row_names_gp = gpar(fontsize = 4),
        column_names_gp = gpar(fontsize = 8),
        column_names_rot = 45))
dev.off()
