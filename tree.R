library(ape)
library(ggtree)

ar_tree <- read.tree("data/gtdbtk.ar122.classify.tree")

ar_tree$tip.label[!str_detect(ar_tree$tip.label, "metabat2bin")] <- ""

ar_tree <- drop.tip(ar_tree, "")
plot(ar_tree)

test <- data.frame(TipLabel = ar_tree$tip.label)
arto <- data.frame(genome = df$genome,
                   order = gsub("c__", "", df$Class),
                   taxonomy = paste(df$Family,
                                    df$Genus,
                                    df$Species,
                                    sep = ";"))

test <- test %>%
  left_join(arto, by = c("TipLabel" = "genome"))
test$taxonomy[is.na(test$taxonomy)] <- ""

ar_tree$tip.label <- test$taxonomy

groupInfo <- split(ar_tree$tip.label, test$order)
names(groupInfo) <- unique(sort(test$order))
ar_tree <- groupOTU(ar_tree, groupInfo)

tree <- ggtree(ar_tree, aes(colour = group), layout = "circular") +
  geom_tiplab(size = 1.9, angle = 0) +
  scale_color_manual(values = c(`Poseidoniia_A` = "blue",
                                Huberarchaeia = "red")) +
  labs(colour = "Class")

ggsave("ar_tree.pdf",
       plot = tree,
       width = 10,
       height = 10,
       device = "pdf",
       path = "results")

## bac

bac_tree <- read.tree("data/gtdbtk.bac120.classify.tree")

bac_tree$tip.label[!str_detect(bac_tree$tip.label, "metabat2bin")] <- ""

bac_tree <- drop.tip(bac_tree, "")
plot(bac_tree)

test <- data.frame(TipLabel = bac_tree$tip.label)
arto <- data.frame(genome = df$genome,
                   order = gsub("p__", "", df$Phylum),
                   taxonomy = paste(df$Family,
                                    df$Genus,
                                    df$Species,
                                    sep = ";"))

test <- test %>%
  left_join(arto, by = c("TipLabel" = "genome"))
test$order[is.na(test$order)] <- "Other"
test$taxonomy[is.na(test$taxonomy)] <- ""

bac_tree$tip.label <- test$taxonomy

groupInfo <- split(bac_tree$tip.label, test$order)
names(groupInfo) <- unique(sort(test$order))
bac_tree <- groupOTU(bac_tree, groupInfo)

tree <- ggtree(bac_tree, aes(colour = group), layout = "circular") +
  geom_tiplab(size = 0.8, aes(angle = angle)) +
  labs(colour = "Phylum")

ggsave("bac_tree.pdf",
       plot = tree,
       width = 7,
       height = 7,
       device = "pdf",
       path = "results")

