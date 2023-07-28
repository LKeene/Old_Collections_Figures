library(tidyverse)
library(ggtree)
library(treeio)

# Practice to figure things out!
practice <- read.beast("BEAST_trees/galbut_RNA1_practice.mcc.tre")

df <- as.data.frame(tree[["data"]][["label"]]) 

nodes <- ggtree(practice, mrsd="2022-01-01", color="grey30", size=0.5) +
  geom_text(aes(label=node), hjust=-.3)
nodes

tree <- ggtree(practice, mrsd="2022-01-01", color="grey30", size=0.5) 

tree2 <- tree +
  theme_tree2() +
  geom_tippoint(aes(shape = Location), size = 3) +
  scale_shape_manual(values = c(0, 8, 16, 4, 17)) +
  geom_highlight(node = 126, fill = "steelblue", alpha = 0.5) +
  geom_highlight(node = 69, fill = "seagreen4", alpha = 0.5) +
  geom_highlight(node = 126, fill = "steelblue", alpha = 0.5) +
  geom_highlight(node = 68, fill = "steelblue", alpha = 0.5) +
  geom_range("height_0.95_HPD", color="gray50", alpha=0.2, size=1)
tree2

rna1_metadata <- read_csv("metadata/galbut_RNA1_metadata.csv")

rna1_metadata <- as.data.frame(rna1_metadata) 

tree3 <- tree2 %<+% rna1_metadata +
  geom_tiplab(aes(label = date, color = factor(date)), show.legend = F, 
              hjust = -0.5, parse = T) +
  labs(x = "Year")
tree3

ggsave("plots/galbut_practice_tree.pdf", width = 15, height = 10)

# Galbut RNA1 Tree
g1 <- read.beast("BEAST_trees/combined_galbut_RNA1.mcc.tre")

n_n <- ggtree(g1, mrsd="2022-01-01", color="grey30", size=0.5) +
  geom_text(aes(label=node), hjust=-.3) +
  geom_tippoint(aes(shape = Host))
n_n

#check to make sure labels are the same
labels <- as.data.frame(n_n[["data"]][["label"]]) 

g1_t1 <- ggtree(g1, mrsd="2022-01-01", color="grey30", size=0.5) 

g1_t2 <- g1_t1 +
  theme_tree2() +
  geom_tippoint(aes(shape = Location), size = 3) +
  scale_shape_manual(values = c(0, 8, 16, 4, 17)) +
  geom_highlight(node = 126, fill = "steelblue", alpha = 0.5) +
  geom_highlight(node = 69, fill = "seagreen4", alpha = 0.5) +
  geom_highlight(node = 126, fill = "steelblue", alpha = 0.5) +
  geom_highlight(node = 68, fill = "steelblue", alpha = 0.5) +
  geom_range("height_0.95_HPD", color="gray50", alpha=0.2, size=1)
g1_t2

rna1_metadata <- read_csv("metadata/galbut_RNA1_metadata.csv")

rna1_metadata <- as.data.frame(rna1_metadata) 

g1_t3 <- g1_t2 %<+% rna1_metadata +
  geom_tiplab(aes(label = date, color = factor(date)), show.legend = F, 
              hjust = -0.5, parse = T) +
  labs(x = "Year")
g1_t3

ggsave("plots/galbut_RNA1_tree.pdf", width = 15, height = 10)











