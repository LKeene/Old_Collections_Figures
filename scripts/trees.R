library(tidyverse)
library(ggtree)
library(treeio)
library(gridExtra)

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
  geom_highlight(node = 51, fill = "deeppink", alpha = 0.5) +
  geom_highlight(node = 53, fill = "deeppink", alpha = 0.5) +
  geom_highlight(node = 52, fill = "deeppink", alpha = 0.5) +
  geom_highlight(node = 68, fill = "deeppink", alpha = 0.5) +
  geom_highlight(node = 54, fill = "deeppink", alpha = 0.5) +
  geom_highlight(node = 13, fill = "deeppink", alpha = 0.5) +
  geom_range("height_0.95_HPD", color="gray50", alpha=0.2, size=1)
g1_t2

rna1_metadata <- read_csv("metadata/galbut_RNA1_metadata.csv")

rna1_metadata <- as.data.frame(rna1_metadata) 

g1_t3 <- g1_t2 %<+% rna1_metadata +
  geom_tiplab(aes(label = date, color = species), 
              hjust = -0.4, parse = T, size = 3.5) +
  scale_color_manual(labels = c(substitute(paste(italic("Diptera"))), 
                                substitute(paste(italic("D. melanogaster"))), 
                                substitute(paste(italic("D. simulans")))),
                     values = c("seagreen", "darkblue", "orange2")) +
  labs(x = "Year", color = "Species")
g1_t3

ggsave("plots/galbut_RNA1_tree.pdf", width = 15, height = 10)


# Galbut RNA2 Tree
g2 <- read.beast("BEAST_trees/Galbut_RNA2_2.mcc.tre")

n_n <- ggtree(g2, mrsd="2022-01-01", color="grey30", size=0.5) +
  geom_text(aes(label=node), hjust=-.3) +
  geom_tippoint(aes(shape = Host))
n_n

#check to make sure labels are the same
labels <- as.data.frame(n_n[["data"]][["label"]]) 

g2_t1 <- ggtree(g2, mrsd="2022-01-01", color="grey30", size=0.5) 

g2_t2 <- g2_t1 +
  theme_tree2() +
  geom_tippoint(aes(shape = Location), size = 3) +
  scale_shape_manual(values = c(0, 8, 18, 4, 17)) +
  geom_highlight(node = 60, fill = "deeppink", alpha = 0.5) +
  geom_highlight(node = 57, fill = "deeppink", alpha = 0.5) +
  geom_highlight(node = 58, fill = "deeppink", alpha = 0.5) +
  geom_highlight(node = 56, fill = "deeppink", alpha = 0.5) +
  geom_highlight(node = 59, fill = "deeppink", alpha = 0.5) +
  geom_highlight(node = 12, fill = "deeppink", alpha = 0.5) +
  geom_highlight(node = 13, fill = "deeppink", alpha = 0.5) +
  geom_range("height_0.95_HPD", color="gray50", alpha=0.2, size=1)
g2_t2

rna2_metadata <- read_csv("metadata/galbut_RNA2_metadata.csv")

rna2_metadata <- as.data.frame(rna2_metadata) 

g2_t3 <- g2_t2 %<+% rna2_metadata +
  geom_tiplab(aes(label = date, color = species), 
              hjust = -0.4, parse = T, size = 3.5) +
  scale_color_manual(labels = c(substitute(paste(italic("D. melanogaster"))), 
                                substitute(paste(italic("D. simulans")))),
                     values = c("darkblue", "orange2")) +
  labs(x = "Year", color = "Species")
g2_t3

ggsave("plots/galbut_RNA2_tree.pdf", width = 15, height = 10)


# Galbut RNA3 Tree
g3 <- read.beast("BEAST_trees/Combined_galbut_RNA3.mcc.tre")

n_n <- ggtree(g3, mrsd="2022-01-01", color="grey30", size=0.5) +
  geom_text(aes(label=node), hjust=-.3) +
  geom_tippoint(aes(shape = Host))
n_n

#check to make sure labels are the same
labels <- as.data.frame(n_n[["data"]][["label"]]) 

g3_t1 <- ggtree(g3, mrsd="2022-01-01", color="grey30", size=0.5) 

g3_t2 <- g3_t1 +
  theme_tree2() +
  geom_tippoint(aes(shape = Location), size = 3) +
  scale_shape_manual(values = c(0, 8, 18, 4, 17)) +
  geom_highlight(node = 55, fill = "deeppink", alpha = 0.5) +
  geom_highlight(node = 12, fill = "deeppink", alpha = 0.5) +
  geom_highlight(node = 53, fill = "deeppink", alpha = 0.5) +
  geom_highlight(node = 54, fill = "deeppink", alpha = 0.5) +
  geom_highlight(node = 52, fill = "deeppink", alpha = 0.5) +
  geom_range("height_0.95_HPD", color="gray50", alpha=0.2, size=1)
g3_t2

rna3_metadata <- read_csv("metadata/galbut_RNA3_metadata.csv")

rna3_metadata <- as.data.frame(rna3_metadata) 

g3_t3 <- g3_t2 %<+% rna3_metadata +
  geom_tiplab(aes(label = date, color = species), 
              hjust = -0.4, parse = T, size = 3.5) +
  scale_color_manual(labels = c(substitute(paste(italic("D. melanogaster"))), 
                                substitute(paste(italic("D. simulans")))),
                     values = c("darkblue", "orange2")) +
  labs(x = "Year", color = "Species")
g3_t3

ggsave("plots/galbut_RNA3_tree.pdf", width = 15, height = 10)


# Galbut RNA1 Tree- NO Old Collections Samples!
g1_n <- read.beast("BEAST_trees/combined_galbut_RNA1_NO_OC.mcc.tre")

n_n <- ggtree(g1_n, mrsd="2022-01-01", color="grey30", size=0.5) +
  geom_text(aes(label=node), hjust=-.3) +
  geom_tippoint(aes(shape = Host))
n_n

#check to make sure labels are the same
labels <- as.data.frame(n_n[["data"]][["label"]]) 

g1_n_t1 <- ggtree(g1_n, mrsd="2022-01-01", color="grey30", size=0.5) 

g1_n_t2 <- g1_n_t1 +
  theme_tree2() +
  geom_tippoint(aes(shape = Location), size = 3) +
  scale_shape_manual(values = c(0, 8, 16, 4, 17)) +
  geom_range("height_0.95_HPD", color="gray50", alpha=0.2, size=1)
g1_n_t2

rna1_metadata <- read_csv("metadata/galbut_RNA1__NO_OC_metadata.csv")

rna1_metadata <- as.data.frame(rna1_metadata) 

g1_n_t3 <- g1_n_t2 %<+% rna1_metadata +
  geom_tiplab(aes(label = date, color = species), 
              hjust = -0.4, parse = T, size = 3.5) +
  scale_color_manual(labels = c(substitute(paste(italic("Diptera"))), 
                                substitute(paste(italic("D. melanogaster"))), 
                                substitute(paste(italic("D. simulans")))),
                     values = c("seagreen", "darkblue", "orange2")) +
  labs(x = "Year", color = "Species")
g1_n_t3

ggsave("plots/galbut_RNA1_NO_OC_tree.pdf", width = 15, height = 10)

# Galbut RNA2 Tree- NO Old Collections Samples
g2_n <- read.beast("BEAST_trees/galbut_RNA2_NO_OC.mcc.tre")

n_n <- ggtree(g2_n, mrsd="2022-01-01", color="grey30", size=0.5) +
  geom_text(aes(label=node), hjust=-.3)
n_n

#check to make sure labels are the same
labels <- as.data.frame(n_n[["data"]][["label"]]) 

g2_n_t1 <- ggtree(g2_n, mrsd="2022-01-01", color="grey30", size=0.5) 

g2_n_t2 <- g2_n_t1 +
  theme_tree2() +
  geom_tippoint(aes(shape = Location), size = 3) +
  scale_shape_manual(values = c(0, 8, 4, 17)) +
  geom_range("height_0.95_HPD", color="gray50", alpha=0.2, size=1)
g2_n_t2

rna2_metadata <- read_csv("metadata/galbut_RNA2_NO_OC_metadata.csv")

rna2_metadata <- as.data.frame(rna2_metadata) 

g2_n_t3 <- g2_n_t2 %<+% rna2_metadata +
  geom_tiplab(aes(label = date, color = species), 
              hjust = -0.4, parse = T, size = 3.5) +
  scale_color_manual(labels = c(substitute(paste(italic("D. melanogaster")))),
                     values = c("darkblue")) +
  labs(x = "Year", color = "Species")
g2_n_t3

ggsave("plots/galbut_RNA2_tree_NO_OC.pdf", width = 15, height = 10)


# Galbut RNA3 Tree- NO Old Collections Samples!
g3_n <- read.beast("BEAST_trees/galbut_RNA3_NO_OC.mcc.tre")

n_n <- ggtree(g3_n, mrsd="2022-01-01", color="grey30", size=0.5) +
  geom_text(aes(label=node), hjust=-.3)
n_n

#check to make sure labels are the same
labels <- as.data.frame(n_n[["data"]][["label"]]) 

g3_n_t1 <- ggtree(g3_n, mrsd="2022-01-01", color="grey30", size=0.5) 

g3_n_t2 <- g3_n_t1 +
  theme_tree2() +
  geom_tippoint(aes(shape = Location), size = 3) +
  scale_shape_manual(values = c(0, 8, 4, 17)) +
  geom_range("height_0.95_HPD", color="gray50", alpha=0.2, size=1)
g3_n_t2

rna3_metadata <- read_csv("metadata/galbut_RNA3_NO_oC_metadata.csv")

rna3_metadata <- as.data.frame(rna3_metadata) 

g3_n_t3 <- g3_n_t2 %<+% rna3_metadata +
  geom_tiplab(aes(label = date, color = species), 
              hjust = -0.4, parse = T, size = 3.5) +
  scale_color_manual(labels = c(substitute(paste(italic("D. melanogaster")))),
                     values = c("darkblue")) +
  labs(x = "Year", color = "Species")
g3_n_t3

ggsave("plots/galbut_RNA3_tree_NO_OC.pdf", width = 15, height = 10)

############

# Trees made with IQ-Tree
galbut_RNA1 <- read.iqtree("IQ-Tree/Galbut_RdRP_alignment_prelim.nex.treefile")

rna1_metadata <- read_csv("metadata/galbut_rna1_ML_metadata.csv")
rna1_metadata <- as.data.frame(rna1_metadata)

n_n <- ggtree(galbut_RNA1, color="grey30", size=0.5) +
  geom_text(aes(label=node), hjust=-.3)
n_n

g_r1 <- ggtree(galbut_RNA1) +
  theme_tree2()
g_r1

g_r1_v2 <- g_r1 %<+% rna1_metadata +
  geom_tiplab(aes(label = name, color = species), 
              hjust = -0.2, parse = T, size = 2) +
  scale_color_manual(labels = c(substitute(paste(italic("Diptera"))), 
                                substitute(paste(italic("D. melanogaster"))), 
                                substitute(paste(italic("D. simulans")))),
                     values = c("seagreen", "darkblue", "orange2")) +
  labs(color = "Species")
g_r1_v2
