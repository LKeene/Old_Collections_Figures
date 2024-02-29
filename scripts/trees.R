library(tidyverse)
library(ggtree)
library(treeio)
library(gridExtra)

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
  scale_color_manual(labels = c(substitute(paste(italic("D. melanogaster"))), 
                                substitute(paste(italic("D. simulans"))),
                                substitute(paste(italic("Diptera")))),
                     values = c("seagreen", "darkblue", "orange2")) +
  labs(color = "Species")
g_r1_v2

#####
# RNA1
galbut_RNA1 <- read.iqtree("IQ-Tree/updated_galbut_RdRp_NoFoCo_alignment.phy.treefile")

rna1_metadata <- read_csv("metadata/galbut_rna1_ML_metadata.csv")
rna1_metadata <- as.data.frame(rna1_metadata)

labels <- as.data.frame(n_n[["data"]][["label"]]) 

n_n <- ggtree(galbut_RNA1, color="grey30", size=0.5) +
  geom_text(aes(label=node), hjust=-.3)
n_n

g_r1 <- ggtree(galbut_RNA1)
g_r1

g_r1_v2 <- g_r1 %<+% rna1_metadata +
  geom_tiplab(aes(label = id, color = species), 
              hjust = -0.2, parse = T, size = 2) +
  scale_color_manual(labels = c(substitute(paste(italic("D. melanogaster"))),
                                substitute(paste(italic("D. simulans"))),
                                substitute(paste(italic("Diptera")))),
                     values = c("seagreen", "darkblue", "orange2")) +
  geom_cladelab(node = 57, label = "Clade A", offset = 0.02) +
  geom_cladelab(node = 95, label = "D. simulans clade", offset = 0.02) +
  geom_cladelab(node = 93, label = "Clade B", offset = 0.01) +
  labs(color = "Species")

g_r1_v2
ggsave("plots/galbut_RdRp_tree.pdf", width = 15, height = 10)
ggsave("plots/galbut_RdRp_tree.jpg", width = 15, height = 10)

# RNA2
galbut_RNA2 <- read.iqtree("IQ-Tree/updated_galbut_capsid_NoFoCo_alignment.phy.treefile")

rna2_metadata <- read_csv("metadata/galbut_rna2_ML_metadata.csv")
rna2_metadata <- as.data.frame(rna2_metadata)


n_n2 <- ggtree(galbut_RNA2, color="grey30", size=0.5) +
  geom_text(aes(label=node), hjust=-.3)
n_n2

labels <- as.data.frame(n_n2[["data"]][["label"]]) 

g_r2 <- ggtree(galbut_RNA2)
g_r2

g_r2_v2 <- g_r2 %<+% rna2_metadata +
  geom_tiplab(aes(label = id, color = species), 
              hjust = -0.2, parse = T, size = 2) +
  scale_color_manual(labels = c(substitute(paste(italic("D. melanogaster"))),
                                substitute(paste(italic("D. simulans")))),
                     values = c("seagreen", "darkblue")) +
  geom_cladelab(node = 46, label = "Clade A", offset = 0.01) +
  geom_cladelab(node = 45, label = "D. simulans clade", offset = 0.01) +
  labs(color = "Species")

g_r2_v2
ggsave("plots/galbut_capsid_tree.pdf", width = 15, height = 10)
ggsave("plots/galbut_capsid_tree.jpg", width = 15, height = 10)

# RNA3
galbut_RNA3 <- read.iqtree("IQ-Tree/galbut_RNA3_NoFoCo_alignment.phy.treefile")

rna3_metadata <- read_csv("metadata/galbut_rna3_ML_metadata.csv")
rna3_metadata <- as.data.frame(rna3_metadata)

n_n3 <- ggtree(galbut_RNA3, color="grey30", size=0.5) +
  geom_text(aes(label=node), hjust=-.3)
n_n3

labels <- as.data.frame(n_n3[["data"]][["label"]])

g_r3 <- ggtree(galbut_RNA3)
g_r3

g_r3_v2 <- g_r3 %<+% rna3_metadata +
  geom_tiplab(aes(label = id, color = species), 
              hjust = -0.2, parse = T, size = 2) +
  scale_color_manual(labels = c(substitute(paste(italic("D. melanogaster"))),
                                substitute(paste(italic("D. simulans")))),
                     values = c("seagreen", "darkblue")) +
  geom_cladelab(node = 45, label = "Clade A", offset = 0.02) +
  geom_cladelab(node = 82, label = "D. simulans clade", offset = 0.02) +
  labs(color = "Species")

g_r3_v2
ggsave("plots/galbut_RNA3_tree.pdf", width = 15, height = 10)
ggsave("plots/galbut_RNA3_tree.jpg", width = 15, height = 10)

## Sigmavirus

sigma <- read.iqtree("IQ-Tree/sigmavirus_outgroup_alignment.phy.treefile")

tree_sig <- ggtree(sigma, color = "grey30", size = 0.5) +
  geom_text(aes(label = node), hjust = -0.3)
tree_sig  

labels <- as.data.frame(tree_sig[["data"]][["label"]])  
labels 

sigma_NO <- read.iqtree("IQ-Tree/sigmavirus_NoOutgroup_alignment.phy.treefile")

sigma_metadata <- read_csv("metadata/sigma_ML_metadata.csv")
sigma_metadata <- as.data.frame(sigma_metadata)  

tree_sig2 <- ggtree(sigma_NO)  
tree_sig2  

tree_sig3 <- tree_sig2 %<+% sigma_metadata +
  geom_tiplab(aes(label = id, color = year), 
              hjust = -0.01, size = 2.5) +
  geom_tippoint(aes(shape = location), size = 2) +
  labs(color = "Year")


tree_sig3
ggsave("plots/sigma_tree.pdf", width = 15, height = 10)
ggsave("plots/sigma_tree.jpg", width = 15, height = 10)



