hline_max1 <- data.frame(group = c("Dry", "Frozen"), target = c("La Jolla virus", 
                                                                "Nora virus",
                                                                "Thika virus"), 
                         hline_1 = c(30.81469, 27.84886, 32.70908, 28.27032,
                                     27.49015, 30.87427))


hline_min1 <- data.frame(group = c("Dry", "Frozen"), target = c("La Jolla virus", 
                                                                "Nora virus", 
                                                                "Thika virus"),
                         hline_2 = c(16.96781, 17.54493, 26.47017, 22.34384,
                                     18.78465, 22.54356))


# Mean line bar of La Jolla, Nora & Thika targets
ggplot(filter(fly_data3, target %in% c("La Jolla virus", "Nora virus", 
                                       "Thika virus"), 
              group != "Fresh"), aes(x=week)) + 
  geom_hline(data = hline_max1, aes(yintercept = hline_1), alpha = 0.5) +
  geom_hline(data = hline_min1, aes(yintercept = hline_2), alpha = 0.5) +
  geom_point(aes(y=ct, fill=target), shape = 21, size = 1.5, stroke = 0.1, 
             color = "black", alpha = 0.75) + 
  scale_fill_manual(values = c("turquoise3", "navyblue", "blue")) + 
  geom_line(aes(y=mean_ct, group=sample_length)) + 
  geom_errorbar(aes(ymin = (mean_ct - sd_ct), ymax = (mean_ct + sd_ct)), 
                width = 0.2, color = "grey50", alpha = 0.25) +
  theme_few(base_size = 11) + 
  facet_grid(target ~group, scales = "free_y") +
  labs(x = "Time Since Sample Collection (weeks)", y = "Mean Ct", fill = "Target")


# remove # to save plot
#ggsave("plots/Mean_LaJollaNoraThika.pdf", units = "in", width = 10, height = 8)


hline_max2 <- data.frame(group = c("Dry", "Dry", "Dry", "Dry", "Frozen", "Frozen",
                                   "Frozen", "Frozen"), 
                         sample_length = c("Galbut_long", "Galbut_short",
                                           "RpL32_long", "RpL32_short"), 
                         hline_1 = c(27.11656, 21.71760, 29.948694, 25.24503, 
                                     25.20976, NA, 23.35923, NA))


hline_min2 <- data.frame(group = c("Dry", "Dry", "Dry", "Dry", "Frozen", "Frozen",
                                   "Frozen", "Frozen"), 
                         sample_length = c("Galbut_long", "Galbut_short",
                                           "RpL32_long", "RpL32_short"),
                         hline_2 = c(18.23504, 16.04849, 20.2549, 20.49891, 
                                     17.29496, NA, 17.70815, NA))

# Mean line bar of galbut and RpL32 short and long
ggplot(filter(fly_data3, sample_length %in% c("Galbut_long", "Galbut_short", 
                                              "RpL32_long", "RpL32_short"), 
              group != "Fresh"), aes(x=week)) + 
  geom_hline(data = hline_max2, aes(yintercept = hline_1), alpha = 0.5) +
  geom_hline(data = hline_min2, aes(yintercept = hline_2), alpha = 0.5) +
  geom_point(aes(y=ct, fill=sample_length), shape = 21, size = 1.5, stroke = 0.1, 
             color = "black", alpha = 0.75) + 
  scale_fill_manual(values = c("turquoise3", "navyblue", "blue", "purple")) + 
  geom_line(aes(y=mean_ct, group=sample_length)) + 
  geom_errorbar(aes(ymin = (mean_ct - sd_ct), ymax = (mean_ct + sd_ct)), 
                width = 0.2, color = "grey50", alpha = 0.25) +
  theme_few(base_size = 11) + 
  facet_grid(sample_length ~ group, scales = "free_y") +
  labs(x = "Time Since Sample Collection (weeks)", y = "Mean Ct" )


# remove # to save plot
#ggsave("plots/Mean_GalbutRpL32.pdf", units = "in", width = 10, height = 8)

# relative change for all fly targets
ggplot(fly_fc, aes(x = week)) +
  geom_point(aes(y = fold_change, fill = group), shape = 21, size = 1, 
             stroke = 0.1, color = "black", alpha = 0.5) +
  scale_fill_manual(values = c("turquoise3", "purple")) +
  geom_line(aes(y = mean_fc, group = group, linetype = group)) +
  scale_y_continuous(trans = "log2") +
  theme_few(base_size = 11) +
  facet_wrap(~factor(target, levels = c("Galbut virus", "La Jolla virus", 
                                        "Nora virus", "Thika virus", 
                                        "RpL32 mRNA")), ncol = 2) +
  labs(x = 'Weeks After Collection', y = "Relative change", fill = "Target", 
       linetype = "Group")

# remove # based on what transformation was used in scale_y_continuous
#ggsave("plots/Relative_fly.pdf", units = "in", width = 10, height = 8)
#ggsave("plots/Relative_fly_log10.pdf", units = "in", width = 10, height = 8)
#ggsave("plots/Relative_fly_log2.pdf", units = "in", width = 10, height = 8)
#ggsave("plots/Relative_fly_sqrt.pdf", units = "in", width = 10, height = 8)


hline_max_mos <- data.frame(group = c("Dry", "Dry", "Dry", "Frozen", "Frozen", 
                                      "Frozen"), 
                            target = c("Verdadero virus", "Rennavirus", 
                                       "Actin mRNA"), 
                            hline_1_mos = c(27.94278, 29.54662, 29.98109,
                                            27.91520, 24.06199, 27.85027))

hline_min_mos <- data.frame(group = c("Dry", "Dry", "Dry", "Frozen", "Frozen", 
                                      "Frozen"), 
                            target = c("Verdadero virus", "Guadeloupe mosquito virus", 
                                       "Actin mRNA"),
                            hline_2_mos = c(22.04536, 19.4869, 24.19849,
                                            22.04304, 14.35967, 23.48702))

# Mean line bar of all Ae. aegypti targets 
ggplot(mosquito_data3, aes(x = week)) + 
  geom_hline(data = hline_max_mos, aes(yintercept = hline_1_mos), alpha = 0.3) +
  geom_hline(data = hline_min_mos, aes(yintercept = hline_2_mos), alpha = 0.3) +
  geom_point(aes(y = ct, fill = target), shape=21, size=1, stroke=0.1, 
             color="black", alpha = 0.5) + 
  scale_fill_manual(values = c("turquoise3", "purple", "green")) + 
  geom_line(aes(y=mean_ct, group=sample_group)) + 
  geom_errorbar(aes(ymin = (mean_ct - sd_ct), ymax = (mean_ct + sd_ct)), 
                width = 0.2, color = "grey50", alpha = 0.25) +
  theme_few(base_size = 11) +
  facet_grid(factor(target, levels = c("Verdadero virus", "Guadeloupe mosquito virus", 
                                       "Actin mRNA")) ~ group, scales = "free_y") +
  labs(x = "Time Since Sample Collection (weeks)", y = "Ct value", 
       fill = "Target")

# remove # to save plot
#ggsave("plots/Mean_Mosquito.pdf", units = "in", width = 10, height = 8)

# Relative change for all Ae. aegypti targets
ggplot(mosquito_data3, aes(x = week)) +
  geom_point(aes(y = fold_change, fill = group), shape = 21, size = 1, 
             stroke = 0.1, color = "black", alpha = 0.5) +
  scale_fill_manual(values = c("turquoise3", "purple")) +
  geom_line(aes(y = mean_fc, group = group, linetype = group)) +
  scale_y_continuous(trans = "log2") +
  theme_few(base_size = 11) +
  facet_wrap(~factor(target, levels = c("Verdadero virus", "Guadeloupe mosquito virus", 
                                        "Actin mRNA")), scales = "free_y",
             ncol = 1) +
  labs(x = 'Weeks After Collection', y = "Relative change", fill = "Target")

# remove # based on what transformation was used in scale_y_continuous
#ggsave("plots/Relative_Mosquito_log10.pdf", units = "in", width = 10, height = 8)
#ggsave("plots/Relative_Mosquito_log2.pdf", units = "in", width = 10, height = 8)
#ggsave("plots/Relative_Mosquito_sqrt.pdf", units = "in", width = 10, height = 8)
#ggsave("plots/Relative_Mosquito.pdf", units = "in", width = 10, height = 8)

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