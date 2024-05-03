#!/usr/bin/env Rscript

# This code block sets up input arguments to either come from the command line
# (if running from the pipeline, so not interactively) or to use expected default values
# (if running in interactive mode, for example in RStudio for troubleshooting
# or development).
#
if (!interactive()) {
  # if running from Rscript
  args = commandArgs(trailingOnly=TRUE)
  tsv_input            = args[1]
  metadata_file        = args[2]
  refseq_metadata_file = args[3]
  R_lib_dir            = args[4]
  R_script_dir         = args[5]
  output_dir           = "./"
} else {
  # if running via RStudio
  tsv_input             = "../results/collected_strand_bias.tsv"
  metadata_file         = "../results/collected_metadata.csv" 
  refseq_metadata_file  = "../../metadata/virus_refseq_metadata.csv"
  R_lib_dir             = "../lib/R/"
  R_script_dir          = "../../scripts/"
  output_dir            = "../results/process/"
}

# these libraries are part of the tidyverse, so will be availabile in the
# tidyverse singularity image we are using (or analogous conda env)
library(tidyverse)

# these libraries are not part of the standard tidyverse, so may have to load it 
# from a specified path
# either from pipeline's R lib dir or from R environment
if (R_lib_dir != "NA") {
  library(rstatix, lib.loc=R_lib_dir)
  library(ggpubr, lib.loc=R_lib_dir)
  
} else {
  library(rstatix)
  library(ggpubr)
}
# read in color info
source(paste0(R_script_dir, "/plot_colors.R"))
fancy_color_scale <- c(old_collection_color, experimental_dried_color, fresh_frozen_color)

# read in the file with the mismatch data for all datasets
datasets <- read.delim(tsv_input, header=F, sep="\t")
colnames(datasets) <- c("sample_id", "sam_file", "refseq", "fwd_count", "rev_count", "fraction_positive_strand")

# datasets %>% group_by(sample_id) %>% summarize()

# only keep datasets with enough mapping reads
datasets <- datasets %>% mutate(total_count = fwd_count + rev_count)
datasets <- datasets %>% filter(total_count > 10)

# read in metadata
metadata <- read.delim(metadata_file, sep=",", header=T)

# virus reference sequence metadata 
refseq_metadata <- read.delim(refseq_metadata_file, sep=",", header=T)

# join in metadata
df <- left_join(datasets, metadata, by="sample_id")
df <- left_join(df, refseq_metadata, by="refseq")
# filter(refseq_metadata, row_number() == 31)
# filter(df, row_number() == 5)


# remove seqs for which virus not defined
missing_virus <- df %>% filter(is.na(virus)) %>% group_by(refseq) %>% summarize()
if (nrow(missing_virus) > 0) {
  message("WARNING: missing virus reference sequence metadata for the following datasets:")
  print (missing_virus, n=50)
}

# relevel categories
df$sample_type <- fct_relevel(df$sample_type, "Old_Collection", "Experimental_dried", "Fresh_frozen")

# rename categories

# keep non-pretty names
df$sample_type_og <- df$sample_type

df$sample_type <- 
  recode(df$sample_type, 
         Old_Collection     = "Old\ncollections", 
         Experimental_dried = "Experimental\ndried",
         Fresh_frozen       = "Fresh\nfrozen")


# averages by group
group_averages <- df %>% 
  group_by(sample_type_og, virus) %>% 
  summarize(mean_pos_strand = mean(fraction_positive_strand),
            median_pos_strand = median(fraction_positive_strand),
            stdev_pos_strand = sd(fraction_positive_strand),
            .groups = "drop")

# create text for paper 

# what fraction of RNA is from +strand: galbut virus
galbut_fresh <- group_averages %>% 
  filter(virus == "Galbut_virus" & sample_type_og == "Fresh_frozen") %>% 
  pull(median_pos_strand)

galbut_dried <- group_averages %>% 
  filter(virus == "Galbut_virus" & sample_type_og == "Experimental_dried") %>% 
  pull(median_pos_strand)

galbut_oc <- group_averages %>% 
  filter(virus == "Galbut_virus" & sample_type_og == "Old_Collection") %>% 
  pull(median_pos_strand)


printf <- function(...) invisible(print(sprintf(...)))

printf("a median of %0.1f%% of galbut virus reads derived from +strand RNA", galbut_fresh * 100)
printf("In experimental dried samples (≤72 weeks old) an average of %0.1f%% of reads were from +strand RNA, and this value fell to %0.1f%% in museum specimens.", galbut_dried * 100, galbut_oc * 100)

# make a combined label virus name and genome type
df <- df %>% mutate(virus_name_type = paste0(virus, "\n(", virus_type, ")"))
df$virus_name_type <- str_replace(df$virus_name_type, "plus_", "+")
df$virus_name_type <- str_replace(df$virus_name_type, "minus_", "-")

# do stats
# only keep groups with > 3 observations
df_stats <- df %>% group_by(virus, sample_type) %>% mutate(n=n()) %>% filter(n>3)
# only keep viruses with >= 2 groups
viruses_with_enough <- df_stats %>% group_by(virus, sample_type) %>% 
  summarize() %>% group_by(virus) %>% 
  summarize (n=n(), .groups="drop") %>% filter(n>1) %>% pull(virus)
df_stats_enough <- df_stats %>% filter(virus %in% viruses_with_enough)
df_enough <- df %>% filter(virus %in% viruses_with_enough)

# is data normally distributed?  Do Shapiro test
df_shapiro_enough <- df_enough %>% group_by(sample_type, virus) %>% shapiro_test(fraction_positive_strand)
filter(df_shapiro_enough, p < 0.05)

# not normally distributed, so should do non-parametric test

# df_t_test <- df_stats_enough %>% group_by(virus) %>% t_test(fraction_positive_strand ~ sample_type)
# df_t_test_plot <- df_t_test %>% add_xy_position(x = "sample_type")

# df_t_test_type <- df_stats_enough %>% group_by(virus_name_type) %>% t_test(fraction_positive_strand ~ sample_type)
# df_t_test_type_plot <- df_t_test_type %>% add_xy_position(x = "sample_type")

df_wilcoxon      <- df_stats_enough %>% 
  group_by(virus) %>% 
  wilcox_test(fraction_positive_strand ~ sample_type)

df_wilcoxon_plot <- df_wilcoxon %>% add_xy_position(x = "sample_type")

df_wilcoxon_type      <- df_stats_enough %>% 
  group_by(virus_name_type) %>% 
  wilcox_test(fraction_positive_strand ~ sample_type)
df_wilcoxon_type_plot <- df_wilcoxon_type %>% add_xy_position(x = "sample_type")


all_p <- ggplot(df) +
  geom_jitter(aes(x=sample_type, y=fraction_positive_strand, fill=sample_type), 
              shape=21, color="black", stroke=0.2, size=1.5, height=0, width=0.25, alpha=0.5) +
  geom_boxplot(aes(x=sample_type, y=fraction_positive_strand, color=sample_type), 
               size=0.25, outlier.shape = NA, fill=NA) +
  theme_bw(base_size=12) +
  scale_y_continuous(limits=c(0,1.1), breaks=seq(0, 1, 0.25)) +
  xlab("") +
  ylab("Fraction of reads from +sense RNA") +
  scale_color_manual(values = fancy_color_scale) +
  scale_fill_manual (values = fancy_color_scale) +
  facet_wrap(~virus_name_type) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme(legend.position = "none") 

all_p

all_p +
  stat_pvalue_manual(df_wilcoxon_type, y.position = c(1.05, 1.09, 1.05),
                     tip.length = 0.01, bracket.shorten=0.1, size=3)
ggsave(paste0(output_dir, "/Fig_X_fraction_plus_strand_all.pdf"), width=10, height=7, units="in")


# plot a plot with viruses that have observations in all 3 categories
enough_p <- ggplot(df_enough) +
  geom_jitter(aes(x=sample_type, y=fraction_positive_strand, fill=sample_type), 
              shape=21, color="black", stroke=0.2, size=2.5, height=0, width=0.25, alpha=0.5) +
  geom_boxplot(aes(x=sample_type, y=fraction_positive_strand, color=sample_type), 
               size=0.25, outlier.shape = NA, fill=NA) +
  theme_bw(base_size=14) +
  scale_y_continuous(limits=c(0,1.1), breaks=seq(0, 1, 0.25)) +
  xlab("") +
  ylab("Fraction of reads from +sense RNA") +
  scale_color_manual(values = fancy_color_scale) +
  scale_fill_manual (values = fancy_color_scale) +
  facet_wrap(~virus) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme(legend.position = "none")  

enough_p
  
enough_p + 
  stat_pvalue_manual(df_wilcoxon, y.position = c(1.05, 1.09, 1.05), tip.length = 0.01, bracket.shorten=0.1, size=3)

ggsave(paste0(output_dir, "/Fig_X_fraction_plus_strand_galbut.pdf"), width=7, height=4, units="in")


