#!/usr/bin/env Rscript

# This code block sets up input arguments to either come from the command line
# (if running from the pipeline, so not interactively) or to use expected default values 
# (if running in interactive mode, for example in RStudio for troubleshooting
# or development).  
#
if (!interactive()) {
  # if running from Rscript
  args = commandArgs(trailingOnly=TRUE)
  # lib_dir=args[1]
  tsv_input       = args[1]
  metadata_input  = args[2]
  R_lib_dir       = args[3]
  R_script_dir    = args[4]
  output_dir      = "./"
} else {
  # if running via RStudio
  # r_lib_dir = "../lib/R/"
  tsv_input       = "../results/process/collected_alfa_counts.tsv"
  metadata_input  = "../results/collect/collected_metadata.csv"
  R_lib_dir       = "../lib/R/"
  R_script_dir    = "../../scripts/"
  output_dir      = "../results/process/"
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
  library(patchwork, lib.loc=R_lib_dir)
  
} else {
  # in this case assuming these will be installed
  library(rstatix)
  library(ggpubr)
  library(patchwork)
}


# read in common color definitions
source(paste0(R_script_dir, "/plot_colors.R"))
fancy_color_scale <- c(fresh_frozen_color, experimental_dried_color, old_collection_color)

# read in metadata
metadata <- read.delim(metadata_input, sep=",", 
                       header=T, stringsAsFactors = F)

# convert to character type in case an integer-like sample ID (e.g. 1004277)
metadata$sample_id <- as.character(metadata$sample_id)

# define subsets of fresh-frozen datasets
metadata <- metadata %>% mutate(ff_type = case_when(
                                str_detect(sample_id, "SRX")  ~ "sra_dataset",
                                str_detect(sample_id, "\\d+Wk_[FM]") ~ "crisitunity",
                                str_detect(sample_id, "Col30") ~ "colony30",
                                sample_type != "Fresh_frozen" ~ NA_character_,
                                .default = "pos_ctrl"))


datasets <- read.delim(file = tsv_input, sep="\t")
colnames(datasets) <- c("sample_id", "type", "counts", "size_in_genome")
# convert to character type in case an integer-like sample ID (e.g. 1004277)
datasets$sample_id <- as.character(datasets$sample_id)

# plot count totals per dataset
datasets_per_sample <- datasets %>% group_by(sample_id) %>% summarize(total_counts = sum(counts))
ggplot(datasets_per_sample) +
  geom_histogram(aes(x=total_counts), bins=60) + 
  scale_x_log10() +
  xlab("# of fly genome mapping reads") +
  ylab("# datasets") +
  theme_this_paper()

ggsave(paste0(output_dir, "fly-mapping_reads_per_dataset.pdf"), units="in", height=6.5, width=9)

# filter out datasets with fewer than 1e7 total counts
too_few_counts <- datasets_per_sample %>% filter(total_counts < 1e7) %>% pull(sample_id)

datasets <- datasets %>% filter(!sample_id %in% too_few_counts)

# confirm that metadata exists for all datasets
dataset_names <- datasets %>% 
  group_by(sample_id) %>% summarize() %>% pull(sample_id)
missing_metadata <- !(dataset_names %in% metadata$sample_id)

if (any(missing_metadata)) {
  message("ERROR: missing metadata for the following datasets:")
  cat(dataset_names[missing_metadata])
}

# merge metadata into counts table
datasets <- left_join(datasets, metadata, by="sample_id")

# filter out non-melanogaster (simulans, etc.) samples since we mapped to D. mel genome
datasets <- datasets %>% filter(str_detect(species, "melan"))

# parse out category,biotype separated by a comma in column #2
type_split <- str_match(datasets$type, "(.*),(.*)")
datasets$category <- type_split[,2]
datasets$biotype  <- type_split[,3]

datasets <- datasets %>% select(-type)

datasets %>% group_by(sample_type) %>% summarize(n=n())

# combine counts for biotypes with >1 category
# sum counts for all categories within each biotype
summed_datasets <- datasets %>% group_by(sample_id, biotype) %>% summarize(counts = sum(counts))

# use summed values for plots, etc
datasets <- left_join(summed_datasets, metadata, by="sample_id")

# calculate normalized counts 
datasets <- datasets %>% 
  group_by(sample_id) %>% 
  mutate(fractional_counts = counts / sum(counts)) 

# reorder sample factors for display
# display samples in order of collection date
datasets$sample_id <- fct_reorder(datasets$sample_id, datasets$date_collected, min)

datasets$sample_type <- fct_relevel(datasets$sample_type, "Fresh_frozen", "Experimental_dried", "Old_Collection")

# don't plot negative control dataset (very few reads)
datasets <- datasets %>% filter(sample_type != "negative")

# recode data values for nicer output
datasets$sample_type <- 
  recode(datasets$sample_type, 
         Old_Collection     = "Museum\nsamples", 
         Experimental_dried = "Experimental\ndried",
         Fresh_frozen       = "Fresh\nfrozen")

datasets$biotype <- 
  recode(datasets$biotype, 
         opposite_strand      = "opposite strand",
         protein_coding       = "protein coding",
         pre_miRNA            = "pre-miRNA",
         transposable_element = "transpos. element")


# -----
# Stats
# -----

# are data normally distributed?

# calculate log fractional counts
datasets <- datasets %>% mutate(log_fractional_counts = log10(fractional_counts))

# this Shapiro tests indicates fractional counts are not normally distributed
df_shapiro     <- datasets %>% group_by(biotype, sample_type) %>% shapiro_test(fractional_counts)
filter(df_shapiro, p<0.05)

# this Shapiro tests indicates fractional counts are not log-normally distributed
df_shapiro_log <- datasets %>% group_by(biotype, sample_type) %>% shapiro_test(log_fractional_counts)
filter(df_shapiro_log, p<0.05)

# Wilcoxon test, since non-normal dist
df_wilcox <- datasets %>% group_by(biotype) %>% wilcox_test(fractional_counts ~ sample_type)
df_wilcox_plot <- df_wilcox %>% add_xy_position(x = "sample_type")

# Wilcoxon test on log10 transformed data to match log10 plot
df_wilcox_log <- datasets %>% group_by(biotype) %>% wilcox_test(log_fractional_counts ~ sample_type)
# df_wilcox_log_plot <- df_wilcox_log %>% add_xy_position(x = "sample_type", scales="free_y")
df_wilcox_log_plot <- df_wilcox_log %>% add_xy_position(x = "sample_type")

dataset_medians <- datasets %>% group_by(sample_type, biotype) %>% summarize(median_fraction = median(fractional_counts),
                                                                             sd_fraction = sd(fractional_counts))

biotype_median <- datasets %>% group_by(biotype) %>% summarize(median_fraction = median(fractional_counts),
                                                                sd_fraction = sd(fractional_counts))

datasets <-  left_join(datasets, biotype_median)

# split biotypes up into those that are common (>5% of RNA) and those that are rare (<5%)
abundant_rna <- filter(datasets, median_fraction >= 0.05)
rarer_rna    <- filter(datasets, median_fraction <  0.05)

# a generic plot function 
abundance_plot <- function(df_to_plot) {
  p <- ggplot(df_to_plot) +
    geom_jitter (aes(x=sample_type, y=fractional_counts, fill=sample_type), shape=21, color="black", stroke=0.2, size=1.5, height=0, width=0.25, alpha=0.5) +
    geom_boxplot(aes(x=sample_type, y=fractional_counts, color=sample_type), size=0.25, outlier.shape = NA, fill=NA) +
    theme_this_paper(base_size = 12) + 
    xlab("") +
    ylab("Fraction of reads mapping\nto indicated feature type") +
    scale_color_manual(values = fancy_color_scale) +
    scale_fill_manual (values = fancy_color_scale) +
    scale_y_continuous(limits=c(0,1)) +
    facet_wrap(~biotype, ncol=4) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    theme(legend.position = "none") +
    theme(strip.background = element_rect(linewidth = 0.5))
  
  p
}


# most abundant RNA types
abundant_p <- abundance_plot(abundant_rna)
abundant_p

# rare types of RNA
rarer_p <- abundance_plot(rarer_rna)
rarer_log_p <- rarer_p + scale_y_log10(limits=c(1e-8, 1e1))
rarer_log_p

# -----------------------------------------
# versions of plots with annotated p values
# -----------------------------------------

# subset df_wilcox_plot based on rarer biotypes
abundant_biotypes <- abundant_rna %>% group_by(biotype) %>% summarize() %>% pull(biotype)
df_wilcox_plot_abundant <- filter(df_wilcox_plot, biotype %in% abundant_biotypes)
y_positions <- rep(c(0.85,0.95,0.85), 4)
abundant_p_stats <- abundant_p + 
  stat_pvalue_manual(df_wilcox_plot_abundant, y.position = y_positions, tip.length = 0.01, bracket.shorten=0.1, size=3, bracket.nudget.y=0.1)
abundant_p_stats

# rarer biotypes
rarer_biotypes <- rarer_rna %>% group_by(biotype) %>% summarize() %>% pull(biotype)
df_wilcox_log_plot_rarer <- filter(df_wilcox_log_plot, biotype %in% rarer_biotypes)
y_positions <- rep(c(6e-2,7e-1,6e-2), 8)

rarer_log_p_stats <- rarer_log_p  +
  stat_pvalue_manual(df_wilcox_log_plot_rarer, y.position = y_positions, tip.length = 0.01, bracket.shorten=0.1, size=3)
rarer_log_p_stats

# plot combined abundant vs rarer biotypes
split_p_stats <- abundant_p_stats + rarer_log_p_stats + 
  plot_layout(ncol = 1, heights=c(1,2) ) +  
  plot_annotation(tag_levels = 'A')
split_p_stats
ggsave(paste0(output_dir, "/Fig_9_biotype_plot.pdf"), width=6.5, height=9, units="in")


# create some output text for paper
# create a file for text output
output_file <-file(paste0(output_dir, "alfa_analysis_text.txt"))

# pull out some p-values for text output
rRNA_oc_v_ff <- filter(df_wilcox, group1 == "Museum\nsamples" & group2 == "Fresh\nfrozen" & biotype == "rRNA") %>% pull(p.adj)
mRNA_oc_v_ff <- filter(df_wilcox, group1 == "Museum\nsamples" & group2 == "Fresh\nfrozen" & biotype == "protein coding") %>% pull(p.adj)
tRNA_oc_v_ff <- filter(df_wilcox, group1 == "Museum\nsamples" & group2 == "Fresh\nfrozen" & biotype == "tRNA") %>% pull(p.adj)
opp_strand_oc_v_ff <- filter(df_wilcox, group1 == "Museum\nsamples" & group2 == "Fresh\nfrozen" & biotype == "opposite strand") %>% pull(p.adj)

output_text <- ""
output_text <- paste0(
  output_text,
  "We also quantified the extent to which different types of non-ribosomal RNA survived in old samples. ",
  "We mapped reads to the D. melanogaster genome and used the ALFA software to quantify the different types ",
  "of RNA present in each dataset [PMID: 30922228]. ",
  "ALFA combines mapping information with genome annotation to assign mapped reads to one of a dozen RNA types (Fig. 9). ",
  "As expected, most reads were categorized as rRNA, though there was relatively less rRNA in older samples (Fig. 9A; p = ",
  sprintf("%0.1e", rRNA_oc_v_ff),
  " for museum vs. fresh samples). ",
  "There was also relatively less protein coding mRNA in older samples (Fig. 9B; p= ",
  sprintf("%0.1e", mRNA_oc_v_ff),
  "). ",
  "Opposite strand RNA levels – that is, reads derived from the RNA strand opposite an annotated feature - ",
  "were elevated in older samples, consistent with the preferential survival of antisense transcripts in RNA duplexes (p = ",
  sprintf("%0.1e", opp_strand_oc_v_ff),
  "). ",
  "Highly structured transfer RNA (tRNA) levels were also elevated in older samples ",
  "(p = ",
  sprintf("%0.1e", tRNA_oc_v_ff),
  ") [PMID: 3905254]. "
)
output_text

# Experimentally dried sample coverage levels were generally intermediate between fresh and museum samples, consistent with their intermediate age (Fig. X, XC).  

# output concatenated text to file
writeLines(output_text, output_file)

# close file 
close(output_file)
