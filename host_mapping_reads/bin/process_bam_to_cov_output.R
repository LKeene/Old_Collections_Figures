#!/usr/bin/env Rscript


# This code block sets up input arguments to either come from the command line
# (if running from the pipeline, so not interactively) or to use expected default values 
# (if running in interactive mode, for example in RStudio for troubleshooting
# or development).  
#
if  (!interactive()) {
  # if running from Rscript
  args = commandArgs(trailingOnly=TRUE)
  tsv_input                   = args[1]
  metadata_input              = args[2]
  interaction_categories_file = args[3]
  R_lib_dir                   = args[4]
  R_script_dir                = args[5]
  dataset_sizes_file          = args[6]
  output_dir                  = "./"
} else {
  # if running via RStudio
  tsv_input                   = "../results/collected_rRNA_locus_coverage.tsv"
  metadata_input              = "../results/collect/collected_metadata.csv"
  interaction_categories_file = "../refseq/interaction_categories.txt"
  R_lib_dir                   = "../lib/R/"
  R_script_dir                = "../../scripts/"
  output_dir                  = "../results/process/"
  dataset_sizes_file          =  "../results/all_read_counts.txt"
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

# create a file for text output
output_file <-file(paste0(output_dir, "rRNA_coverage_text.txt"))
output_text <- ""

# read in common color definitions
source(paste0(R_script_dir, "/plot_colors.R"))
fancy_color_scale <- c(fresh_frozen_color, experimental_dried_color, old_collection_color)

# read in metadata
metadata <- read.delim(metadata_input, sep=",", 
                       header=T, stringsAsFactors = F)

# reorder metadata sample types so they display in desired order
metadata$sample_type <- fct_relevel(metadata$sample_type, "Fresh_frozen", "Experimental_dried", "Old_Collection")

# reorder sample factors for display
metadata$sample_id <- fct_reorder(metadata$sample_id, metadata$date_collected, min)
metadata$sample_id_in_paper <- fct_reorder(metadata$sample_id_in_paper, metadata$date_collected, min)

# change category labels for plotting 
metadata$sample_type <- 
  recode(metadata$sample_type, 
         Old_Collection     = "Museum\nsamples", 
         Experimental_dried = "Experimental\ndried",
         Fresh_frozen       = "Fresh\nfrozen")

# read in dataset sizes file
dataset_sizes <- read.delim(dataset_sizes_file, sep="\t", header=T)
dataset_sizes <- dataset_sizes %>% filter(count_type == "post_trimming")
dataset_sizes <- dataset_sizes %>% mutate(total_reads = count) %>% select(-count_type, -count)

# read in rRNA coverage info
coverage <- read.delim(file = tsv_input, sep="\t", header=F)
# name columns
colnames(coverage) <- c("sample_id", "refseq",  
                      "position", "fwd_cov", "rev_cov")

# don't keep position 0 - undefined in 1-index sequences
coverage <- filter(coverage, position > 0)

# confirm that metadata exists for all datasets
dataset_names <- coverage %>% 
  group_by(sample_id) %>% summarize() %>% pull(sample_id)
missing_metadata <- !(dataset_names %in% metadata$sample_id)

if (any(missing_metadata)) {
  message("WARNING: missing metadata for the following datasets:")
  cat(dataset_names[missing_metadata])
}

# merge metadata into coverage table
coverage <- left_join(coverage, metadata, by="sample_id")

# read in interaction categories data: a file that has whether individual bases
# in the rRNA reference sequence are predicted to be involved in base pairing, pseudoknots, etc
interaction_categories <- read.delim(interaction_categories_file, sep="\t", header = F)
colnames(interaction_categories) <- c("refseq", "position", "interaction_category")

# relevel interaction categories
interaction_categories$interaction_category <- 
  fct_relevel(interaction_categories$interaction_category, 
              "base_paired", "not_paired", "pseudo_knot", "not_in_structure")

# what fraction of bases are in each category?
bp_in_refseq <- nrow(interaction_categories) 
category_fractions <- interaction_categories %>% group_by(interaction_category) %>% summarize (fraction = n() / bp_in_refseq,
                                                                                               n_in_cat = n())
category_fractions


# join in interaction categories to coverage
coverage <- left_join(coverage, interaction_categories, by=c("refseq", "position"))

# different types of strand-specific library prep and different types of adapters
# produce stranded RNA-seq libraries with different read1 orientations relative to
# the RNA strand of origin.  The R1 can be in the orientation of the original RNA
# or it can be the opposite.  See, for instance:
# https://www.idtdna.com/pages/support/faqs/can-the-xgen-unique-dual-index-umi-adapters-be-used-for-rna-seq
# We need to account for this possible difference.  To do this, the metadata
# spreadsheet has a R1_strand column which indicates the R1 orientation with respect
# to the original RNA..

# if R1_strand set as "forward", need to swap fwd & rev mapping reads because 
# of expectation of R1 orientation
coverage <- coverage %>% mutate(fwd_cov_temp = fwd_cov, 
                                fwd_cov = if_else(R1_strand == "forward", rev_cov, fwd_cov),
                                rev_cov = if_else(R1_strand == "forward", fwd_cov_temp, rev_cov)) %>%
  select(-fwd_cov_temp)

# calculate ratio of Fwd/Rev coverage
coverage <- coverage %>% mutate(cov_ratio = fwd_cov/rev_cov,
                                total_cov = fwd_cov + rev_cov)

# what is distribution of total average coverage?
coverage_averages <- 
  coverage %>% 
  group_by(sample_id) %>% 
  summarize(median_total_cov    = median(total_cov, na.rm = T),
            mean_total_cov      = mean  (total_cov, na.rm = T),
            sd_total_cov        = sd    (total_cov, na.rm = T),
            median_fwd_cov      = median(fwd_cov, na.rm = T),
            mean_fwd_cov        = mean  (fwd_cov, na.rm = T),
            sd_fwd_cov          = sd    (fwd_cov, na.rm = T),
            median_rev_cov      = median(rev_cov, na.rm = T),
            mean_rev_cov        = mean  (rev_cov, na.rm = T),
            sd_rev_cov          = sd    (rev_cov, na.rm = T))

# join in averages into main coverage df
coverage <- left_join(coverage, coverage_averages)

coverage_averages <- left_join(coverage_averages, metadata)

# filter out datasets for which no defined sample type metadata: neg ctrl samples
coverage_averages <- coverage_averages %>% filter (!is.na(sample_type))
            
# plot median coverage
coverage_levels <- ggplot(coverage_averages) + 
  geom_boxplot(aes(x=sample_type, y=median_total_cov, color=sample_type)) +
    theme_bw() +
  scale_y_log10() +
  scale_color_manual(values = fancy_color_scale) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

coverage_levels
ggsave(paste0(output_dir, "median_coverage_levels.pdf"), width=10, height=7.5, units="in" )

# filter out samples with < 10x median coverage
min_sufficient_coverage <- 10
samples_with_sufficient_coverage <- coverage_averages %>% filter(median_total_cov > min_sufficient_coverage) %>% pull(sample_id)
coverage <- filter(coverage, sample_id %in% samples_with_sufficient_coverage)
coverage_averages <- filter(coverage_averages, sample_id %in% samples_with_sufficient_coverage)

# filter out non-melanogaster datasets
melanogaster_samples <- metadata %>% filter(str_detect(species, "melano")) %>% pull(sample_id)
coverage <- filter(coverage, sample_id %in% melanogaster_samples)
coverage_averages <- filter(coverage_averages, sample_id %in% melanogaster_samples)

# calculate coverage at each position realtive to average coverage across whole refseq
coverage <- coverage %>%
  mutate(relative_total_cov = total_cov / median_total_cov,
         relative_fwd_cov   = fwd_cov   / median_fwd_cov,
         relative_rev_cov   = rev_cov   / median_rev_cov)
# calculate coverage values normalized to total read counts
# an RPM-like normalization
coverage <- left_join(coverage, dataset_sizes)
coverage <- coverage %>%
  mutate(total_cov_rpm = total_cov * 1e6 / total_reads,
         fwd_cov_rpm   = fwd_cov   * 1e6 / total_reads,
         rev_cov_rpm   = rev_cov   * 1e6 / total_reads)

# create coverage windows 
window_size <- 20
coverage <- coverage %>% 
  mutate(window_start = window_size * (position %/% window_size),
         window_mid   = window_start + (window_size / 2))

# calculate % base paired in each window
coverage_windows <- coverage %>% group_by(sample_id, refseq, window_mid) %>% 
  summarize(mean_fwd_cov = mean(fwd_cov), 
            mean_rev_cov = mean(rev_cov),
            mean_relative_fwd_cov = mean(relative_fwd_cov),
            mean_relative_rev_cov = mean(relative_rev_cov),
            mean_fwd_cov_rpm      = mean(fwd_cov_rpm),
            mean_rev_cov_rpm      = mean(rev_cov_rpm),
            .groups="drop")

# merge metadata into counts table
coverage_windows <- left_join(coverage_windows, metadata, by="sample_id")


# -------------------------------
# RATIOS OF + to - strand mapping
# -------------------------------

ratio_averages <- 
  filter(coverage, is.finite(cov_ratio)) %>%
  group_by(sample_id, refseq) %>%
  summarize(median_cov_ratio = median(cov_ratio),
            mean_cov_ratio   = mean(cov_ratio, na.rm = T),
            sd_cov_ratio     = sd(cov_ratio), 
            .groups = "drop")

ratio_averages <- left_join(ratio_averages, metadata)

# run stats 

# normally distributed?
ratio_shapiro <- ratio_averages %>% group_by(sample_type) %>% shapiro_test(median_cov_ratio)
filter(ratio_shapiro, p<0.05)

# not normally distributed so use non-parametric test 
df_ratio_wilcox      <- ratio_averages %>% wilcox_test(median_cov_ratio ~ sample_type)
df_ratio_wilcox_plot <- df_ratio_wilcox %>% add_xy_position(x = "sample_type")

y_positions <- c(3.5,3.7,3.5)

# make a combined strand ratio plot for supp figures
rRNA_ratio_p <- ggplot(ratio_averages) +
  geom_boxplot(aes(x=sample_type, y=median_cov_ratio, color=sample_type), fill=NA, size=0.5, outlier.shape=NA) +
  geom_jitter(aes(x=sample_type, y=median_cov_ratio, fill=sample_type), shape=21, color="black", stroke=0.2, size=2, height=0, width=0.25) +
  scale_y_log10() +
  theme_this_paper(base_size = 14) + 
  scale_color_manual(values = fancy_color_scale) +
  scale_fill_manual(values = fancy_color_scale) +
  # ggtitle("Surviving ribosomal RNA is preferentially anti-sense", subtitle="consistent with the hypothesis that old RNA is double-stranded") +
  theme(plot.title = element_text(size = 14),
        plot.subtitle = element_text(size = 13)) +
  # ylab("Median ratio of +strand to -strand\nrRNA-mapping reads\nin individual datasets") +
  ylab("Median ratio of +strand \nto -strand\nrRNA-mapping reads") + 
  xlab("") +
  stat_pvalue_manual(df_ratio_wilcox_plot, y.position = y_positions, 
                     tip.length = 0.01, bracket.shorten=0.1, size=3, bracket.nudget.y=0.1)


rRNA_ratio_p
ggsave(paste0(output_dir, "Fig_6_rRNA_strandedness_plot.pdf"), width=7, height=4, units="in")

# normally distributed?
fwd_cov_shapiro <- coverage_averages %>% group_by(sample_type) %>% shapiro_test(median_fwd_cov)
filter(fwd_cov_shapiro, p<0.05)

# not normally distributed so use non-parametric test 
df_fwd_cov_wilcox      <- coverage_averages %>% wilcox_test(median_fwd_cov ~ sample_type)
df_fwd_cov_wilcox_plot <- df_fwd_cov_wilcox %>% add_xy_position(x = "sample_type")

y_positions <- c(6.0,6.2,6.0)

rRNA_fwd_p <- ggplot(coverage_averages) +
  geom_boxplot(aes(x=sample_type, y=median_fwd_cov, color=sample_type), fill=NA, size=0.5, outlier.shape=NA) +
  geom_jitter (aes(x=sample_type, y=median_fwd_cov, fill=sample_type), shape=21, color="black", stroke=0.2, size=2, height=0, width=0.25) +
  scale_y_log10() +
  theme_this_paper(base_size = 14) +
  scale_color_manual(values = fancy_color_scale) +
  scale_fill_manual(values = fancy_color_scale) +
  # ggtitle("Similar +sense rRNA coverage levels in older samples") +
  ylab("Median +strand rRNA\ncoverage depth") + 
  xlab("") +
  stat_pvalue_manual(df_fwd_cov_wilcox_plot, y.position = y_positions, 
                     tip.length = 0.01, bracket.shorten=0.1, size=3, bracket.nudget.y=0.1)
  
rRNA_fwd_p

# normally distributed?
rev_cov_shapiro <- coverage_averages %>% group_by(sample_type) %>% shapiro_test(median_rev_cov)
filter(rev_cov_shapiro, p<0.05)

# not normally distributed so use non-parametric test 
df_rev_cov_wilcox      <- coverage_averages %>% wilcox_test(median_rev_cov ~ sample_type)
df_rev_cov_wilcox_plot <- df_rev_cov_wilcox %>% add_xy_position(x = "sample_type")

y_positions <- c(4.0,4.2,4.0)

rRNA_rev_p <- ggplot(coverage_averages) +
  geom_boxplot(aes(x=sample_type, y=median_rev_cov, color=sample_type), fill=NA, size=0.5, outlier.shape=NA) +
  geom_jitter (aes(x=sample_type, y=median_rev_cov, fill=sample_type), shape=21, color="black", stroke=0.2, size=2, height=0, width=0.25) +
  scale_y_log10() +
  theme_this_paper(base_size = 14) +
  scale_color_manual(values = fancy_color_scale) +
  scale_fill_manual(values = fancy_color_scale) +
  # ggtitle("There is relatively more antisense ribosomal RNA in older samples") + 
  ylab("Median -strand rRNA\ncoverage depth") + 
  xlab("") +
  stat_pvalue_manual(df_fwd_cov_wilcox_plot, y.position = y_positions, 
                     tip.length = 0.01, bracket.shorten=0.1, size=3, bracket.nudget.y=0.1)

rRNA_rev_p

# combined plot
combined_rRNA_plot <- 
  # apply_plot_theme(rRNA_ratio_p) + 
  # apply_plot_theme(rRNA_fwd_p) + 
  # apply_plot_theme(rRNA_rev_p) + 
  rRNA_ratio_p + 
  rRNA_fwd_p + 
  rRNA_rev_p + 
  plot_layout(ncol = 1) + 
  plot_annotation(tag_levels = 'A')
combined_rRNA_plot 
ggsave(paste0(output_dir, "Fig_SX_rRNA_strandedness_plot_all.pdf"), width=7, height=9, units="in")

# output text for paper

ratio_averages_by_type <- ratio_averages %>% group_by(sample_type) %>% summarize(median_ratio = median(median_cov_ratio))
ff_ratio_avg <- ratio_averages_by_type %>% filter(sample_type == "Fresh\nfrozen") %>% pull(median_ratio)
ed_ratio_avg <- ratio_averages_by_type %>% filter(sample_type == "Experimental\ndried") %>% pull(median_ratio)
oc_ratio_avg <- ratio_averages_by_type %>% filter(sample_type == "Museum\nsamples") %>% pull(median_ratio)

ed_ratio_pval <- df_ratio_wilcox %>% filter(group1 == "Experimental\ndried" & group2 == "Fresh\nfrozen") %>% pull(p.adj)
oc_ratio_pval <- df_ratio_wilcox %>% filter(group1 == "Museum\nsamples" & group2 == "Fresh\nfrozen") %>% pull(p.adj)

output_text <- paste0(
  output_text,
  "There was evidence that double-stranded rRNA preferentially survived in old samples. ",
  "In fresh frozen samples, most rRNA was +strand: the ratio of +strand to -strand RNA was ",
  sprintf("%0.0f", ff_ratio_avg), ":1. ",
  "In dried and museum sample, the +strand:-strand ratio dropped to ",
  sprintf("%0.0f", ed_ratio_avg), ":1 and ", 
  sprintf("%0.0f", oc_ratio_avg), ":1 (p = ",
  sprintf("%0.1e", ed_ratio_pval), " and ",
  sprintf("%0.1e", oc_ratio_pval), " respectively; Fig. XA). ",
  "This decreased ratio was driven by a relative decrease in +strand RNA and a relative ",
  "increase in -strand RNA in older samples (Figs. XB & XC & Supp Fig X). ",
  "negative-sense rRNA could derive from anti-sense transcription of rRNA genes or psuedogenes. ",
  "Indeed, oppsite-strand transcripts levels were elevated in old samples (Fig. X), ",
  "consistent with preferential survival of sense:antisense duplexes. ",
  "\n\n"
)

writeLines(output_text, output_file)

# an example of -strand rRNA transcription:
# http://flybase.org/jbrowse/?data=data%2Fjson%2Fdmel&loc=X%3A23276863..23284877&tracks=Gene_span%2Ctopoview_NCBI_aggregate_100&highlight=

# ---------------------------------------------------
# Coverage by base type (base-paired, unpaired, etc)
# ---------------------------------------------------

category_coverage_averages <- coverage %>% 
  group_by(sample_id, refseq, interaction_category) %>%
  summarize(median_category_fwd_cov = median(fwd_cov), .groups="drop")

sample_coverage_averages <- coverage %>% 
  group_by(sample_id, refseq) %>%
  summarize(median_refseq_fwd_cov = median(fwd_cov), .groups="drop")

category_coverage_averages <- left_join(category_coverage_averages, sample_coverage_averages)

category_coverage_averages <- mutate(category_coverage_averages, relative_category_fwd_cov = median_category_fwd_cov / median_refseq_fwd_cov)

category_coverage_averages <- left_join(category_coverage_averages, metadata)

# total median coverage for each category / sample-type
coverage_averages_by_category <- category_coverage_averages %>%
  group_by(sample_type, interaction_category) %>%
  summarize(median_relative_category_fwd_cov = median(relative_category_fwd_cov),
            .groups="drop")


# change category labels for plotting 
category_coverage_averages$interaction_category <- 
  recode(category_coverage_averages$interaction_category, 
         not_in_structure   = "Bases not in 3D structure",
         base_paired        = "Paired bases",
         not_paired         = "Unpaired bases",
         pseudo_knot        = "Bases in pseudo-knot\ntype structure")

# do stats

# is normal?
category_shapiro <- category_coverage_averages %>% 
  group_by(sample_type, interaction_category) %>% shapiro_test(relative_category_fwd_cov)
filter(category_shapiro, p < 0.05)

# not normally distributed so do Wilcoxon
df_cat_wilcox      <- category_coverage_averages %>% group_by(interaction_category) %>% wilcox_test(relative_category_fwd_cov ~ sample_type)
df_cat_wilcox_plot <- df_cat_wilcox %>% add_xy_position(x = "sample_type")

y_positions <- rep(c(4.05,4.2,4.05), 4)
       
ggplot(category_coverage_averages) +
  geom_boxplot(aes(x=sample_type, y=relative_category_fwd_cov, color=sample_type), 
               fill=NA, size=0.3, outlier.shape=NA) +
  geom_jitter (aes(x=sample_type, y=relative_category_fwd_cov, fill =sample_type), 
               shape=21, color="black", stroke=0.2, size=2, height=0, width=0.25, alpha=0.5) +
  geom_hline  (aes(yintercept = 1), size = 0.25, linetype=2, color="grey50") +
  facet_wrap(~interaction_category, nrow=1) +
  scale_color_manual(values = fancy_color_scale) +
  scale_fill_manual(values = fancy_color_scale) +
  ylab("Median coverage of bases relative to overall median coverage") + 
  xlab("") +
  theme_this_paper(base_size = 12) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  stat_pvalue_manual(df_cat_wilcox_plot, y.position = y_positions, 
                     # label="p.adj",
                     tip.length = 0.01, bracket.shorten=0.1, size=3, bracket.nudget.y=0.1)


ggsave(paste0(output_dir, "Fig_8C_ribosomal_coverage_by_base_categories.pdf"), width=7.5, height=5.5, units="in")

# text for paper
output_text <- paste0(
  output_text,
  "In addition to different ratios of +strand to -strand rRNA, rRNA coverage patterns varied in old samples (Fig. X). ",
  "rRNA coverage in fresh sample datasets was relatively even. In contrast coverage levels in older samples flucuated across rRNAs, ",
  " with some regions exhibiting reproducibly lower than average coverage. ",
  "For instance, coverage around base 3000 of the 28S rRNA had consistently lower coverage in most old samples (Fig. X). ",
  "\n\n",
  "To investigate why some rRNA regions might be differentially represented in old samples ",
  "we took advantage of a 3-dimensional structure of the D. melanogaster ribosome [PMID: 23636399]. ",
  "This structure includes ribosomal proteins and RNA (Fig. XA). ",
  "It is possible to see whether individual RNA bases in the structure ",
  "are interacting with other bases, and whether individual bases are actually present in the structure (Fig. XB). ",
  "\n\n"
  )

# text for paper
output_text <- paste0(
  output_text,
  "Using RNApdbee software, we categorized rRNA bases into one of 4 categories: ",
  "base-paired, ",  
  "unpaired, present in a higher-order secondary structure like pseudoknots, or missing from the 3D structure ",
  "[PMID:  15941360, 29718468, 23636399]. ",
  sprintf ("%0.1f",
           category_fractions %>% filter(interaction_category == "base_paired") %>% pull(fraction) * 100) ,
  "% of rRNA bases were categorized as paired with another rRNA base, ",
  sprintf ("%0.1f",
           category_fractions %>% filter(interaction_category == "not_paired") %>% pull(fraction) * 100) ,
  "% as unpaired, ",
  sprintf ("%0.1f",
           category_fractions %>% filter(interaction_category == "pseudo_knot") %>% pull(fraction) * 100) ,
  "% in higher order secondary structures, and ",
  sprintf ("%0.1f",
           category_fractions %>% filter(interaction_category == "not_in_structure") %>% pull(fraction) * 100) ,
  "% were not present in the 3D structure (Fig. X)." ,
  "\n\n"
)

# pull out certain p values
oc_v_f_bp_p <- df_cat_wilcox %>% 
  filter(interaction_category == "Paired bases" & group1 == "Museum\nsamples" & group2 == "Fresh\nfrozen") %>%
  pull(p.adj)

oc_v_f_up_p <- df_cat_wilcox %>% 
  filter(interaction_category == "Unpaired bases" & group1 == "Museum\nsamples" & group2 == "Fresh\nfrozen") %>%
  pull(p.adj)

oc_v_f_pk_p <- df_cat_wilcox %>% 
  filter(interaction_category == "Bases in pseudo knot-type structures" & group1 == "Museum\nsamples" & group2 == "Fresh\nfrozen") %>%
  pull(p.adj)

oc_v_f_np_p <- df_cat_wilcox %>% 
  filter(interaction_category == "Bases not in 3D structure" & group1 == "Museum\nsamples" & group2 == "Fresh\nfrozen") %>%
  pull(p.adj)

oc_v_f_np_fd <- coverage_averages_by_category %>% 
  filter(interaction_category == "not_in_structure" & sample_type == "Museum\nsamples") %>%
  pull(median_relative_category_fwd_cov) 


output_text <- paste0(
  output_text,
 "We calculated the median coverage of bases in each of these 4 categories ",
 "relative to total median coverage (Fig X). ",
 "Paired bases had higher average coverage in older samples (p = ",
 sprintf("%0.1e", oc_v_f_bp_p),
 " in museum samples vs. fresh ones). ",
 "Similarly, bases in higher-order secondary structures had higher coverage in old samples (p = ",
 sprintf("%0.1e ", oc_v_f_pk_p),
 "relative to fresh). ",
 "In contrast, unpaired bases had lower average coverage in old samples (p = ",
 sprintf("%0.1e", oc_v_f_up_p),
 "vs. fresh samples). ",
 "Bases that were not captured in the 3D structure were the least well reprented in older samples, ",
 " with coverage levels ",
 sprintf("%0.0f", 1/oc_v_f_np_fd),
 "x lower than overall average coverage levels in museum samples ",
 "(p = ",
 sprintf("%0.1e ", oc_v_f_np_p),
 "relative to Fresh frozen). ",
 "Experimentally dried sample coverage levels were generally intermediate between ",
 "fresh and museum samples, consistent with their intermediate age (Fig. X). ",
 "\n\n"
 )
 
# new paragraph summarizing
output_text <- paste0(
  output_text,
 "These patterns suggested that the molecular environment immediately surrounding RNA bases influenced the likelihood that they ",
 "would survive over long periods.  Being base-paired or in a higher-order RNA secondary structure protected bases, ",
 "enabling them to survive over long periods. ",
 "In contrast, unpaired bases or bases not present in the 3D structure ",
 "- presumably because they were outside of the protective environment of the ribosome - ",
 "were less likely to surivive. ",
 "\n\n"
)

# Figure legend
output_text <- paste0(
  output_text,
  "Figure X: Certain types of bases tended to survive in old rRNA. ",
  "The median coverage level of bases in the indicated categories relative to total median coverage is plotted. ",
  "Each point represents a dataset from an individual fly. ",
  "Signifcance levels of Wilcoxon test adjusted pvalues are indicated. ",
  "\n\n"
)
       
       
# Methods text:   What do significance asterisks mean?  
# from: http://rpkgs.datanovia.com/ggpubr/reference/stat_compare_means.html
output_text <- paste0(
  output_text,
  "P-values adjusted for multiple testing are shown on plots using the following significance levels:  ",
  "ns: p > 0.05; ",
  "*: p <= 0.05; ",
  "**: p <= 0.01; ",
  "***: p <= 0.001; ",
  "****: p <= 0.0001; ",
  "\n\n"
)  

## ggpubr citation:
output_text <- paste0(output_text,
       "Kassambara A (2023). ggpubr: 'ggplot2' Based Publication Ready Plots. R package version 0.6.0, https://rpkgs.datanovia.com/ggpubr/.\n\n")

print(output_text)
writeLines(output_text, output_file)

# ---------------------------------
# Plot coverage by position
# ---------------------------------

# relabel refseq names
# A2/B5 refer to 3D structure chains in PDB
coverage_windows$refseq    <- recode(coverage_windows$refseq, B2 = "18S SSU rRNA", A5 = "28S LSU rRNA")

# reorder sample factors for display
coverage_windows$sample_id <- fct_reorder(coverage_windows$sample_id, coverage_windows$date_collected, min)

# plot coverage with individual datasets as individual lines
ggplot(coverage_windows) + 
  geom_line(aes(x=window_mid, y=mean_fwd_cov_rpm, group=sample_id, color=sample_type), linewidth=0.25) +
  scale_color_manual(values = fancy_color_scale) +
  scale_y_log10() +
  facet_grid(sample_type~refseq, scales="free_x") +
  theme_this_paper(base_size = 14)  +
  theme(strip.text.y = element_text(angle = 0))+
  theme(legend.position = "none") +
  ylab(paste0("Mean +strand coverage per million reads")) +
  xlab("Position in rRNA") 

ggsave(paste0(output_dir, "Fig_7_rRNA_coverage_traces.pdf"), width=10, height=7.5, units="in")

# plot -strand coverage with individual datasets as individual lines
ggplot(coverage_windows) + 
  geom_line(aes(x=window_mid, y=mean_rev_cov_rpm, group=sample_id, color=sample_type), linewidth=0.25) +
  scale_color_manual(values = fancy_color_scale) +
  scale_y_log10() +
  facet_grid(sample_type~refseq, scales="free_x") +
  theme_this_paper(base_size = 14)  +
  theme(strip.text.y = element_text(angle = 0))+
  theme(legend.position = "none") +
  ylab(paste0("Mean -strand coverage per million reads")) +
  xlab("Position in rRNA") 

ggsave(paste0(output_dir, "Fig_SX_neg_strand_rRNA_coverage.pdf"), width=7.5, height=7.5, units="in")

# plot +strand vs. -strand coverage as scatter plot
ggplot(coverage_windows) + 
  geom_point(aes(x=mean_fwd_cov_rpm, y=mean_rev_cov_rpm, fill=sample_type), shape=21, color="black", size=1, stroke=0.25, alpha=0.25) +
  geom_abline(slope=1, intercept=0, color="slateblue", linewidth=0.5, linetype=2, alpha=0.5) +
  scale_color_manual(values = fancy_color_scale) +
  scale_fill_manual(values = fancy_color_scale) +
  scale_x_log10() +
  scale_y_log10() +
  coord_equal()+
  facet_grid(sample_type~refseq) +
  theme_this_paper(base_size = 14)  +
  theme(strip.text.y = element_text(angle = 0))+
  theme(legend.position = "none") +
  ylab(paste0("Mean -strand coverage per million reads")) +
  ylab(paste0("Mean +strand coverage per million reads")) 

ggsave(paste0(output_dir, "Fig_SX_both_strand_rRNA_coverage_scatter_plot.pdf"), width=7.5, height=7.5, units="in")

# is +strand and -strand coverage correlated?
cov_cor <- coverage %>% group_by(sample_id, refseq) %>% cor_test(fwd_cov, rev_cov)
cov_cor <- left_join(cov_cor, metadata)

# how well correlated are +strand and -strand
ggplot(cov_cor) +
  geom_boxplot(aes(x=sample_type, y=cor, color=sample_type),
               fill=NA, size=0.3) +
  scale_color_manual(values = fancy_color_scale) +
  facet_wrap(~refseq, nrow = 1) + 
  theme_this_paper(base_size = 14)  +
  theme(strip.text.y = element_text(angle = 0))+
  theme(legend.position = "none") +
  ylab("Correlation coefficient of +strand coverage to -strand coverage") +
  xlab("") 

ggsave(paste0(output_dir, "Fig_SX_both_strand_rRNA_coverage_correlation.pdf"), width=7.5, height=7.5, units="in")

# a function to plot coverage by position for one or more datasets
# with datasets faceted
plot_coverage_by_position <- function (df_to_plot, facet_var="sample_id_in_paper"){
  ggplot(df_to_plot) + 
    geom_line(aes(x=window_mid, y=mean_fwd_cov_rpm, color=sample_type)) +
    scale_color_manual(values = fancy_color_scale) +
    scale_y_log10() +
    facet_grid( vars(.data[[facet_var]]), vars(refseq), scales="free_x") +
    theme_this_paper(base_size = 14)  +
    theme(strip.text.y = element_text(angle = 0))+
    theme(legend.position = "none") +
    ylab("+strand rRNA coverage per million reads") +
    xlab("Position in rRNA") 
}

plot_coverage_by_position(filter(coverage_windows, sample_type == "Museum\nsamples"))
ggsave(paste0(output_dir, "Fig_SX_OC_sample_coverage.pdf"), width=7.5, height=7.5, units="in")

plot_coverage_by_position(filter(coverage_windows, sample_type == "Experimental\ndried"), facet_var="sample_id")
ggsave(paste0(output_dir, "Fig_SX_ED_sample_coverage.pdf"), width=7.5, height=7.5, units="in")

plot_coverage_by_position(filter(coverage_windows, sample_type == "Fresh\nfrozen"), facet_var="sample_id")
ggsave(paste0(output_dir, "Fig_SX_FF_sample_coverage.pdf"), width=7.5, height=7.5, units="in")


# close file 
close(output_file)
