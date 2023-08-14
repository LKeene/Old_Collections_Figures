#!/usr/bin/env Rscript


# This code block sets up input arguments to either come from the command line
# (if running from the pipeline, so not interactively) or to use expected default values 
# (if running in interactive mode, for example in RStudio for troubleshooting
# or development).  
#
if (!interactive()) {
  # if running from Rscript
  args = commandArgs(trailingOnly=TRUE)
  tsv_input                   = args[1]
  metadata_input              = args[2]
  interaction_categories_file = args[3]
} else {
  # if running via RStudio
  tsv_input                   = "../results/collected_rRNA_locus_coverage.tsv"
  metadata_input              = "../refseq/metadata.csv"
  interaction_categories_file = "../refseq/interaction_categories.txt"
}

# quit early while developing
if (!interactive()) {
  quit(status = 0)
}

library(tidyverse)
library(patchwork)

# read in metadata
metadata <- read.delim(metadata_input, sep=",", 
                       header=T, stringsAsFactors = F)

# reorder metadata sample types so they display in desired order
metadata$sample_type <- fct_relevel(metadata$sample_type, "Old_Collection", "Experimental_dried", "Fresh_frozen")

coverage <- read.delim(file = tsv_input, sep="\t", header=F)
coverage_columns <- c("sample_id", "refseq", 
                      "position", "fwd_cov", "rev_cov")
colnames(coverage) <- coverage_columns

# confirm that metadata exists for all datasets
dataset_names <- coverage %>% 
  group_by(sample_id) %>% summarize() %>% pull(sample_id)
missing_metadata <- !(dataset_names %in% metadata$sample_id)

if (any(missing_metadata)) {
  message("ERROR: missing metadata for the following datasets:")
  cat(dataset_names[missing_metadata])
}

# merge metadata into counts table
coverage <- left_join(coverage, metadata, by="sample_id")

# read in interaction categories data: a file that has whether individual bases
# in the rRNA reference sequence are predicted to be involved in base pairing, pseudoknots, etc
interaction_categories <- read.delim(interaction_categories_file, sep="\t", header = F)
colnames(interaction_categories) <- c("refseq", "position", "interaction_category")
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

coverage_averages <- left_join(coverage_averages, metadata)
            
# plot median coverage
coverage_levels <- ggplot(coverage_averages) + 
  geom_col(aes(x=sample_id, y=median_total_cov)) +

    theme_bw() +
  scale_y_log10() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

coverage_levels
ggsave("median_coverage_levels.pdf", width=10, height=7.5, units="in" )

# filter out samples with < 10x median coverage
min_sufficient_coverage <- 10
samples_with_sufficient_coverage <- coverage_averages %>% filter(median_total_cov > min_sufficient_coverage) %>% pull(sample_id)
coverage <- filter(coverage, sample_id %in% samples_with_sufficient_coverage)
coverage_averages <- filter(coverage_averages, sample_id %in% samples_with_sufficient_coverage)

# filter out D. simulans samples
non_simulans_samples <- metadata %>% filter(species != "simulans") %>% pull(sample_id)
coverage <- filter(coverage, sample_id %in% non_simulans_samples)
coverage_averages <- filter(coverage_averages, sample_id %in% non_simulans_samples)

# reorder sample factors for display
coverage$sample_id <- fct_reorder(coverage$sample_id, coverage$date_collected, min)


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

#RColorBrewer Dark2 colors
dark2_lightblue   <- rgb(166/255,206/255,227/255)
dark2_blue        <- rgb(31/255,120/255,180/255)
dark2_lightgreen  <- rgb(178/255,223/255,138/255)
dark2_green       <- rgb(51/255,160/255,44/255)
dark2_pink        <- rgb(251/255,154/255,153/255)
dark2_red         <- rgb(227/255,26/255,28/255)
dark2_lightorange <- rgb(253/255,191/255,111/255)
dark2_orange      <- rgb(255/255,127/255,0/255)
dark2_lightpurple <- rgb(202/255,178/255,214/255)
dark2_purple      <- rgb(106/255,61/255,154/255)
dark2_yellow      <- rgb(255/255,255/255,153/255)
dark2_brown       <- rgb(177/255,89/255,40/255)

fancy_color_scale <- c(dark2_orange, dark2_blue, dark2_green)

# change category labels for plotting 
ratio_averages$sample_type <- 
  recode(ratio_averages$sample_type, 
         Old_Collection     = "Old\ncollections", 
         Experimental_dried = "Experimental\ndried",
         Fresh_frozen       = "Fresh\nfrozen")

rRNA_ratio_p <- ggplot(ratio_averages) +
  geom_boxplot(aes(x=sample_type, y=median_cov_ratio, color=sample_type), fill=NA, size=0.5, outlier.shape=NA) +
  geom_jitter(aes(x=sample_type, y=median_cov_ratio, fill=sample_type), shape=21, color="black", stroke=0.2, size=2, height=0, width=0.25) +
  scale_y_log10() +
  theme_classic(base_size = 13) +
  theme(legend.position = "none") +
  scale_color_manual(values = fancy_color_scale) +
  scale_fill_manual(values = fancy_color_scale) +
  # ggtitle("Surviving ribosomal RNA is preferentially anti-sense", subtitle="consistent with the hypothesis that old RNA is double-stranded") +
  theme(plot.title = element_text(size = 14),
        plot.subtitle = element_text(size = 13)) +
  ylab("Median ratio of +strand to -strand\nrRNA-mapping reads\nin individual datasets") +
  xlab("")

rRNA_ratio_p
ggsave("rRNA_ratios_figure.pdf", width=5, height=6.5, units="in")

# change category labels for plotting 
coverage_averages$sample_type <- 
  recode(coverage_averages$sample_type, 
         Old_Collection     = "Old\ncollections", 
         Experimental_dried = "Experimental\ndried",
         Fresh_frozen       = "Fresh\nfrozen")

  
rRNA_fwd_p <- ggplot(coverage_averages) +
  geom_boxplot(aes(x=sample_type, y=median_fwd_cov, color=sample_type), fill=NA, size=0.5, outlier.shape=NA) +
  geom_jitter (aes(x=sample_type, y=median_fwd_cov, fill=sample_type), shape=21, color="black", stroke=0.2, size=2, height=0, width=0.25) +
  scale_y_log10() +
  theme_classic(base_size = 13) +
  theme(legend.position = "none") +
  scale_color_manual(values = fancy_color_scale) +
  scale_fill_manual(values = fancy_color_scale) +
  # ggtitle("Similar +sense rRNA coverage levels in older samples") +
  ylab("Median +strand\ncoverage depth") + 
  xlab("")
  
rRNA_fwd_p
ggsave("rRNA_fwd_coverage_figure.pdf", width=5, height=6.5, units="in")

rRNA_rev_p <- ggplot(coverage_averages) +
  geom_boxplot(aes(x=sample_type, y=median_rev_cov, color=sample_type), fill=NA, size=0.5, outlier.shape=NA) +
  geom_jitter (aes(x=sample_type, y=median_rev_cov, fill=sample_type), shape=21, color="black", stroke=0.2, size=2, height=0, width=0.25) +
  scale_y_log10() +
  theme_classic(base_size = 13) +
  theme(legend.position = "none") +
  scale_color_manual(values = fancy_color_scale) +
  scale_fill_manual(values = fancy_color_scale) +
  # ggtitle("There is relatively more antisense ribosomal RNA in older samples") + 
  ylab("Median -strand\ncoverage depth") + 
  xlab("")

rRNA_rev_p
ggsave("rRNA_rev_coverage_figure.pdf", width=5, height=6.5, units="in")

# combined plot
combined_rRNA_plot <- rRNA_ratio_p + rRNA_fwd_p + rRNA_rev_p + plot_layout(ncol = 1)
combined_rRNA_plot
ggsave("combined_rRNA_plot.pdf", width=4, height=10, units="in")

# The drosophila genome contains antisense rRNA pseudogenes


# -----------------------
# INTERACTION CATEGORIES 
# -----------------------

# the oldest sample
oldest <- filter(coverage, sample_id =="1004284")
# a +control 
control <- filter(coverage, sample_id =="PosCtrl_Pool1")

coverage_long <- coverage %>% pivot_longer(cols = c(total_cov, fwd_cov, rev_cov), names_to = "cov_type", names_pattern = "(.*)_.*", values_to="cov")

# calculate coverage averages, but break down by RNApdbee interaction category
coverage_averages_by_category <- 
  coverage_long %>% 
  group_by(sample_id, interaction_category, cov_type) %>% 
  summarize(median_cov    = median(cov, na.rm = T),
            .groups="drop")

coverage_averages_by_category_wide <- 
  coverage_averages_by_category %>% pivot_wider(names_from = interaction_category, values_from=median_cov)

coverage_averages_by_category_ratios <-
  coverage_averages_by_category_wide %>% 
  mutate(
    not_in_structure = not_in_structure / base_paired,
    not_paired       = not_paired       / base_paired,
    pseudo_knot      = pseudo_knot      / base_paired)  %>% 
  pivot_longer(cols=c(not_in_structure, not_paired, pseudo_knot), names_to = "interaction_category", values_to = "rel_base_paired")



coverage_averages_by_category_ratios <- left_join(coverage_averages_by_category_ratios, metadata)

# change category labels for plotting 
coverage_averages_by_category_ratios$sample_type <- 
  recode(coverage_averages_by_category_ratios$sample_type, 
         Old_Collection     = "Old\ncollections", 
         Experimental_dried = "Experimental\ndried",
         Fresh_frozen       = "Fresh\nfrozen")
       
# change category labels for plotting 
coverage_averages_by_category_ratios$interaction_category <- 
  recode(coverage_averages_by_category_ratios$interaction_category, 
         not_in_structure   = "Bases not in 3D structure",
         not_paired         = "Bases in structure but not paired",
         pseudo_knot        = "Bases in pseudo knot-type structures")
       
coverage_averages_by_category_ratios$interaction_category <- 
  fct_relevel(coverage_averages_by_category_ratios$interaction_category, 
              "Bases not in 3D structure",
              "Bases in structure but not paired",
              "Bases in pseudo knot-type structures")

ggplot(filter(coverage_averages_by_category_ratios, cov_type == "total")) +
  geom_boxplot(aes(x=sample_type, y=rel_base_paired, color=sample_type), fill=NA, size=0.5, outlier.shape=NA) +
  geom_jitter (aes(x=sample_type, y=rel_base_paired, fill =sample_type), shape=21, color="black", stroke=0.2, size=3, height=0, width=0.25) +
  geom_hline  (aes(yintercept = 1), size = 0.25, linetype=2, color="grey50") +
  # scale_y_log10() +
  theme_classic(base_size = 13) +
  theme(legend.position = "none") +
  scale_color_manual(values = fancy_color_scale) +
  scale_fill_manual(values = fancy_color_scale) +
  # ggtitle("There is relatively more antisense ribosomal RNA in older samples") + 
  # facet_grid(cov_type ~ interaction_category) +
  facet_wrap(~interaction_category) +
  ylab("Median coverage levels relative to coverage of bases\nthat are base-paired in ribosome cryo-EM structure") +
  xlab("")

ggsave("coverage_by_base_categories.pdf", width=10, height=6.5, units="in")


# interaction categories in the oldest sample
ggplot() +
  geom_boxplot(aes(x=interaction_category, y=fwd_coverage)) +
  scale_y_log10() + 
  theme_bw()

# interaction categories in the oldest sample
ggplot(filter(coverage, !is.na(interaction_category))) +
  geom_boxplot(aes(x=interaction_category, y=fwd_coverage)) +
  facet_wrap(~sample_id) +
  scale_y_log10() + 
  theme_bw()

# interaction categories in the oldest sample
ggplot(filter(coverage, !is.na(interaction_category))) +
  geom_boxplot(aes(x=interaction_category, y=rev_coverage)) +
  facet_wrap(~sample_id) +
  scale_y_log10() + 
  theme_bw()

# create coverage windows 
window_size <- 20
coverage <- coverage %>% mutate(window_start = window_size * (position %/% window_size),
                         window_mid   = window_start + (window_size / 2))

# calculate % base paired in each window
coverage_windows <- coverage %>% group_by(sample_id, refseq, window_mid) %>% 
  summarize(mean_fwd_cov = mean(fwd_coverage), .groups="drop")

# merge metadata into counts table
coverage_windows <- left_join(coverage_windows, metadata, by="sample_id")


bp <- bp %>% mutate(window_start = window_size * (position %/% window_size),
                    window_mid   = window_start + (window_size / 2))

# calculate % base paired in each window
bp_window_values <- bp %>% group_by(window_mid) %>% summarize(pct_bp = sum(in_bp)/window_size)

# flesh out values in bp_windows: need NA values for windows that have no data
possible_windows <- seq(from = 0 + (window_size / 2), to = max_x, by=window_size)

# create a blank DF with all possible window midpoints
bp_windows <- data.frame(matrix(ncol = 2, nrow = length(possible_windows)))
colnames(bp_windows) <- c("window_mid", "percent_base_paired")
bp_windows$window_mid <- possible_windows
bp_windows$percent_base_paired <- NA_real_

bp_windows <- left_join(bp_windows, bp_window_values, by="window_mid")

# use NA or the real value if available
bp_windows <- bp_windows %>% mutate(percent_base_paired = 
                        if_else(is.na(pct_bp), 
                                percent_base_paired, 
                                pct_bp)) %>% select(-pct_bp)


# merge bp & cov
cov_bp_windows <- left_join(coverage_windows, bp_windows, by="window_mid")



oldest <- filter(coverage, sample_id =="1004284")

ggplot(coverage) +
  geom_line(aes(x=position, y=fwd_coverage, color=sample_id)) +
  scale_y_log10() +
  facet_wrap(~sample_id) +
  theme_bw()  +
  theme(legend.position = "none")

# some_dataset_ids <- c("Albany1902_1", "1004279", "1004283", "Davidson2006_1", "FoCo17_Pos", "PosCtrl_Pool1")
some_dataset_ids <- c("Albany1902_1", "1004283","PosCtrl_Pool1")
some_datasets <- filter(coverage, sample_id %in% some_dataset_ids)
cov_plot <- ggplot(some_datasets) +
  geom_line(aes(x=position, y=fwd_coverage, color=sample_id)) +
  scale_y_log10() +
  # facet_wrap(~sample_id, ncol=1) +
  facet_wrap(~sample_id) +
  xlab("") + 
  ylab ("Coverage of sense mapping reads") +
  theme_bw(base_size=14) +
  theme(legend.position = "none")
  
cov_plot
  
bp_plot <- ggplot (bp_windows) +
  geom_line(aes(x=window_mid, y=percent_base_paired)) +
  theme_bw(base_size = 14) +
  ylab(paste0("Fraction of ", window_size, " nt windows\nin structure\nor base-paired")) +
  xlab("Position in pre-rRNA (nt)") 

bp_plot

cov_plot + bp_plot  + plot_layout(ncol=1, heights=c(8,2))

ggsave("coverage_bp_plot.pdf", units="in", height=7.5, width=10)

# scatter plot of cov vs. %basepaired
some_datasets <- filter(cov_bp_windows, sample_id %in% some_dataset_ids)

ggplot(some_datasets) +
  geom_point(aes(x=percent_base_paired, y=mean_fwd_cov, fill=sample_id), 
             shape=21, size=1, color="black", stroke=0.25) +
  geom_smooth(aes(x=percent_base_paired, y=mean_fwd_cov)) +
  theme_bw() +
  # scale_x_log10() +
  scale_y_log10() +
  xlab("Percent of window base paired or in structure") +
  ylab("Mean sense coverage(x)") +
  facet_wrap(~sample_id)

ggplot(coverage) +
  geom_line(aes(x=position, y=rev_coverage, color=sample_id)) +
  scale_y_log10() +
  facet_wrap(~sample_id) +
  theme_bw() 

# individual base data
individual_positions <- left_join(coverage, bp)
individual_positions <- individual_positions %>% filter(!is.na(in_bp))

some_individual_positions <- filter(individual_positions, sample_id %in% some_dataset_ids)

ggplot(some_individual_positions) +
  geom_point(aes(x=in_bp, y=fwd_coverage)) +
  geom_boxplot(aes(x=in_bp, y=fwd_coverage, group=in_bp)) +
  theme_bw() + 
  scale_y_log10() +
  facet_wrap(~sample_id)
