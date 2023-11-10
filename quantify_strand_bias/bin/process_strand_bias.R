#!/usr/bin/env Rscript

library(tidyverse)
# library(rstatix)
# library(ggpubr)

# sessionInfo()

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
  output_dir           = "./"
} else {
  # if running via RStudio
  tsv_input             = "../results/collected_strand_bias.tsv"
  metadata_file         = "../refseq/metadata.csv"
  refseq_metadata_file  = "../virus_refseq/refseq_metadata.csv"
  output_dir            = "../results/process/"
}

# read in the file with the mismatch data for all datasets
datasets <- read.delim(tsv_input, header=F, sep="\t")
colnames(datasets) <- c("sample_id", "sam_file", "refseq", "fwd_count", "rev_count", "fraction_positive_strand")

datasets %>% group_by(sample_id) %>% summarize()

# only keep datasets with enough mapping reads
datasets <- datasets %>% mutate(total_count = fwd_count + rev_count)
datasets <- datasets %>% filter(total_count > 10)



# read in metadata
metadata <- read.delim(metadata_file, sep=",", header=T)
refseq_metadata <- read.delim(refseq_metadata_file, sep=",", header=T)

# join in metadata
df <- left_join(datasets, metadata, by="sample_id")
df <- left_join(df, refseq_metadata, by="refseq")

# remove seqs for which virus not defined
df %>% filter(is.na(virus))

# relevel categories
df$sample_type <- fct_relevel(df$sample_type, "Old_Collection", "Experimental_dried", "Fresh_frozen")

# rename categories
df$sample_type <- 
  recode(df$sample_type, 
         Old_Collection     = "Old\ncollections", 
         Experimental_dried = "Experimental\ndried",
         Fresh_frozen       = "Fresh\nfrozen")

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

all_p <- ggplot(df) +
  geom_jitter(aes(x=sample_type, y=fraction_positive_strand, fill=sample_type), 
              shape=21, color="black", stroke=0.2, size=1.5, height=0, width=0.25, alpha=0.5) +
  geom_boxplot(aes(x=sample_type, y=fraction_positive_strand, color=sample_type), 
               size=0.25, outlier.shape = NA, fill=NA) +
  theme_bw() +
  ylim(c(0,1)) +
  xlab("") +
  ylab("Fraction of reads from +sense RNA") +
  scale_color_manual(values = fancy_color_scale) +
  scale_fill_manual (values = fancy_color_scale) +
  facet_wrap(~virus) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme(legend.position = "none")

all_p
ggsave(paste0(output_dir, "/fraction_plus_strand_all.pdf"), width=10, height=7, units="in")

# galbut only
galbut_df <- filter(df,virus == "galbut")
galbut_p <- ggplot(galbut_df) + 
  geom_jitter(aes(x=sample_type, y=fraction_positive_strand, fill=sample_type), 
              shape=21, color="black", stroke=0.2, size=1.5, height=0, width=0.25, alpha=0.5) +
  geom_boxplot(aes(x=sample_type, y=fraction_positive_strand, color=sample_type), 
               size=0.25, outlier.shape = NA, fill=NA) +
  theme_bw(base_size=14) +
  ylim(c(0,1)) +
  xlab("") +
  scale_color_manual(values = fancy_color_scale) +
  scale_fill_manual (values = fancy_color_scale) +
  ylab("Fraction of reads from +sense RNA") +
  # theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme(legend.position = "none")

galbut_p
ggsave(paste0(output_dir, "/fraction_plus_strand_galbut.pdf"), width=6.5, height=5, units="in")

# sig different?
# galbut_t_test <- galbut_df %>% t_test(fraction_positive_strand ~ sample_type, paired=F)
# galbut_t_test

# galbut_p + stat_pvalue_manual(galbut_t_test, label = "p", y.position = 0.1)

vera_df <- df %>% 
  filter(virus == "vera" & sample_type != "Experimental\ndried") 
# get rid of unused factor levels
vera_df$sample_type <- fct_drop(vera_df$sample_type)
# vera_df %>% t_test(fraction_positive_strand ~ sample_type)

# df for stats
# only keep groups with > 3 observations
df_stats <- df %>% group_by(virus, sample_type) %>% mutate(n=n()) %>% filter(n>3)
# only keep viruses with >= 2 groups
viruses_with_enough <- df_stats %>% group_by(virus, sample_type) %>% 
  summarize() %>% group_by(virus) %>% 
  summarize (n=n(), .groups="drop") %>% filter(n>1) %>% pull(virus)
df_stats <- df_stats %>% filter(virus %in% viruses_with_enough)

# df_t_test <- df_stats %>% group_by(virus) %>% t_test(fraction_positive_strand ~ sample_type)
# write.table(df_t_test, file=paste0(output_dir, "fraction_plus_strand_t_test.txt", 
                                       # quote=F, sep="\t", 
                                       # col.names=F, row.names=F))


# code for plotting dsRNA as a fxn of week
dry_datasets <- df %>% filter(str_detect(sample_id, "Dry"))
dry_metadata <- str_match(dry_datasets$sample_id, "(\\d+)WkDry_([MF])_(\\d+)")
dry_datasets$week <- as.numeric(dry_metadata[,2])
dry_datasets$replicate <- dry_metadata[,4]

# only keep datasets with enough mapping reads
dry_datasets <- dry_datasets %>% mutate(total_count = fwd_count + rev_count)
dry_datasets <- dry_datasets %>% filter(total_count > 10)

# # average replicates
dry_datasets_avg <- dry_datasets %>%
  group_by(week, refseq) %>%
  summarize(mean_frac_pos = mean(fraction_positive_strand),
            sd_frac_pos = sd(fraction_positive_strand),
            .groups="drop")

# plot fraction + strand as a 
ggplot(dry_datasets_avg) +
  geom_point(aes(x=week, y=mean_frac_pos, fill=refseq), shape=21) +
  geom_smooth(aes(x=week, y=mean_frac_pos)) +
  geom_errorbar(aes(x=week, ymin=mean_frac_pos-sd_frac_pos, ymax=mean_frac_pos+sd_frac_pos, color=refseq), width=1.5, size=0.5) +
  theme_classic() +
  # theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust=0.5)) +
  facet_wrap(~refseq) +
  # scale_y_log10() +
  # ylim(c(0,1)) +
  xlab("Weeks samples dried") +
  ylab("Fraction of reads +strand mapping")

ggsave("strand_bias_in_experimental_dry_samples.pdf", width=10, height=7.5, units="in")
