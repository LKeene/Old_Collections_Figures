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
  refseq_metadata_file        = args[3]
  R_lib_dir                   = args[4]
  R_script_dir                = args[5]
  output_dir                  = "./"
} else {
  # if running via RStudio
  tsv_input                   = "../results/collected_virus_coverage.tsv"
  metadata_input              = "../results/collect/collected_metadata.csv"
  refseq_metadata_file  = "../../metadata/virus_refseq_metadata.csv"
  R_lib_dir                   = "../lib/R/"
  R_script_dir                = "../../scripts/"
  output_dir                  = "../results/plot/"
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

# virus reference sequence metadata 
refseq_metadata <- read.delim(refseq_metadata_file, sep=",", header=T)

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

# merge in virus refseq metadata
coverage <- left_join(coverage, refseq_metadata, by="refseq")

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

# calculate total coverage = fwd + rev cov
coverage <- coverage %>% mutate(total_cov = fwd_cov + rev_cov)

# print(coverage %>% group_by(sample_id) %>% summarize(), n=100)

# calculate coverage averages
coverage_averages <- 
  coverage %>% 
  group_by(sample_id) %>% 
  summarize(median_total_cov    = median(total_cov, na.rm = T),
            mean_total_cov      = mean  (total_cov, na.rm = T),
            sd_total_cov        = sd    (total_cov, na.rm = T))

coverage_averages <- left_join(coverage_averages, metadata)

# filter out datasets for which no defined sample type metadata: neg ctrl samples
coverage_averages <- coverage_averages %>% filter (!is.na(sample_type))
            
# create coverage windows 
window_size <- 5
coverage <- coverage %>% 
  mutate(window_start = window_size * (position %/% window_size),
         window_mid   = window_start + (window_size / 2))

# calculate % base paired in each window
coverage_windows <- coverage %>% group_by(sample_id, refseq, window_mid) %>% 
  summarize(mean_fwd_cov = mean(fwd_cov), 
            mean_rev_cov = mean(rev_cov),
            mean_total_cov = mean(total_cov),
            .groups="drop")

# merge metadata into counts table
coverage_windows <- left_join(coverage_windows, metadata, by="sample_id")
coverage_windows <- left_join(coverage_windows, refseq_metadata, by="refseq")

# focus on galbut virus 

# -----------------------
# plot windowed coverage
# -----------------------
galbut_cov_windows <- filter(coverage_windows, virus=="Galbut_virus" & sample_type == "Museum\nsamples") 

galbut_cov_windows_longer <- galbut_cov_windows %>% pivot_longer(cols=c(mean_fwd_cov, mean_rev_cov), names_to="cov_strand", values_to="mean_cov")

# rename cov to something nicer
galbut_cov_windows_longer$cov_strand <- case_match(galbut_cov_windows_longer$cov_strand, "mean_fwd_cov" ~ "+strand", "mean_rev_cov" ~ "-strand", .default = NA)
galbut_cov_windows_longer$cov_strand <- fct_relevel(galbut_cov_windows_longer$cov_strand, "+strand", "-strand")

# -------------
# non-windowed 
# -------------
galbut_cov <- filter(coverage, virus=="Galbut_virus" & sample_type == "Museum\nsamples") 
galbut_cov_longer <- galbut_cov %>% pivot_longer(cols=c(fwd_cov, rev_cov), names_to="cov_strand", values_to="cov")

# rename cov to something nicer
galbut_cov_longer$cov_strand <- case_match (galbut_cov_longer$cov_strand, "fwd_cov" ~ "+strand", "rev_cov" ~ "-strand", .default = NA)
galbut_cov_longer$cov_strand <- fct_relevel(galbut_cov_longer$cov_strand, "+strand", "-strand")

left_1600  <- 364
right_1601 <- 722
left_1948  <- 1089
right_1949 <- 1165

# the location of the galbut virus-targeting primers that we routinely use to detect
# galbut virus RNA1 and their location in RNA1 (with common coordinates for all sequences)
# note to achieve common coordinates I put N bases at the beginning of several galbut virus RNA1 sequences
primer_loc <- tribble(
  ~name,        ~ymin,  ~ymax, ~xmin, ~xmax, ~segment, ~cov_strandxx,
   "1600_1601", 1e0,   1e3,   left_1600, right_1601,   "RNA_1",  "+strand",
   "1948_1949", 1e0,   1e3,   left_1948, right_1949,   "RNA_1",  "+strand"
)
primer_loc

# a color-blind friendly qualitative color scale
# from: https://colorbrewer2.org/#type=diverging&scheme=RdYlBu&n=11
color_scale <- c(
rgb(165,0,38, maxColorValue = 255),
rgb(244,109,67, maxColorValue = 255),
rgb(253,174,97, maxColorValue = 255),
rgb(171,217,233, maxColorValue = 255),
rgb(69,117,180, maxColorValue = 255),
rgb(49,54,149, maxColorValue = 255))

# plot galbut virus coverage - non-windowed
ggplot(galbut_cov_longer) +
  geom_line(aes(x=position, y=cov, color=sample_id_in_paper), linewidth = 0.5) + 
  geom_rect(data=primer_loc, 
            aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="slateblue", alpha=0.25) +
  theme_this_paper() +
  scale_color_manual(values=color_scale) +
  scale_y_log10(limits=c(1e0, 1e3)) +
  theme(legend.position="right", legend.title=element_blank()) +
  facet_grid(cov_strand~segment) +
  xlab("Position in genome segment") +
  ylab("Mean coverage in 20 base window on indicated strand")


# plot galbut virus coverage - windowed
ggplot(galbut_cov_windows_longer) +
  geom_line(aes(x=window_mid, y=mean_cov, color=sample_id_in_paper), linewidth=0.35) + 
  geom_rect(data=primer_loc, 
            aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="slateblue", alpha=0.25) +
  theme_this_paper() +
  scale_color_manual(values=color_scale) +
  scale_y_log10(limits=c(1e0, 1e3)) +
  scale_x_continuous(breaks=c(0,1500) ) +
  theme(legend.position="right", legend.title=element_blank()) +
  facet_grid(cov_strand~segment) +
  xlab("Position in genome segment") +
  ylab("Mean coverage in 20 base window on indicated strand")

ggsave(paste0(output_dir, "Fig_SX_galbut_virus_coverage.pdf"), units="in", width=7.5, height=6)
