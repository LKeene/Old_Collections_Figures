


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
} else {
  # if running via RStudio
  # r_lib_dir = "../lib/R/"
  tsv_input      = "../results/process/collected_alfa_counts.tsv"
  metadata_input = "../refseq/metadata.csv"
}

library(tidyverse)

# read in metadata
metadata <- read.delim(metadata_input, sep=",", 
                       header=T, stringsAsFactors = F)

datasets <- read.delim(file = tsv_input, sep="\t")
colnames(datasets) <- c("sample_id", "type", "counts", "size_in_genome")

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

# parse out category,biotype separated by a comma in column #2
type_split <- str_match(datasets$type, "(.*),(.*)")
datasets$category <- type_split[,2]
datasets$biotype  <- type_split[,3]

datasets <- datasets %>% select(-type)

# calculate normalized counts 
datasets <- datasets %>% 
  group_by(sample_id) %>% 
  mutate(fractional_counts = counts / sum(counts)) 

# reorder sample factors for display
datasets$sample_id <- fct_reorder(datasets$sample_id, datasets$date_collected, min)

# don't plot negative control dataset (very few reads)
datasets <- datasets %>% filter(control_type != "Negative")

# biotypes - fill by control type
ggplot(datasets) +
  geom_col(aes(x=sample_id, y=fractional_counts, fill=control_type)) +
  scale_fill_manual(values=c("grey", "cornflowerblue", "Orange")) +
  theme_bw() + 
  xlab("Samples ordered by collection date") +
  ylab("Percent reads mapping to indicated feature type") +
  # scale_y_log10() +
  facet_wrap(~biotype) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

ggsave("alfa_biotype_plot.pdf", width=10, height=7, units="in")

# biotypes - focus on opposite strand
ggplot(filter(datasets, biotype == "opposite_strand" )) +
  geom_col(aes(x=sample_id, y=fractional_counts, fill=control_type)) +
  scale_fill_manual(values=c("grey", "cornflowerblue", "Orange")) +
  theme_bw() + 
  xlab("Samples ordered by collection date") +
  ylab("Percent mapping reads") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

ggsave("alfa_biotype_plot_opposite_strand.pdf", width=10, height=7, units="in")

# biotypes - focus on opposite strand
ggplot(filter(datasets, biotype == "opposite_strand" )) +
  geom_col(aes(x=sample_id, y=fractional_counts, fill=storage_type)) +
  scale_fill_manual(values=c("cornflowerblue", "orange", "grey")) +
  theme_bw() + 
  xlab("Samples ordered by collection date") +
  ylab("Percent mapping reads") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

ggsave("alfa_biotype_plot_opposite_strand_color_by_storage.pdf", width=10, height=7, units="in")


ggsave("alfa_biotype_plot_opposite_strand.pdf", width=10, height=7, units="in")
# biotypes - fill by species
ggplot(datasets) + 
  geom_col(aes(x=sample_id, y=fractional_counts, fill=species)) +
  # scale_fill_manual(values=c("grey", "cornflowerblue")) +
  theme_bw() + 
  # scale_y_log10() +
  facet_wrap(~biotype) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

ggsave("alfa_biotype_plot_color_by_species.pdf", width=10, height=7, units="in")
