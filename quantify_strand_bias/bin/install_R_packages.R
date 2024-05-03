#!/usr/bin/env Rscript 
#
# This script installs R packages needed for the pipeline
# from downloaded tarballs
#
# This is useful when pipeline R scripts require packages
# that aren't supplied as part of a base singularity R image
#
# Mark Stenglein 2/29/2024 
#

#
# This code block sets up input arguments to either come from the command line
# (if running from the pipeline, so not interactively) or to use expected default values 
# (if running in interactive mode, for example in RStudio for troubleshooting
# or development).  
#
if (!interactive()) {
  # if running from Rscript
  args = commandArgs(trailingOnly=TRUE)
  r_lib_dir=args[1]
  packages=args[-1]
} else {
  # if running via RStudio
  r_lib_dir = "../lib/R/"
  packages="rstatix"
}

# for each package
for (package in packages) {
  # install package from tarball into specified library directory
  install.packages(package, lib=r_lib_dir)

  # confirm it worked
  library(package, character.only=T, lib.loc=r_lib_dir)
}


