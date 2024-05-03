## Analyis of strand bias in virus-mapping reads

This directory contains code used to characterize ratios of +strand and -strand virus-mapping reads in NGS datasets.

This is implemented as a nextflow workflow.  The main entry point to this workflow is the run_strand_bias_workflow shell script

### Software dependencies

These analyses are implemented in [nextflow](https://www.nextflow.io/docs/latest/).  Dependencies are handled through singularity so installation of other software besides nextflow and singularity shouldn't be necessary.  

