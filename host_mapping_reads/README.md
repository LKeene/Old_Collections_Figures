This directory contains code used to analyze host-mapping (fly-mapping) reads in Old Collection datasets.  

The main questions we sought to address are:

- What kinds of host RNA are present in old samples?
- Is there evidence of preferential survival of dsRNA?
- If dsRNA preferentially survives, what is the source of dsRNA?
- Is there evidence of RNA damage?

This directory contains a self-contained nextflow workflow to address these questions.  The main entry point to this workflow is the run_datasets shell script.


### Software dependencies

These analyses are implemented in [nextflow](https://www.nextflow.io/docs/latest/).  Dependencies are handled through singularity so installation of other software besides nextflow and singularity shouldn't be necessary.  

### Reference sequence dependnecies

You will need to download the *D. melanogaster* [reference genome](https://ftp.ensembl.org/pub/release-108/fasta/drosophila_melanogaster/dna/Drosophila_melanogaster.BDGP6.32.dna.toplevel.fa.gz) and [annotation](https://ftp.ensembl.org/pub/release-108/gtf/drosophila_melanogaster/Drosophila_melanogaster.BDGP6.32.108.gtf.gz) from the indicated links.  The pipeline expects these files to be in the refseq directory.
