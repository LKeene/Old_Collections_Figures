process PROCESS_STRAND_BIAS_OUTPUT {
  label 'lowmem_non_threaded'

  // singularity info for this process
  if (workflow.containerEngine == 'singularity'){
      container "docker://rocker/tidyverse:4.2.2"
  }     

  input:
  path (txt)
  path (metadata)
  path (refseq_metadata)

  output:
  path "*.pdf"         , emit: pdf
  path (txt)           , emit: txt, includeInputs: true
  // path "versions.yml"  , emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:

  def args             = task.ext.args ?: ''

  """
   Rscript ${params.bin_dir}/process_strand_bias.R $txt $metadata $refseq_metadata
  """

}
