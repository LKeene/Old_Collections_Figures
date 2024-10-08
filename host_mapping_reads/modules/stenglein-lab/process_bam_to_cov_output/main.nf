process PROCESS_BAM_TO_COV_OUTPUT {
  label 'lowmem_threaded'

  // singularity info for this process
  if (workflow.containerEngine == 'singularity'){
      container "docker://rocker/tidyverse:4.2.2"
  }     

  input:
  path (txt)
  path (metadata)

  output:
  path "*.pdf"         , emit: pdf, optional: true
  path "*.txt"         , emit: txt, optional: true, includeInputs: true
  // path "versions.yml"  , emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:

  def args             = task.ext.args ?: ''

  """
   Rscript ${params.bin_dir}/process_bam_to_cov_output.R $txt $metadata
  """
}
