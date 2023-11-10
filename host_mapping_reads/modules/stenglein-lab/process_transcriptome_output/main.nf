process PROCESS_TRANSCRIPTOME_OUTPUT {
  label 'lowmem_non_threaded'

  // singularity info for this process
  if (workflow.containerEngine == 'singularity'){
      container "docker://rocker/tidyverse:4.2.2"
  }     

  input:
  path (txt)
  path (metadata)

  output:
  // path "*.pdf"         , emit: pdf
  path "*.txt"         , emit: txt, includeInputs: true
  // path "versions.yml"  , emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:

  def args             = task.ext.args ?: ''

  """
   Rscript ${params.bin_dir}/process_transcriptome_counts.R $txt $metadata
  """

}
