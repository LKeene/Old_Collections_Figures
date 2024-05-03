process COLLECT_METADATA {
  label 'lowmem_non_threaded'

  // singularity info for this process
  if (workflow.containerEngine == 'singularity'){
      container "docker://rocker/tidyverse:4.3.2"
  }     

  input:
  path (shared_script_dir)
  path (metadata_dir)

  output:
  path "*.csv"         , emit: collected_metadata

  when:
  task.ext.when == null || task.ext.when

  script:

  def args             = task.ext.args ?: ''

  """
   Rscript ${shared_script_dir}/collect_metadata.R $metadata_dir 
  """

}
