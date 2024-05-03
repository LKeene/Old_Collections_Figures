process PROCESS_ALFA_OUTPUT {
  label 'lowmem_threaded'

  // singularity info for this process
  if (workflow.containerEngine == 'singularity'){
      container "docker://rocker/tidyverse:4.3.2"
  }     

  input:
  path (tsv)
  path (metadata)
  path (R_lib_dir)
  path (R_script_dir)

  output:
  path "*.pdf"         , emit: pdf
  path "*.tsv"         , emit: tsv, includeInputs: true
  // path "versions.yml"  , emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:

  def args             = task.ext.args ?: ''

  """
   Rscript ${params.bin_dir}/process_alfa_output.R $tsv $metadata $R_lib_dir $R_script_dir
  """

}
