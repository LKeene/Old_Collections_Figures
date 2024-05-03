process PROCESS_MISMATCH_OUTPUT {
  label 'lowmem_non_threaded'

  // singularity info for this process
  if (workflow.containerEngine == 'singularity'){
      container "docker://rocker/tidyverse:4.3.2"
  }     

  input:
  path (txt)
  path (metadata)
  path (R_lib_dir)
  path (R_shared_script_dir)


  output:
  path "*.pdf"         , emit: pdf
  path "*.txt"         , emit: txt, includeInputs: true
  // path "versions.yml"  , emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:

  def args             = task.ext.args ?: ''

  """
   # Rscript ${params.bin_dir}/analyze_mismatch_types.R $txt $metadata $R_lib_dir $R_shared_script_dir
   analyze_mismatch_types.R $txt $metadata $R_lib_dir $R_shared_script_dir
  """

}
