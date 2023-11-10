process TABULATE_MAPPING_STATS {
  tag "all"
  label 'lowmem_non_threaded'

  // singularity info for this process
  if (workflow.containerEngine == 'singularity'){
      container "docker://rocker/tidyverse:4.2.2"
  }

  input:
  path (sam_stats)

  output:
  path ("summarized_mapping_stats.txt")  , emit: txt

  when:
  task.ext.when == null || task.ext.when

  script:
  """
  Rscript ${params.bin_dir}/tabulate_mapping_stats.R $sam_stats 
  """

}
