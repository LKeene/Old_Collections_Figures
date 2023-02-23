process QUANTIFY_MISMATCHES_FROM_MAPPING {
  label 'lowmem_threaded'

  // singularity info for this process
  if (workflow.containerEngine == 'singularity'){
      container "https://depot.galaxyproject.org/singularity/perl:5.26.2"
  }     

  input:
  path (fasta)
  tuple val(meta), path (sam)

  output:
  tuple val(meta), path ("*.txt")  , emit: txt
  // path "versions.yml"  , emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:

  def args             = task.ext.args ?: ''

  // ${params.bin_dir}/tabulate_mismatches_from_sam $args -q ${params.min_mismatch_map_q_score} -f $fasta $sam > ${meta.id}.mismatches.txt

  """
  ${params.bin_dir}/tabulate_mismatches_from_sam $args -f $fasta $sam > ${meta.id}.mismatches.txt
  """

}
