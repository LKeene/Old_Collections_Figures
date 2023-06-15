process EXTRACT_STRAND_BIAS {
  label 'lowmem_threaded'
  tag "$meta.id"

  // singularity info for this process
  if (workflow.containerEngine == 'singularity'){
      container "https://depot.galaxyproject.org/singularity/perl:5.26.2"
  }     

  input:
  tuple val(meta), path (sam), val(R1_orientation)

  output:
  tuple val(meta), path ("*.txt")  , emit: txt
  // path "versions.yml"  , emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:

  def args             = task.ext.args ?: ''

  // a flag indicating strand-specific orientation of library
  def rev_param = R1_orientation == "reverse" ? "-r" : ""

  """
  ${params.bin_dir}/extract_strand_bias_from_sam $rev_param $args $sam > ${meta.id}.strand_bias.txt
  """

}
