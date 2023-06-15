process FAST_UNIQ {
  tag "$meta.id"
  label 'process_low'

  conda (params.enable_conda ? 'bioconda::fastuniq=1.1' : null)
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
    container "https://depot.galaxyproject.org/singularity/fastuniq:1.1--h470a237_1"
  } else {
    container "quay.io/biocontainers/fastuniq:1.1--h470a237_1"
  }

  input:
  tuple val(meta), path(reads)

  output:
  tuple val(meta), path('*.trim.uniq.fastq.gz') , emit: reads
  path "versions.yml"                      , emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:

  def args             = task.ext.args ?: ''
  def prefix           = task.ext.prefix ?: "${meta.id}"

  def r1 = reads[0]
  def r2 = reads[1]

  def r1_out = r1.replaceAll(".fastq.gz$", ".trim.fastq.gz"
  def r2_out = r2 ? r2.replaceAll(".fastq.gz$", ".trim.fastq.gz" : ""

  def input_files    = reads[1] ? "-i $r1 $r2" 
                                : "-i $r1"
  def output_files   = reads[1] ? "-o $r1_out -p $r2_out"
                                : "-o $r1_out"


  """
  fastuniq \
   $args
   $input_files \
   $output_files 

   cat <<-END_VERSIONS > versions.yml
   "${task.process}":
       cd-hit-dup: 4.8
   END_VERSIONS
  """

}
