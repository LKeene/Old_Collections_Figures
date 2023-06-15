process BAM_TO_SAM {
  tag "$meta.id"
  label 'lowmem_threaded'

 conda "bioconda::samtools=1.16.1"
 container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
   'https://depot.galaxyproject.org/singularity/samtools:1.16.1--h6899075_1' :
   'quay.io/biocontainers/samtools:1.16.1--h6899075_1' }"

  input:
  tuple val(meta), path (bam)

  output:
  tuple val (meta), path ("*.sam")  , emit: sam
  path "versions.yml"               , emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
  def args = task.ext.args ?: ''

  def sam  = bam.getName().replace(".bam", ".sam")

  """
  # -F 4 will not output unmapped reads
  samtools view $bam -h -F 4 --threads $task.cpus -o $sam

  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
  END_VERSIONS
  """

}
