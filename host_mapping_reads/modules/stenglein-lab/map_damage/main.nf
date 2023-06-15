process MAP_DAMAGE {
  tag "$meta.id"
  label 'lowmem_threaded'

 conda "bioconda::mapdamage2=2.2.1"
 container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
   'https://depot.galaxyproject.org/singularity/mapdamage2:2.2.1--pyr40_0' :
   'quay.io/biocontainers/mapdamage2:2.2.1--pyr40_0' }"

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
  samtools view $bam -h -F 20 --threads $task.cpus -o $sam

  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
  END_VERSIONS
  """

}
