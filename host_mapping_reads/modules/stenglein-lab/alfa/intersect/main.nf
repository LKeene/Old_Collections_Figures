process ALFA_INTERSECT {
  tag "$meta.id"
  label 'lowmem_threaded'

  conda (params.enable_conda ? 'bioconda::alfa=1.1.1' : null)
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
    container "https://depot.galaxyproject.org/singularity/alfa:1.1.1--pyh5e36f6f_0"
  } else {
    container "quay.io/biocontainers/alfa:1.1.1--pyh5e36f6f_0"
  }

  input:
  tuple val(meta), path (bam), val(R1_orientation)
  val  (index_name)
  path (index_files)

  output:
  tuple val (meta), path ("*.pdf")       , emit: pdf
  tuple val (meta), path ("*.tsv")       , emit: tsv
  tuple val (meta), path ("*.bedgraph")  , emit: bedgraph
  path "versions.yml"                    , emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
  def args             = task.ext.args ?: ''

  // a flag indicating strand-specific orientation of library
  // def rev_param = R1_orientation == "reverse" ? "reverse" : "forward"

  """
  alfa -g $index_name --bam $bam $meta.id -s $R1_orientation --pdf ${meta.id}.pdf -p $task.cpus

  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
      alfa: 1.1.1 
  END_VERSIONS
  """

}
