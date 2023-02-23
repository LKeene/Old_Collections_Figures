process ALFA_INDEX {
  // tag "$meta.id"
  label 'lowmem_threaded'

  conda (params.enable_conda ? 'bioconda::alfa=1.1.1' : null)
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
    container "https://depot.galaxyproject.org/singularity/alfa:1.1.1--pyh5e36f6f_0"
  } else {
    container "quay.io/biocontainers/alfa:1.1.1--pyh5e36f6f_0"
  }

  input:
  path(gtf)

  output:
  val  "alfa_index"       , emit: index_name
  path "alfa_index*"      , emit: index_files
  path "versions.yml"     , emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:

  def args             = task.ext.args ?: ''

  """
  # sort GTF, per ALFA manual
  sort -k1,1 -k4,4n -k5,5nr $gtf > ${gtf}.sorted.gtf

  # make alfa index(es)
  alfa -a ${gtf}.sorted.gtf -g alfa_index

  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
      alfa: 1.1.1 
  END_VERSIONS
  """

}
