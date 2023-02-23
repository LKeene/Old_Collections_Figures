process BAM_TO_COV {
  tag "$meta.id"
  label 'lowmem_threaded'

  conda (params.enable_conda ? 'bioconda::bamtocov=2.7.0' : null)
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
    container "https://depot.galaxyproject.org/singularity/bamtocov:2.7.0--hbd632db_1"
  } else {
    container "quay.io/biocontainers/bamtocov:2.7.0--hbd632db_1"
  }

  input:
  tuple val(meta), path (bam)

  output:
  tuple val (meta), path ("*.coverage_report.txt")  , emit: coverage_report
  tuple val (meta), path ("*.coverage.txt")         , emit: coverage
  path "versions.yml"                               , emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
  def args             = task.ext.args ?: ''

  """
  bamtocov --stranded --threads $task.cpus -o ${meta.id}.coverage_report.txt $bam  > ${meta.id}.coverage.txt

  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
      bamtocov: 2.7.0 
  END_VERSIONS
  """

}
