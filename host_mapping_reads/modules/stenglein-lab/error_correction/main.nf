process ERROR_CORRECTION {
  label 'lowmem_threaded'

  // singularity info for this process
  if (workflow.containerEngine == 'singularity'){
      container "https://depot.galaxyproject.org/singularity/bfc:r181--he4a0461_9"
  }     

  input:
  tuple val(meta) , path(reads, stageAs: 'uncorrected/*')

  output:
  tuple val(meta), path("*.fastq.gz"), emit: reads
  path "versions.yml"  , emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:

  def args             = task.ext.args ?: ''

  """
    bfc -t ${task.cpus} -s 180m ${reads.name} | gzip > ${reads.fileName.name}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bfc: \$(echo \$(bfc -v 2>&1))
    END_VERSIONS
  """

}
