process EXTRACT_READ_MAPPING_COUNTS {
    tag "$bam"
    label 'process_low'

    conda (params.enable_conda ? 'bioconda::samtools=1.15.1' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.15.1--h1170115_0' :
        'quay.io/biocontainers/samtools:1.15.1--h1170115_0' }"

    input:
    tuple val(meta), path(bam)                                                

    output:
    tuple val(meta),  path ("${bam}.txt"),   emit: tsv
    path "${bam}.txt" ,                      emit: read_count_table
    path "versions.yml" ,                    emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    samtools coverage $args $bam | grep -v ^# | cut -f 1,4 > ${bam}.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
