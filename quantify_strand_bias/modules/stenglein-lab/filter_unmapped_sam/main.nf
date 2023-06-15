workflow FILTER_BAM {
  take:
    input
  main:
    ch_versions = Channel.empty()
    
    // input.subscribe{ println "FILTER_BAM input: $it\n" }

    FILTER_UNMAPPED_READS(input)
    ch_versions = FILTER_UNMAPPED_READS.out.versions

    FILTER_EMPTY_BAM(FILTER_UNMAPPED_READS.out.bam)

  emit:
    versions = ch_versions
    bam = FILTER_EMPTY_BAM.out.bam

}

process FILTER_UNMAPPED_READS {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::samtools=1.16.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.16.1--h6899075_1' :
        'quay.io/biocontainers/samtools:1.16.1--h6899075_1' }"

    input:
    tuple val(meta), path (input)

    output:
    tuple val(meta), path("*.filt_unmapped.bam"),  emit: bam
    path  "versions.yml",            emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    samtools \\
        view \\
        --exclude-flags 0x4 \\
        --threads ${task.cpus-1} \\
        -o ${meta.id}.filt_unmapped.bam \\
        $input 

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}

process FILTER_EMPTY_BAM {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::samtools=1.16.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.16.1--h6899075_1' :
        'quay.io/biocontainers/samtools:1.16.1--h6899075_1' }"

    input:
    tuple val(meta), path (input)

    output:
    tuple val(meta), path("*.filt.bam"),  emit: bam

    when:
    task.ext.when == null || task.ext.when

    shell:
    '''
    mapping_read_count=$(samtools stats !{input} | grep "^SN" | grep "reads mapped:" | cut -f 3)

    # only create a link to bam if there are >0 mapping reads in bam file
    if [[ $mapping_read_count -gt 0 ]]
    then
      ln !{input} !{meta.id}.filt.bam
    else
      touch !{meta.id}.filt.bam
    fi

    '''

}
