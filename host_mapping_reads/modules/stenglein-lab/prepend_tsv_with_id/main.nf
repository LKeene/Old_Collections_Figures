process PREPEND_TSV_WITH_ID {
    tag "$tsv"
    label 'process_low'

    // we just need a base linux environment for this module
    // which is an assumption of the nextflow pipeline

    input:
    tuple val(meta), path(tsv)                                                
    // val (prefix)

    output:
    tuple val(meta), path ("*.prepended.tsv") , emit: tsv
    path ("*prepended.tsv") ,                   emit: tsv_file
    path "versions.yml" ,                       emit: versions

    when:
    task.ext.when == null || task.ext.when


    shell:

    // def prefix   = task.ext.prefix ?: "${meta.id}"                                
    // def tsv_base = tsv.getBaseName()

    // prepend tsv output with a column containing sample ID (from meta.id)
    // ignore lines beginning with # (comment lines)
    '''
    grep -v -e "^#" !{tsv} | awk '{print "!{meta.id}" "\t" $0}' > "!{tsv}.prepended.tsv"
    
    cat <<-END_VERSIONS > versions.yml
    "!{task.process}":
    END_VERSIONS
    '''
}
