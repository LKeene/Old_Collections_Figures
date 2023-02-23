process COMPRESS_INPUT {
    tag "$meta.id"
    label 'process_low'

    // we just need a base linux environment with gzip for this module
    // which is an assumption of the nextflow pipeline
    // so don't need to define a singularity image or conda environment
    // I think (?)

    input:
    tuple val(meta), path(files)                                                

    output:
    tuple val(meta), path ("*.gz") , emit: compressed_files
    path "versions.yml" ,             emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    if (meta.single_end) {                                                      
      """
      gzip -q $files

      cat <<-END_VERSIONS > versions.yml
      "${task.process}":
              gzip: \$( gzip --version | sed -e "s/FastQC v//g" )             
      END_VERSIONS
      """
    } else {
      """
      gzip -q $files[0]
      gzip -q  $files[1]

      cat <<-END_VERSIONS > versions.yml
      "${task.process}":
              gzip: \$( gzip --version | sed -e "s/FastQC v//g" )             
      END_VERSIONS
      """
    }
}
