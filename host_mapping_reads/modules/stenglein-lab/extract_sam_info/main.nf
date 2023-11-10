process EXTRACT_SAM_INFO {
  tag "$meta.id"
  label 'lowmem_non_threaded'

  input:
  tuple val(meta), path (sam)

  output:
  tuple val (meta), path ("*.sam_stats.txt")  , emit: txt
  // path "versions.yml"               , emit: versions

  when:
  task.ext.when == null || task.ext.when

  shell:
  '''
  extract_pct_identity_from_sam < !{sam} | awk '{ print "!{meta.id}" "\t" $0; }' >  !{sam}.sam_stats.txt
  '''
}
