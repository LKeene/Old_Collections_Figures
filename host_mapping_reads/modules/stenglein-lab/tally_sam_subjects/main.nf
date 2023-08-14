process TALLY_SAM_SUBJECTS {
  tag "$meta.id"
  label 'lowmem_non_threaded'

  input:
  tuple val(meta), path (sam)

  output:
  tuple val (meta), path ("*.txt")  , emit: txt
  // path "versions.yml"               , emit: versions

  when:
  task.ext.when == null || task.ext.when

  shell:
  '''
  # awk will ignore header lines (starting w/ @) and unmapped reads (field 3 == *)
  awk '( $0 !~ /^@/ ) && ( $3 != "*" ) { print }' !{sam} |   
  cut -f 3 |                 # cut out 3rd field (subject ID)
  sort |                     # sort subject IDs
  uniq -c |                  # collapse to unique list of subject IDs and tally counts
  sort -nr |                 # sort numerically by tally, largest first
  awk '{print $2 "\t" $1}' >  !{sam}.subjects.txt # reverse column order
  '''

}
