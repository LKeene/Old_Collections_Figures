include { ALFA_INDEX                  } from '../../modules/stenglein-lab/alfa/index/main'
include { ALFA_INTERSECT              } from '../../modules/stenglein-lab/alfa/intersect/main'

workflow ALFA_ANALYSIS {

 take:
  bam                    // [meta, bam, R1_orientation]
  genome_annotation_gtf  // 

 main:

  // define some empty channels for keeping track of stuff
  ch_versions     = Channel.empty()                                               

  // prepare ALFA index
  ALFA_INDEX(genome_annotation_gtf)
  ch_versions = ch_versions.mix ( ALFA_INDEX.out.versions )      

  // prepare ALFA index
  // the filter on bam avoids running empty bam files
  ALFA_INTERSECT(bam.filter{ it[1].size() > 0 } , ALFA_INDEX.out.index_name, ALFA_INDEX.out.index_files)
  ch_versions = ch_versions.mix ( ALFA_INTERSECT.out.versions )      

 emit: 
  versions      = ch_versions
  pdf           = ALFA_INTERSECT.out.pdf       
  tsv           = ALFA_INTERSECT.out.tsv       
  bedgraph      = ALFA_INTERSECT.out.bedgraph  

}
