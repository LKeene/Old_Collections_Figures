include { BWA_MEM                     } from '../../modules/nf-core/bwa/mem/main'

workflow MAP_TO_GENOME {

 take:
  reads
  index

 main:

  // define some empty channels for keeping track of stuff
  ch_versions     = Channel.empty()                                               

  // map reads 
  def sort_bam = true 
  BWA_MEM (reads, index, sort_bam)
  ch_versions = ch_versions.mix ( BWA_MEM.out.versions )      

 emit: 
  versions      = ch_versions
  bam           = BWA_MEM.out.bam

}
