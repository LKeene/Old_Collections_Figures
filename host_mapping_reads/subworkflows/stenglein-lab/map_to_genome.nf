include { BWA_MEM                                } from '../../modules/nf-core/bwa/mem/main'
include { FILTER_BAM                             } from '../../modules/stenglein-lab/filter_unmapped_sam/main'

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

  // filter out unmapped reads and remove empty bam (no mapped reads)
  FILTER_BAM(BWA_MEM.out.bam)
  ch_versions = ch_versions.mix ( FILTER_BAM.out.versions )      

 emit: 
  versions      = ch_versions
  bam           = FILTER_BAM.out.bam

}
