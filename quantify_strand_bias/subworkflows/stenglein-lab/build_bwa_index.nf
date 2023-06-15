include { BWA_INDEX                   } from '../../modules/nf-core/bwa/index/main'

workflow BUILD_BWA_INDEX {

 take:

  genome_fasta  // [meta, fasta]  

 main:

  // define some empty channels for keeping track of stuff
  ch_versions     = Channel.empty()                                               

  BWA_INDEX (genome_fasta)

  ch_versions = ch_versions.mix ( BWA_INDEX.out.versions )      
  
  
 emit: 
  versions      = ch_versions
  index         = BWA_INDEX.out.index

}
