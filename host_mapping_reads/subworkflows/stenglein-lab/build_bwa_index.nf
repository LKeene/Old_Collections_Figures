include { BWA_INDEX                   } from '../../modules/nf-core/bwa/index/main'

workflow BUILD_BWA_INDEX {

 take:

  genome_fasta

 main:

  // define some empty channels for keeping track of stuff
  ch_versions     = Channel.empty()                                               

  // prepare a genome index 
  // the nf-core bwa index module has a tuple input val(meta), path(fasta), so make 
  // some a metadata string to make it happy
  ch_index = Channel.value("bwa-index-meta").map{it -> [it, file(genome_fasta)]}

  BWA_INDEX (ch_index)
  ch_versions = ch_versions.mix ( BWA_INDEX.out.versions )      
  
  
 emit: 
  versions      = ch_versions
  index         = BWA_INDEX.out.index

}
