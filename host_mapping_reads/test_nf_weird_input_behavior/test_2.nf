#!/usr/bin/env nextflow

nextflow.enable.dsl=1


index_path = Channel.empty()

index_path = Channel.fromPath("./test.nf").collect()
// Channel.fromPath("./test.nf").set{index_path}
// index_path = Channel.fromPath("./test.nf")
reads_path = Channel.fromFilePairs("../../2023_1_OC_datasets/concatenated_datasets_R1_only/*.fastq.gz", size: -1)

// reads2.subscribe {println "reads2 value:  $it" }
// reads_path.subscribe {println "reads_path:  $it" }
// index_path.subscribe {println "index_path:  $it" }

/*
process output_components {                                                     
                                                                                
   input:                                                                       
   val(reads)                                                                
   val(index)                                                                
                                                                                
   script:                                                                      
   """                                                                          
   echo "$reads - $index" > ${reads}_${index}_abc.txt                           
   """                                                                          
}           
*/

process output_components {                                                     
                                                                                
   input:                                                                       
   tuple val (sample_id), path (fastq) from reads_path
   path(index) from index_path
                                                                                
   script:                                                                      
   def rbn = fastq.getBaseName() 
   def bn = index.getBaseName() 
   """                                                                          
   echo "$rbn" > ${rbn}_${bn}_abc.txt                           
   """                                                                          
}           
