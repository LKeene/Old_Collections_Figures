process output_components {                                                     
                                                                                
   input:                                                                       
   val(reads)                                                                
   tuple val(meta), path(index)                                                                
   // tuple val(reads) val(index)                                                  
                                                                                
   script:                                                                      
   def bn = index.getBaseName()
   """                                                                          
   echo "$reads" > ${reads}_${bn}_abc.txt                           
   """                                                                          
}           
