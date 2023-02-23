process prepare_path {                                                     
                                                                                
   input:                                                                       
   val(a_path)                                                                

   output:                                                                       
   tuple val ("abc"), path("*.result.txt"), emit: result_path
                                                                                
   script:                                                                      
   // def bn = a_path.getBaseName()
   """                                                                          
   echo $a_path > ${a_path}.result.txt
   """                                                                          
}           
