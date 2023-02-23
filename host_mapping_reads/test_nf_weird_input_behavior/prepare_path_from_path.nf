process prepare_path_from_path {                                                     
                                                                                
   input:                                                                       
   path(a_path)                                                                

   output:                                                                       
   tuple val ("abc"), path("*.result.txt"), emit: result_path
                                                                                
   script:                                                                      
   def bn = a_path.getBaseName()
   """                                                                          
   echo $bn > ${bn}.result.txt
   """                                                                          
}           
