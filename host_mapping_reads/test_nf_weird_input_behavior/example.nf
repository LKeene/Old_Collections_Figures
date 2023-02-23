#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// queue channel
numbers = Channel.of(1,2,3,4)

params.input_path="/home/datasets/a_file.txt"

process process_with_val_input {                                                     
                                                                               
  input:                                                                       
  val(a_file_name)                                                                

  output:                                                                       
  path("result.txt"), emit: result_path
                                                                               
  script:                                                                      
  """                                                                          
  echo "abc" > result.txt
  """                                                                          
}           

process process_with_path_input {                                                     

  input:                                                                       
  path(a_path)                                                                

  output:                                                                       
  path("result.txt"), emit: result_path
                                                                               
  script:                                                                      
  """                                                                          
  echo "abc" > result.txt
  """                                                                          
}           

workflow CREATE_INDEX_FROM_FILENAME {
  take:
   some_path

  main:
   ch_path = process_with_val_input(some_path).result_path

  emit:
   index = ch_path
}

workflow CREATE_INDEX_FROM_PATH {
  take:
   some_path

  main:
   ch_path = process_with_path_input(some_path).result_path

  emit:
   index = ch_path
}

process do_something_with_index {                                                     
                                                                               
  input:                                                                       
   val(number)                                                                
   path(index)                                                                
                                                                               
  script:                                                                      
   """                                                                          
   echo "def" > result.txt
   """                                                                          
}           

workflow DO_SOMETHING_WITH_INDEX {
  take:
   numbers
   index

  main:
   do_something_with_index(numbers, index)
}

workflow DO_SOMETHING_WITH_INDEX_2 {
  take:
   numbers
   index

  main:
   do_something_with_index(numbers, index)
}

workflow MAIN_WORKFLOW {

  CREATE_INDEX_FROM_FILENAME(params.input_path)
  DO_SOMETHING_WITH_INDEX(numbers, CREATE_INDEX_FROM_FILENAME.out.index)

  // path_ch = Channel.fromPath(params.input_path)
  path_ch = file(params.input_path)
  CREATE_INDEX_FROM_PATH(path_ch)
  DO_SOMETHING_WITH_INDEX_2(numbers, CREATE_INDEX_FROM_PATH.out.index)
}

workflow {
  MAIN_WORKFLOW ()
}

