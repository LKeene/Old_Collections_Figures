#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// queue channel
numbers = Channel.of(1,2,3,4)
letters = Channel.of("A", "B", "C")

process make_files {
  input: 
   val(v)

  output:
   path("*.txt")

  script:
  """
   echo $v > ${v}.txt
  """
}

workflow FLOW1 {
  make_files(numbers)
}

workflow FLOW2 {
  make_files(letters)
}

workflow {
  FLOW1 ()
  FLOW2 ()
}

