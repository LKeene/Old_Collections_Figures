include { BUILD_BWA_INDEX                            } from '../../subworkflows/stenglein-lab/build_bwa_index'
include { MAP_TO_GENOME                              } from '../../subworkflows/stenglein-lab/map_to_genome'
include { BAM_TO_SAM                                 } from '../../modules/stenglein-lab/bam_to_sam'
include { PREPEND_TSV_WITH_ID                        } from '../../modules/stenglein-lab/prepend_tsv_with_id'
                                                                                
include { EXTRACT_STRAND_BIAS                        } from '../../modules/stenglein-lab/extract_strand_bias'
include { PROCESS_STRAND_BIAS_OUTPUT                 } from '../../modules/stenglein-lab/process_strand_bias'
                                                                                
workflow QUANTIFY_STRAND_BIAS {                                            
                                                                                
  take:
  input // [meta, fastq, bwa_index, R1_orientation]
 
  main:                                                                         

  // input.subscribe{ println "input: $it\n" }

  // make input channels in typical nf-core format of [meta, files]
  reads = input.map{ meta, fastq, bwa, ori -> [meta, fastq] }
  fasta = input.map{ meta, fastq, bwa, ori -> [meta, bwa]   }
  ori = input.map{ meta, fastq, bwa, ori -> [meta, ori]   }

  MAP_TO_GENOME(reads, fasta)

  BAM_TO_SAM(MAP_TO_GENOME.out.bam)                                             

  BAM_TO_SAM.out.sam.join(ori).set{ori_ch}

  EXTRACT_STRAND_BIAS(ori_ch)
  PREPEND_TSV_WITH_ID(EXTRACT_STRAND_BIAS.out.txt)                              
                                                                                
  PROCESS_STRAND_BIAS_OUTPUT(PREPEND_TSV_WITH_ID.out.tsv.collectFile(name: "collected_strand_bias.tsv"){it[1]}, params.metadata)
}         
