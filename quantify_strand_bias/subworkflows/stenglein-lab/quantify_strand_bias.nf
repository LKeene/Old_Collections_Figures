include { COLLECT_METADATA                           } from '../../modules/stenglein-lab/collect_metadata'
include { BUILD_BWA_INDEX                            } from '../../subworkflows/stenglein-lab/build_bwa_index'
include { MAP_TO_GENOME                              } from '../../subworkflows/stenglein-lab/map_to_genome'
include { BAM_TO_SAM                                 } from '../../modules/stenglein-lab/bam_to_sam'
include { PREPEND_TSV_WITH_ID                        } from '../../modules/stenglein-lab/prepend_tsv_with_id'
                                                                                
include { EXTRACT_STRAND_BIAS                        } from '../../modules/stenglein-lab/extract_strand_bias'
include { PROCESS_STRAND_BIAS_OUTPUT                 } from '../../modules/stenglein-lab/process_strand_bias'

include { SAVE_OUTPUT_FILE as SAVE_OUTPUT_FILE_METADATA  } from '../../modules/stenglein-lab/save_output_file'
include { SAVE_OUTPUT_FILE                           } from '../../modules/stenglein-lab/save_output_file'
                                                                                
workflow QUANTIFY_STRAND_BIAS {                                            
                                                                                
  take:
  input     // [meta, fastq, bwa_index, R1_orientation]
  R_lib_dir // path to a directory with local R packages installed
 
  main:                                                                         

  // input.subscribe{ println "input: $it\n" }

  // make input channels in typical nf-core format of [meta, files]
  reads = input.map{ meta, fastq, bwa, ori -> [meta, fastq] }
  fasta = input.map{ meta, fastq, bwa, ori -> [meta, bwa]   }
  ori = input.map{ meta, fastq, bwa, ori -> [meta, ori]   }

  MAP_TO_GENOME(reads, fasta)

  BAM_TO_SAM(MAP_TO_GENOME.out.bam.filter{it[1].size() > 0})                                             

  BAM_TO_SAM.out.sam.join(ori).set{ori_ch}

  EXTRACT_STRAND_BIAS(ori_ch)
  PREPEND_TSV_WITH_ID(EXTRACT_STRAND_BIAS.out.txt)                              
                                                                                
  // a directory with additional R scripts
  R_script_dir_ch = Channel.fromPath(params.R_shared_script_dir)

  // directory with metadata files
  metadata_dir_ch = Channel.fromPath(params.metadata_dir)
  COLLECT_METADATA(R_script_dir_ch, metadata_dir_ch)

  // save to main results dir
  SAVE_OUTPUT_FILE_METADATA(COLLECT_METADATA.out.collected_metadata)

  // file with virus refseq metadata
  virus_refseq_ch = Channel.fromPath(params.refseq_metadata)

  PROCESS_STRAND_BIAS_OUTPUT(PREPEND_TSV_WITH_ID.out.tsv.collectFile(name: "collected_strand_bias.tsv"){it[1]}, COLLECT_METADATA.out.collected_metadata, virus_refseq_ch, R_lib_dir, R_script_dir_ch)

  // this will save main output file
  SAVE_OUTPUT_FILE(PROCESS_STRAND_BIAS_OUTPUT.out.txt)

}         
