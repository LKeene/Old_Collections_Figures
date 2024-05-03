include { PARSE_MAPPING_SAMPLESHEET       } from '../../subworkflows/stenglein-lab/parse_mapping_samplesheet'
include { MARSHALL_FASTQ                  } from '../../subworkflows/stenglein-lab/marshall_fastq'
include { SETUP_R_DEPENDENCIES            } from '../../modules/stenglein-lab/setup_R_dependencies'
include { BUILD_BWA_INDEX                 } from '../../subworkflows/stenglein-lab/build_bwa_index'
include { QUANTIFY_STRAND_BIAS            } from '../../subworkflows/stenglein-lab/quantify_strand_bias'

workflow STRANDEDNESS_OF_DRIED_RNA {                                                    

  // R_tar_dir_ch = Channel.fromPath(params.R_tar_dir, checkIfExists: true)
  SETUP_R_DEPENDENCIES(params.R_packages)

  MARSHALL_FASTQ(params.fastq_dir, params.fastq_pattern, params.collapse_duplicate_reads)

  PARSE_MAPPING_SAMPLESHEET(params.mapping_samplesheet)

  // pull out necessary info for building BWA Index 
  PARSE_MAPPING_SAMPLESHEET.out.sample_sheet.map{ meta, fasta, R1_orientation -> [meta, fasta] }.set{fasta_ch}

  // keep track of R1 orientation for strandedness for different datasets
  PARSE_MAPPING_SAMPLESHEET.out.sample_sheet.map{ meta, fasta, R1_orientation -> [meta, R1_orientation] }.set{ori_ch}

  BUILD_BWA_INDEX(fasta_ch)
  
  // combine matching fastq and fasta from previous outputs into a single channel
  ch_mapping_1 = MARSHALL_FASTQ.out.reads
    .join(BUILD_BWA_INDEX.out.index)
  
  // join in orientation info
  ch_mapping = ch_mapping_1.join(ori_ch)

  QUANTIFY_STRAND_BIAS(ch_mapping, SETUP_R_DEPENDENCIES.out.R_lib_dir)
}

