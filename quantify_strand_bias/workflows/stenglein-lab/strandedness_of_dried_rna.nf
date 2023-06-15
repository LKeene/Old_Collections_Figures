include { PARSE_MAPPING_SAMPLESHEET                  } from '../../subworkflows/stenglein-lab/parse_mapping_samplesheet'
include { PREPROCESS_READS                           } from '../../subworkflows/stenglein-lab/preprocess_reads'
include { QUANTIFY_STRAND_BIAS                       } from '../../subworkflows/stenglein-lab/quantify_strand_bias'
include { BUILD_BWA_INDEX                            } from '../../subworkflows/stenglein-lab/build_bwa_index'

workflow STRANDEDNESS_OF_DRIED_RNA {                                                    

  PREPROCESS_READS(params.fastq_dir, params.fastq_pattern, params.collapse_duplicate_reads)

  PARSE_MAPPING_SAMPLESHEET(params.mapping_samplesheet)

  // pull out necessary info for building BWA Index 
  PARSE_MAPPING_SAMPLESHEET.out.sample_sheet.map{ meta, fasta, R1_orientation -> [meta, fasta] }.set{fasta_ch}

  // keep track of R1 orientation for strandedness for different datasets
  PARSE_MAPPING_SAMPLESHEET.out.sample_sheet.map{ meta, fasta, R1_orientation -> [meta, R1_orientation] }.set{ori_ch}

  BUILD_BWA_INDEX(fasta_ch)
  
  // combine matching fastq and fasta from previous outputs into a single channel
  ch_mapping_1 = PREPROCESS_READS.out.reads
    .join(BUILD_BWA_INDEX.out.index)
  
  // join in orientation info
  ch_mapping = ch_mapping_1.join(ori_ch)

  QUANTIFY_STRAND_BIAS(ch_mapping)
}

