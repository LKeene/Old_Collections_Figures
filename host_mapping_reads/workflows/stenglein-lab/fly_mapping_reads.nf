include { MARSHALL_FASTQ                             } from '../../subworkflows/stenglein-lab/marshall_fastq' 
include { SAVE_OUTPUT_FILE as SAVE_COUNTS_FILE       } from '../../modules/stenglein-lab/save_output_file'
include { PARSE_MAPPING_SAMPLESHEET                  } from '../../subworkflows/stenglein-lab/parse_mapping_samplesheet'

include { SETUP_R_DEPENDENCIES                       } from '../../modules/stenglein-lab/setup_R_dependencies' 

//genome mapping
include { BUILD_BWA_INDEX                            } from '../../subworkflows/stenglein-lab/build_bwa_index'
include { MAP_TO_GENOME                              } from '../../subworkflows/stenglein-lab/map_to_genome'

// metadata
include { COLLECT_METADATA                           } from '../../modules/stenglein-lab/collect_metadata'
include { SAVE_OUTPUT_FILE as SAVE_OUTPUT_FILE_METADATA  } from '../../modules/stenglein-lab/save_output_file'

// ALFA analysis
include { ALFA_ANALYSIS                              } from '../../subworkflows/stenglein-lab/alfa_analysis'
include { PREPEND_TSV_WITH_ID as PREPEND_ALFA_OUTPUT } from '../../modules/stenglein-lab/prepend_tsv_with_id'
include { PROCESS_ALFA_OUTPUT                        } from '../../modules/stenglein-lab/process_alfa_output'

// transcriptome mapping
include { BUILD_BWA_INDEX as BUILD_BWA_INDEX_RNA     } from '../../subworkflows/stenglein-lab/build_bwa_index'
include { MAP_TO_GENOME as MAP_TO_TRANSCRIPTOME      } from '../../subworkflows/stenglein-lab/map_to_genome'
include { BAM_TO_SAM as BAM_TO_SAM_TRANSCRIPTOME     } from '../../modules/stenglein-lab/bam_to_sam'
include { TALLY_SAM_SUBJECTS                         } from '../../modules/stenglein-lab/tally_sam_subjects'
include { PREPEND_TSV_WITH_ID as PREPEND_TALLY_OUTPUT} from '../../modules/stenglein-lab/prepend_tsv_with_id'
include { PROCESS_TRANSCRIPTOME_OUTPUT               } from '../../modules/stenglein-lab/process_transcriptome_output'

// rRNA locus mapping
include { BUILD_BWA_INDEX as BUILD_BWA_INDEX_RRNA    } from '../../subworkflows/stenglein-lab/build_bwa_index'
include { MAP_TO_GENOME as MAP_TO_RRNA_LOCUS         } from '../../subworkflows/stenglein-lab/map_to_genome'
include { BAM_TO_COV                                 } from '../../modules/stenglein-lab/bam_to_cov'
include { PREPEND_TSV_WITH_ID as PREPEND_BTC_OUTPUT  } from '../../modules/stenglein-lab/prepend_tsv_with_id'
include { PROCESS_BAM_TO_COV_OUTPUT                  } from '../../modules/stenglein-lab/process_bam_to_cov_output'

// RNA damage?
include { BAM_TO_SAM                                 } from '../../modules/stenglein-lab/bam_to_sam'
include { ERROR_CORRECTION                           } from '../../modules/stenglein-lab/error_correction'
include { MAP_TO_GENOME as MAP_FOR_DAMAGE            } from '../../subworkflows/stenglein-lab/map_to_genome'
include { QUANTIFY_MISMATCHES_FROM_MAPPING           } from '../../modules/stenglein-lab/quantify_mismatches_from_mapping'
include { PREPEND_TSV_WITH_ID as PREPEND_MISMATCH_OUTPUT } from '../../modules/stenglein-lab/prepend_tsv_with_id'
include { PROCESS_MISMATCH_OUTPUT                    } from '../../modules/stenglein-lab/process_mismatch_output'

workflow FLY_MAPPING_READS {

  MARSHALL_FASTQ(params.fastq_dir, params.fastq_pattern, params.collapse_duplicate_reads)

  SAVE_FASTQ_OUTPUT(ch_processed_reads.map{meta, reads -> reads}.flatten())

  SAVE_COUNTS_FILE(MARSHALL_FASTQ.out.fastq_counts.collectFile(name: "all_fastq_counts.txt"))

  PARSE_MAPPING_SAMPLESHEET(params.mapping_samplesheet)

  SETUP_R_DEPENDENCIES(params.R_packages)

  /* 
    ------------------
    Map to host genome
    ------------------
   */
  BUILD_BWA_INDEX(params.genome_fasta)

  MAP_TO_GENOME(MARSHALL_FASTQ.out.reads, BUILD_BWA_INDEX.out.index)

  // ch_genome_gtf = Channel.fromPath(params.genome_annotation_gtf, checkIfExists: true)
  ch_genome_gtf = file(params.genome_annotation_gtf)

  // a directory with additional R scripts
  R_script_dir_ch = Channel.fromPath(params.R_shared_script_dir)

  // directory with metadata files
  metadata_dir_ch = Channel.fromPath(params.metadata_dir)
  COLLECT_METADATA(R_script_dir_ch, metadata_dir_ch)

  // save to main results dir
  SAVE_OUTPUT_FILE_METADATA(COLLECT_METADATA.out.collected_metadata)

  // keep track of R1 orientation for strandedness for different datasets
  // PARSE_MAPPING_SAMPLESHEET.out.sample_sheet.map{ meta, R1_orientation -> [meta, R1_orientation] }.set{ori_ch}

  // join in orientation info into mapping results
  ch_mapping = MAP_TO_GENOME.out.bam.join(PARSE_MAPPING_SAMPLESHEET.out.sample_sheet)

  /*
    -------------
    ALFA analysis
    -------------
   */
  ALFA_ANALYSIS (ch_mapping, ch_genome_gtf)

  // a directory with additional R scripts
  R_script_dir_ch = Channel.fromPath(params.R_shared_script_dir)

  PREPEND_ALFA_OUTPUT(ALFA_ANALYSIS.out.tsv)
  PROCESS_ALFA_OUTPUT(PREPEND_ALFA_OUTPUT.out.tsv.collectFile(name: "collected_alfa_counts.tsv"){it[1]}, COLLECT_METADATA.out.collected_metadata, SETUP_R_DEPENDENCIES.out.R_lib_dir, R_script_dir_ch)

  /*
    ----------------------
    Transcriptome mapping
    ----------------------
   */

  // map to transcriptome
  BUILD_BWA_INDEX_RNA(params.transcriptome_fasta)
  MAP_TO_TRANSCRIPTOME(MARSHALL_FASTQ.out.reads, BUILD_BWA_INDEX_RNA.out.index)
  BAM_TO_SAM_TRANSCRIPTOME(MAP_TO_TRANSCRIPTOME.out.bam.filter{it[1].size() > 0})
  TALLY_SAM_SUBJECTS(BAM_TO_SAM_TRANSCRIPTOME.out.sam)
  PREPEND_TALLY_OUTPUT(TALLY_SAM_SUBJECTS.out.txt)
  PROCESS_TRANSCRIPTOME_OUTPUT(PREPEND_TALLY_OUTPUT.out.tsv.collectFile(name: "transcriptome_tallies.txt"){it[1]}, COLLECT_METADATA.out.collected_metadata)

  /*
    ----------------------
    Ribosomal RNA mapping
    ----------------------
   */

  // map to rRNA
  BUILD_BWA_INDEX_RRNA(params.rRNA_fasta)
  MAP_TO_RRNA_LOCUS(MARSHALL_FASTQ.out.reads, BUILD_BWA_INDEX_RRNA.out.index)
  BAM_TO_COV(MAP_TO_RRNA_LOCUS.out.bam.filter{it[1].size() > 0})
  PREPEND_BTC_OUTPUT(BAM_TO_COV.out.per_base_coverage)
  PROCESS_BAM_TO_COV_OUTPUT(PREPEND_BTC_OUTPUT.out.tsv.collectFile(name: "collected_rRNA_locus_coverage.tsv"){it[1]}, COLLECT_METADATA.out.collected_metadata, , SETUP_R_DEPENDENCIES.out.R_lib_dir, R_script_dir_ch)

  /*
    -------------------
    RNA damage analysis  
    -------------------
   */
  // Is RNA damaged?

  MAP_FOR_DAMAGE(MARSHALL_FASTQ.out.reads, BUILD_BWA_INDEX_RRNA.out.index)
  BAM_TO_SAM(MAP_FOR_DAMAGE.out.bam.filter{it[1].size() > 0})
  QUANTIFY_MISMATCHES_FROM_MAPPING(params.rRNA_fasta, BAM_TO_SAM.out.sam)
  PREPEND_MISMATCH_OUTPUT(QUANTIFY_MISMATCHES_FROM_MAPPING.out.txt)
  PROCESS_MISMATCH_OUTPUT(PREPEND_MISMATCH_OUTPUT.out.tsv.collectFile(name: "all_mismatch_counts.txt"){it[1]}, COLLECT_METADATA.out.collected_metadata, SETUP_R_DEPENDENCIES.out.R_lib_dir, R_script_dir_ch)

}
