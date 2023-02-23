include { PREPROCESS_READS                           } from '../../subworkflows/stenglein-lab/preprocess_reads'

//genome mapping
include { BUILD_BWA_INDEX                            } from '../../subworkflows/stenglein-lab/build_bwa_index'
include { MAP_TO_GENOME                              } from '../../subworkflows/stenglein-lab/map_to_genome'

// ALFA analysis
include { ALFA_ANALYSIS                              } from '../../subworkflows/stenglein-lab/alfa_analysis'
include { PREPEND_TSV_WITH_ID as PREPEND_ALFA_OUTPUT } from '../../modules/stenglein-lab/prepend_tsv_with_id'
include { PROCESS_ALFA_OUTPUT                        } from '../../modules/stenglein-lab/process_alfa_output'

// transcriptome mapping
include { BUILD_BWA_INDEX as BUILD_BWA_INDEX_RNA     } from '../../subworkflows/stenglein-lab/build_bwa_index'
include { MAP_TO_GENOME as MAP_TO_TRANSCRIPTOME      } from '../../subworkflows/stenglein-lab/map_to_genome'

// rRNA locus mapping
include { BUILD_BWA_INDEX as BUILD_BWA_INDEX_RRNA    } from '../../subworkflows/stenglein-lab/build_bwa_index'
include { MAP_TO_GENOME as MAP_TO_RRNA_LOCUS         } from '../../subworkflows/stenglein-lab/map_to_genome'
include { BAM_TO_COV                                 } from '../../modules/stenglein-lab/bam_to_cov'
include { PREPEND_TSV_WITH_ID as PREPEND_BTC_OUTPUT  } from '../../modules/stenglein-lab/prepend_tsv_with_id'
include { PROCESS_BAM_TO_COV_OUTPUT                  } from '../../modules/stenglein-lab/process_bam_to_cov_output'

// RNA damage?
include { BAM_TO_SAM                                 } from '../../modules/stenglein-lab/bam_to_sam'
include { QUANTIFY_MISMATCHES_FROM_MAPPING           } from '../../modules/stenglein-lab/quantify_mismatches_from_mapping'
include { PREPEND_TSV_WITH_ID as PREPEND_MISMATCH_OUTPUT } from '../../modules/stenglein-lab/prepend_tsv_with_id'
include { PROCESS_MISMATCH_OUTPUT                    } from '../../modules/stenglein-lab/process_mismatch_output'

// virus mapping for strand bias 
include { BUILD_BWA_INDEX as BUILD_BWA_INDEX_GALBUT  } from '../../subworkflows/stenglein-lab/build_bwa_index'
include { MAP_TO_GENOME as MAP_TO_GALBUT             } from '../../subworkflows/stenglein-lab/map_to_genome'
// include { BAM_TO_SAM                                 } from '../../modules/stenglein-lab/bam_to_sam'
// include { PROCESS_STRAND_BIAS_OUTPUT                 } from '../../modules/stenglein-lab/process_strand_bias'

workflow FLY_MAPPING_READS {                                                    

  PREPROCESS_READS()

  BUILD_BWA_INDEX(params.genome_fasta)

  MAP_TO_GENOME(PREPROCESS_READS.out.reads, BUILD_BWA_INDEX.out.index)

  // ch_genome_gtf = Channel.fromPath(params.genome_annotation_gtf, checkIfExists: true)
  ch_genome_gtf = file(params.genome_annotation_gtf)

  ALFA_ANALYSIS (MAP_TO_GENOME.out.bam, ch_genome_gtf)

  PREPEND_ALFA_OUTPUT(ALFA_ANALYSIS.out.tsv)
  PROCESS_ALFA_OUTPUT(PREPEND_ALFA_OUTPUT.out.tsv.collectFile(name: "collected_alfa_counts.tsv"){it[1]}, params.metadata)

  // map to transcriptome
  BUILD_BWA_INDEX_RNA(params.transcriptome_fasta)
  MAP_TO_TRANSCRIPTOME(PREPROCESS_READS.out.reads, BUILD_BWA_INDEX_RNA.out.index)

  // map to rRNA
  BUILD_BWA_INDEX_RRNA(params.rRNA_fasta)
  MAP_TO_RRNA_LOCUS(PREPROCESS_READS.out.reads, BUILD_BWA_INDEX_RRNA.out.index)
  BAM_TO_COV(MAP_TO_RRNA_LOCUS.out.bam)
  PREPEND_BTC_OUTPUT(BAM_TO_COV.out.coverage)
  // TODO: fix R script here
  //PROCESS_BAM_TO_COV_OUTPUT(PREPEND_BTC_OUTPUT.out.tsv.collectFile(name: "collected_rRNA_locus_coverage.tsv"){it[1]}, params.metadata)

  // Is RNA damaged?
  // quantify mismatches in rRNA-mapping reads
  BAM_TO_SAM(MAP_TO_RRNA_LOCUS.out.bam)
  // rRNA_ch = Channel.fromPath(params.rRNA_fasta)
  QUANTIFY_MISMATCHES_FROM_MAPPING(params.rRNA_fasta, BAM_TO_SAM.out.sam)
  PREPEND_MISMATCH_OUTPUT(QUANTIFY_MISMATCHES_FROM_MAPPING.out.txt)
  PROCESS_MISMATCH_OUTPUT(PREPEND_MISMATCH_OUTPUT.out.tsv.collectFile(name: "all_mismatch_counts.txt"){it[1]}, params.metadata)

  // Strand bias?
  // map to galbut virus
  // BUILD_BWA_INDEX_GALBUT(params.galbut_virus_fasta)
  // MAP_TO_GALBUT(PREPROCESS_READS.out.reads, BUILD_BWA_INDEX_GALBUT.out.index)
  // BAM_TO_SAM(MAP_TO_GALBUT.out.bam)
  // PROCESS_STRAND_BIAS_OUTPUT(BAM_TO_SAM.out.sam)
  // PROCESS_BAM_TO_COV_OUTPUT(PREPEND_BTC_OUTPUT.out.tsv.collectFile(name: "collected_rRNA_locus_coverage.tsv"){it[1]}, params.metadata)

}

