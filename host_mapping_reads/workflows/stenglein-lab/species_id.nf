include { PREPROCESS_READS                           } from '../../subworkflows/stenglein-lab/preprocess_reads'
include { BUILD_BWA_INDEX                            } from '../../subworkflows/stenglein-lab/build_bwa_index'
include { MAP_TO_GENOME as MAP_TO_REFERENCE          } from '../../subworkflows/stenglein-lab/map_to_genome'
include { BAM_TO_SAM                                 } from '../../modules/stenglein-lab/bam_to_sam'
include { TALLY_SAM_SUBJECTS                         } from '../../modules/stenglein-lab/tally_sam_subjects'
include { PREPEND_TSV_WITH_ID                        } from '../../modules/stenglein-lab/prepend_tsv_with_id'
include { PROCESS_SPECIES_OUTPUT                     } from '../../modules/stenglein-lab/process_species_output'

workflow SPECIES_ID {                                                    

  PREPROCESS_READS(params.fastq_dir, params.fastq_pattern, params.collapse_duplicate_reads)

  BUILD_BWA_INDEX(params.species_id_fasta)

  MAP_TO_REFERENCE(PREPROCESS_READS.out.reads, BUILD_BWA_INDEX.out.index)

  BAM_TO_SAM(MAP_TO_REFERENCE.out.bam.filter{it[1].size() > 0})

  TALLY_SAM_SUBJECTS(BAM_TO_SAM.out.sam)

  PREPEND_TSV_WITH_ID(TALLY_SAM_SUBJECTS.out.txt)

  PROCESS_SPECIES_OUTPUT(PREPEND_TSV_WITH_ID.out.tsv.collectFile(name: "all_species_tallies.txt"){it[1]}, params.metadata)
}

