include { FASTQC   as FASTQC_PRE           } from '../../modules/nf-core/fastqc/main'
include { FASTQC   as FASTQC_POST_TRIM     } from '../../modules/nf-core/fastqc/main'
include { FASTQC   as FASTQC_POST_COLLAPSE } from '../../modules/nf-core/fastqc/main'

include { MULTIQC                     } from '../../modules/nf-core/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../../modules/nf-core/custom/dumpsoftwareversions/main'

include { CHECK_FASTQ_COMPRESSED      } from '../../modules/stenglein-lab/check_fastq_compressed/main'
include { EXTRACT_READ_MAPPING_COUNTS } from '../../modules/stenglein-lab/extract_read_mapping_counts/main'
include { PREPEND_TSV_WITH_ID         } from '../../modules/stenglein-lab/prepend_tsv_with_id/main'
include { TRIM_READS                  } from '../../modules/stenglein-lab/cutadapt/main'
include { CD_HIT_EST                  } from '../../modules/stenglein-lab/cd_hit_est/main'

include { COUNT_FASTQ as COUNT_FASTQ_INITIAL       } from '../../modules/stenglein-lab/count_fastq/main'
include { COUNT_FASTQ as COUNT_FASTQ_POST_TRIM     } from '../../modules/stenglein-lab/count_fastq/main'
include { COUNT_FASTQ as COUNT_FASTQ_POST_COLLAPSE } from '../../modules/stenglein-lab/count_fastq/main'
include { COUNT_FASTQ as COUNT_FASTQ_POST_HOST     } from '../../modules/stenglein-lab/count_fastq/main'

include { PROCESS_FASTQ_COUNTS        } from '../../modules/local/process_fastq_counts/main'


workflow PREPROCESS_READS {

 main:

  // define some empty channels for keeping track of stuff
  ch_versions     = Channel.empty()                                               
  ch_fastq_counts = Channel.empty()                                               

  /*
   These fastq files are the main input to this workflow
  */
  Channel
  .fromFilePairs("${params.fastq_dir}/${params.fastq_pattern}", size: -1, checkIfExists: true, maxDepth: 1)
  .map{ name, reads ->

         // define a new empty map named meta for each sample
         // and populate it with id and single_end values
         // for compatibility with nf-core module expected parameters
         // reads are just the list of fastq
         def meta        = [:]

         // this map gets rid of any _S\d+ at the end of sample IDs but leaves fastq
         // names alone.  E.g. strip _S1 from the end of a sample ID..  
         // This is typically sample #s from Illumina basecalling.
         // could cause an issue if sequenced the same sample with 
         // multiple barcodes so was repeated on a sample sheet. 
         meta.id         = name.replaceFirst( /_S\d+$/ ,"")

         // if 2 fastq files then paired end data, so single_end is false
         meta.single_end = reads[1] ? false : true

         // this last statement in the map closure is the return value
         [ meta, reads ] }

  .set { ch_reads }

  CHECK_FASTQ_COMPRESSED( ch_reads ) 

  // count numbers of initial reads
  COUNT_FASTQ_INITIAL ( ch_reads.map{ meta, reads -> [ meta, reads, "initial"] } )
  ch_fastq_counts = ch_fastq_counts.mix(COUNT_FASTQ_INITIAL.out.count_file)
  ch_versions = ch_versions.mix ( COUNT_FASTQ_INITIAL.out.versions )      
  
  // run fastqc on input reads
  FASTQC_PRE ( ch_reads )
  ch_versions = ch_versions.mix ( FASTQC_PRE.out.versions )      

  // trim low quality bases and adapters from reads
  TRIM_READS ( ch_reads )
  ch_versions = ch_versions.mix ( TRIM_READS.out.versions ) 

  // count numbers of reads after trimming
  COUNT_FASTQ_POST_TRIM ( TRIM_READS.out.reads.map{ meta, reads -> [ meta, reads, "post_trimming"] } )
  ch_fastq_counts = ch_fastq_counts.mix(COUNT_FASTQ_POST_TRIM.out.count_file)

  // run fastqc on post trimmed reads
  FASTQC_POST_TRIM ( TRIM_READS.out.reads )

  ch_trimmed_reads = Channel.empty()
  ch_trimmed_reads = TRIM_READS.out.reads

  /*
  // collapse duplicates reads 
  CD_HIT_EST ( TRIM_READS.out.reads ) ch_versions = ch_versions.mix ( CD_HIT_EST.out.versions ) 

  // run fastqc on post collapsed reads
  FASTQC_POST_COLLAPSE ( CD_HIT_EST.out.reads )

  // count reads after collapsing
  COUNT_FASTQ_POST_COLLAPSE ( CD_HIT_EST.out.reads.map{ meta, reads -> [ meta, reads, "post_collapse"] } )
  ch_fastq_counts = ch_fastq_counts.mix(COUNT_FASTQ_POST_COLLAPSE.out.count_file)
  */

  //                                                                          
  // MODULE: Pipeline reporting                                               
  //                                                                          
  /*
  CUSTOM_DUMPSOFTWAREVERSIONS (                                               
      ch_versions.unique().collectFile(name: 'collated_versions.yml')         
  )    
  */

  //                                                                          
  // MODULE: MultiQC                                                          
  //                                                                          
  // workflow_summary    = WorkflowRnaseq.paramsSummaryMultiqc(workflow, summary_params)
  // ch_workflow_summary = Channel.value(workflow_summary)                   

  // collect files that will be input to multiqc
  ch_multiqc_files = Channel.empty()
  // ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
  ch_multiqc_files = ch_multiqc_files.mix(FASTQC_PRE.out.zip.collect{it[1]}.ifEmpty([]))
  ch_multiqc_files = ch_multiqc_files.mix(FASTQC_POST_TRIM.out.zip.collect{it[1]}.ifEmpty([]))
  ch_multiqc_files = ch_multiqc_files.mix(TRIM_READS.out.log.collect{it[1]}.ifEmpty([]))
  // ch_multiqc_files = ch_multiqc_files.mix(FASTQC_POST_COLLAPSE.out.zip.collect{it[1]}.ifEmpty([]))

  // run multiqc
  /*
  MULTIQC (
      ch_multiqc_files.collect(), [], [], []
  )
  multiqc_report = MULTIQC.out.report.toList()
  ch_versions    = ch_versions.mix(MULTIQC.out.versions)
  */


  /*
     Process all fastq counts and create PDF output
   */
  // PROCESS_FASTQ_COUNTS(ch_fastq_counts.collectFile(name: "all_fastq_counts.txt"))



 emit: 
  versions      = ch_versions
  multiqc_files = ch_multiqc_files 
  fastq_counts  = ch_fastq_counts 
  reads         = ch_trimmed_reads

}
