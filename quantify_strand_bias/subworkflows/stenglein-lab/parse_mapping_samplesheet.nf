//
// Retreive sample sheet that maps sampleIds to reference Fastas
//

// Function to get list of [ sample_id, refseq_fasta ]
def parse_samplesheet(LinkedHashMap row) {
    //Check if samplesheet contains required expected columns
    if (row.sampleID == null || row.referenceFasta == null) {
        exit 1, "ERROR: Please check input samplesheet -> Column 'sampleID' and 'referenceFasta' are required but not detected."
    }

    // iterate through sample sheet
    def array = []
    if (!file(row.referenceFasta).exists()) {
      exit 1, "ERROR: Please check input samplesheet -> Reference fasta file ${row.referenceFasta} does not exist.\n"
    }

    // make first element of output nf-core-style meta
    def meta        = [:]
    meta.id         = row.sampleID
    meta.single_end = true

    // return a list of sample_id and the path to the reference Fasta
    array = [ meta , file(row.referenceFasta), row.R1_strand ]

    return array
}

workflow PARSE_MAPPING_SAMPLESHEET {
    take:
    input // a tsv sample sheet containing two columns: sample_id and the path of a fasta file containing reference (virus) sequences 

    main:
      // extracts read files from TSV and distribute into channels
      Channel
        .fromPath(input, checkIfExists: true)
        .splitCsv(header:true, sep:'\t')
        .map { parse_samplesheet(it) }
        .set { ch_sample_sheet }

    emit:
    sample_sheet   = ch_sample_sheet
}
