//
// Retreive sample sheet that maps sampleIds to reference Fastas
//

// Function to get list of [ sample_id, R1_orientation ]
def parse_samplesheet(LinkedHashMap row) {
    //Check if samplesheet contains required expected columns
    if (row.sampleID == null || row.R1_strand == null) {
        println("row: " + row)
        exit 1, "ERROR: Please check input samplesheet -> Column 'sampleID' and 'R1_strand' are required but not detected."
    }

    // iterate through sample sheet
    def array = []

    // make first element of output nf-core-style meta
    def meta        = [:]
    meta.id         = row.sampleID
    // this pipeline designed to work with single end datasets
    meta.single_end = true

    // return a list of sample_id and the path to the reference Fasta
    array = [ meta , row.R1_strand ]

    return array
}

workflow PARSE_MAPPING_SAMPLESHEET {
    take:
    input // a tsv sample sheet containing two columns: sample_id and the orientation of the R1 read (for strand-specific lib prep)

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
