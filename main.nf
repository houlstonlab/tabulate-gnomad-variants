#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { COORDINATES } from './modules/coordinates.nf'
include { SUBSET }      from './modules/subset.nf'
include { FILTER }      from './modules/filter.nf'
include { EXTRACT }     from './modules/extract.nf'
include { AGGREGATE }   from './modules/aggregate.nf'
include { REPORT }      from './modules/report.nf'
include { CLINVAR }     from './modules/clinvar.nf'

// Define input channels
genes_coords_ch =  Channel.of (1 .. 22, 'X', 'Y')
    | map { [ "chr${it}", params.genome, params.style ] }

gnomad_ch = Channel.fromFilePairs(params.gnomad_files, flat: true)

clinvar_ch = Channel.fromFilePairs(params.clinvar_files, flat: true)

category_ch = Channel.of( 'Pathogenic', 'Damaging', 'Splicing', 'NotScored', 'Scored' )
variables_ch = Channel.of( 'frequency' )

workflow  {
    // Extract clinvar variants
    genes_coords_ch
        | COORDINATES
        | combine(clinvar_ch)
        | combine(category_ch)
        | CLINVAR

    // Filter gnomad variants
    gnomad_ch
        | combine(COORDINATES.out)
        | SUBSET
        | combine(category_ch)
        | map { [it[0], it[3], it[1], it[2]] }
        | combine(CLINVAR.out, by: [0, 1])
        | FILTER
        | combine(variables_ch)
        | EXTRACT
        | multiMap { it ->
            cat: it
            all: [it[0], 'ALL', it[2], it[3]]
        }
        | set { frequency }

    // Aggregate by category, and ALL
    frequency.cat 
        | concat(frequency.all)
        | filter { it[2] == 'frequency' }
        | AGGREGATE

    // Generate report
    FILTER.out | REPORT   
    
    // Collect and store results
    AGGREGATE.out
        | concat(REPORT.out)
        | collectFile (
            keepHeader: true,
            storeDir: "${params.output_dir}/summary",
        )
        { it -> [ "gnomad.${it[1]}.${it[2]}.tsv", it.last() ] } 
}
