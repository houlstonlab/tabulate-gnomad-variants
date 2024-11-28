#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { COORDINATES } from './modules/coordinates.nf'
include { SUBSET }      from './modules/subset.nf'
include { FILTER }      from './modules/filter.nf'
include { COMBINE }     from './modules/combine.nf'
include { EXTRACT }     from './modules/extract.nf'
include { AGGREGATE }   from './modules/aggregate.nf'
include { REPORT }      from './modules/report.nf'
include { CLINVAR }     from './modules/clinvar.nf'

// Define input channels
genes_coords_ch =  Channel.of (1 .. 22, 'X', 'Y')
    | map { [ "chr${it}", params.genome, params.style ] }

gnomad_ch = Channel.fromFilePairs(params.gnomad_files, flat: true)
    | map { it ->
        def key = it[0].split('\\.')[-1]
        [key, it[1], it[2]]
    }

clinvar_ch = Channel.fromFilePairs(params.clinvar_files, flat: true)

category_ch = Channel.of( 'Pathogenic', 'Damaging', 'Splicing',
                        //   'NotScored', 'Scored',
                          'High', 'PTV', 'Stop' )
clinvar_category_ch = Channel.of( 'Pathogenic', 'Benign' )
variables_ch = Channel.of( 'frequency' )

workflow  {
    // Extract clinvar variants
    genes_coords_ch
        | COORDINATES
        | combine(clinvar_ch)
        | combine(clinvar_category_ch)
        | CLINVAR
        | branch {
            pathogenic: it[1] == 'Pathogenic'
            benign: it[1] == 'Benign'
        }
        | set { clinvar }

    // // Filter gnomad variants
    gnomad_ch
        | combine(COORDINATES.out, by: 0)
        | SUBSET
        | combine(category_ch)
        | combine(clinvar.pathogenic, by: 0)
        | combine(clinvar.benign, by: 0)
        | FILTER
        | groupTuple(by: 1)
        | COMBINE
        | combine(variables_ch)
        | EXTRACT
        | multiMap { it ->
            cat: it
            all: ['ALL', it[1], it[2]]
        }
        | set { frequency }

    // Aggregate by category, and ALL
    frequency.cat 
        | concat(frequency.all)
        | filter { it[1] == 'frequency' }
        | groupTuple(by: [0,1])
        | AGGREGATE

    // Generate report
    COMBINE.out | REPORT   
}
