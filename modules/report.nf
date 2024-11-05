process REPORT {
    tag "${chrom}:${category}"

    label 'simple'

    container params.bcftools

    publishDir("${params.output_dir}/reports", mode: 'copy')

    input:
    tuple val(chrom), val(category), 
          path(file), path(index), path(variants)

    output:
    tuple val(chrom), val(category), val("report"),
          path("${chrom}.${category}.report.tsv")

    script:
    """
    #!/bin/bash
    calculate_stats.sh \
        ${chrom} ${category} \
        ${file} \
        > ${chrom}.${category}.report.tsv
    """
}
