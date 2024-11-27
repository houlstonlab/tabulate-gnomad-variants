process REPORT {
    tag "${chrom}:${category}"

    label 'simple'

    container params.bcftools

    publishDir("${params.output_dir}/reports", mode: 'copy')

    input:
    tuple val(category), 
          path(file), path(index), path(variants)

    output:
    tuple val(category), val("report"),
          path("${category}.report.tsv")

    script:
    """
    #!/bin/bash
    calculate_stats.sh \
        ${category} \
        ${file} \
        > ${category}.report.tsv
    """
}
