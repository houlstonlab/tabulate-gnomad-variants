process COMBINE {
    tag "${category}"

    label 'simple'

    container params.bcftools

    publishDir("${params.output_dir}/combined/", mode: 'copy')

    input:
    tuple val(chrom), val(category),
          path(file), path(index),
          path(variants)

    output:
    tuple val(category),
          path("${category}.vcf.gz"),
          path("${category}.vcf.gz.tbi"),
          path("${category}.tsv")
        
    script:
    """
    #!/bin/bash
    # Combine vcfs
    bcftools concat \
        --naive \
        ${file} \
        --threads ${task.cpu} \
        -Oz -o ${category}.vcf.gz
    tabix ${category}.vcf.gz

    # Combine variants
    cat ${variants} > ${category}.tsv
    """
}
