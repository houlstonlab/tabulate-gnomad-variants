process EXTRACT {
    tag "${category}:${variable}"

    label 'simple'

    container params.bcftools

    publishDir("${params.output_dir}/variants", mode: 'copy')

    input:
    tuple val(category),
          path(file), path(index), path(variants),
          val(variable)

    output:
    tuple val(category), val(variable),
          path("${category}.${variable}.tsv")
 
    script:
    if (variable == 'frequency') {
        """
        #!/bin/bash
        # Get frequency
        bcftools +split-vep \
            -a vep \
            -c SYMBOL \
            -s worst \
            -f '%SYMBOL\t%CHROM:%POS:%REF:%ALT\t%AC_nfe\t%AN_nfe\t%AF_nfe\t%nhomalt_nfe\n' \
            ${file} \
            > ${category}.${variable}.tsv
        """
    }
}
