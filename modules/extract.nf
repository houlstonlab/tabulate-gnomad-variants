process EXTRACT {
    tag "${chrom}:${category}"

    label 'simple'

    container params.bcftools

    publishDir("${params.output_dir}/frequency", mode: 'copy')

    input:
    tuple val(chrom), val(category),
          path(file), path(index), path(variants),
          val(variable)

    output:
    tuple val(chrom), val(category), val(variable),
          path("${chrom}.${category}.${variable}.tsv")
 
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
            > ${chrom}.${category}.frequency.tsv
        """
    }
}
