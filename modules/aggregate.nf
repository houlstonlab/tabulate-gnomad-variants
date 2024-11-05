process AGGREGATE {
    tag "${category}"

    label 'simple'

    publishDir("${params.output_dir}/aggregated", mode: 'copy')

    input:
    tuple val(chrom), val(category), val(variable), path(file)

    output:
    tuple val(chrom), val(category), val("aggregate"),
          path("${chrom}.${category}.aggregate.tsv")
 
    script:
    """
    #!/bin/bash
    echo -e "gene\tnvar\tac\tan\taf\tnhom" > ${chrom}.${category}.aggregate.tsv
    cat ${file} | \
    sort -u | \
    awk '
    {
        key = \$1
        count[key][\$2]++
        for (i = 3; i <= NF; i++) {
            sum[key][i] += \$i
        }
    }
    END {
        for (key in sum) {
            printf "%s\t", key
            unique_count = length(count[key])
            printf "%d\t", unique_count
            for (i = 3; i <= NF; i++) {
                if (sum[key][i] == int(sum[key][i])) {
                    printf "%d\t", sum[key][i]
                } else {
                    printf "%f\t", sum[key][i]
                }
            }
            printf "\\n"
        }
    }
    ' \
    >> ${chrom}.${category}.aggregate.tsv
    """
}