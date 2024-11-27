process AGGREGATE {
    tag "${category}"

    label 'simple'

    publishDir("${params.output_dir}/aggregate", mode: 'copy')

    input:
    tuple val(category), val(variable), path(file)

    output:
    tuple val(category), val("aggregate"),
          path("gnomad.${category}.aggregate.tsv")
 
    script:
    """
    #!/bin/bash
    echo -e "gene\tnvar\tac\tan\taf\tnhom" > gnomad.${category}.aggregate.tsv
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
    >> gnomad.${category}.aggregate.tsv
    """
}