process SUBSET {
    tag "${chrom}"

    label 'simple'

    container params.bcftools

    publishDir("${params.output_dir}/subsets/", mode: 'copy')

    input:
    tuple val(chrom), path(gnomad_file), path(gnomad_index),
          path(gene_coords)

    output:
    tuple val(chrom),
          path("${chrom}.vcf.gz"),
          path("${chrom}.vcf.gz.tbi")
    
    script:
    """
    #!/bin/bash
    # Create bed file
    cat ${gene_coords} > region.bed
    cat ${gene_coords} | sed 's/^chr//' >> region.bed

    # Subset gnomad
    bcftools view -R region.bed ${gnomad_file} | \
    bcftools view -i 'AC_nfe > 0' | \
    bcftools annotate --set-id '%CHROM:%POS:%REF:%ALT' | \
    bcftools view --threads ${task.cpu} -Oz -o ${chrom}.vcf.gz

    # Index
    tabix ${chrom}.vcf.gz
    """
}
