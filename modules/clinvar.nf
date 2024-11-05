process CLINVAR {
    tag "${chrom}:${category}"

    label 'simple'

    container params.bcftools

    publishDir("${params.output_dir}/clinvar/", mode: 'copy')

    input:
    tuple val(chrom), path(gene_coords),
          val(clinvar), path(clinvar_file), path(clinvar_index),
          val(category)

    output:
    tuple val(chrom), val(category),
          path("${chrom}.${category}.txt")
    
    script:
    if ( category == 'Pathogenic' ) {
        """
        #!/bin/bash
        # Create bed file
        cat ${gene_coords} > region.bed
        cat ${gene_coords} | sed 's/^chr//' >> region.bed

        # Filter clinvar
        bcftools view \
            -R region.bed \
            -i 'ORIGIN="1" && (CLNSIG~"Pathogenic" || CLNSIG~"Likely_pathogenic")' \
            ${clinvar_file} | \
        bcftools annotate \
            --set-id '%CHROM:%POS:%REF:%ALT' | \
        bcftools view \
            -Oz -o ${chrom}.${category}.vcf.gz

        tabix ${chrom}.${category}.vcf.gz

        bcftools query -f 'chr%ID\n' ${chrom}.${category}.vcf.gz > ${chrom}.${category}.txt
        """
    } else if ( category == 'Damaging' || category == 'Splicing' ) {
        """
        #!/bin/bash
        # Create bed file
        cat ${gene_coords} > region.bed
        cat ${gene_coords} | sed 's/^chr//' >> region.bed

        # Filter clinvar
        bcftools view \
            -R region.bed \
            -i 'ORIGIN="1" && (CLNSIG~"Benign" || CLNSIG~"Likely_benign")' \
            ${clinvar_file} | \
        bcftools annotate \
            --set-id '%CHROM:%POS:%REF:%ALT' | \
        bcftools view \
            -Oz -o ${chrom}.${category}.vcf.gz

        tabix ${chrom}.${category}.vcf.gz

        bcftools query -f 'chr%ID\n' ${chrom}.${category}.vcf.gz > ${chrom}.${category}.txt
        """
    }
}
