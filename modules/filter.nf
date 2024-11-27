process FILTER {
    tag "${chrom}:${category}"

    label 'simple'

    container params.bcftools

    publishDir("${params.output_dir}/filtered/", mode: 'copy')

    input:
    tuple val(chrom), path(gnomad_file), path(gnomad_index),
          val(category),
          val(pathogenic), path(pathogenic_file),
          val(benign), path(benign_file)
          

    output:
    tuple val(chrom), val(category),
          path("${chrom}.${category}.vcf.gz"),
          path("${chrom}.${category}.vcf.gz.tbi"),
          path("${chrom}.${category}.tsv")
    
    script:
    if ( category == 'Pathogenic' ) {
        """
        #!/bin/bash
        # Filter gnomad
        bcftools view -e 'AF_grpmax > 0.01 || AF_nfe > 0.01' ${gnomad_file} | \
        bcftools view -i ID==@${pathogenic_file} | \
        bcftools view --threads ${task.cpu} -Oz -o ${chrom}.${category}.vcf.gz

        tabix ${chrom}.${category}.vcf.gz

        bcftools query -f '%ID\n' ${chrom}.${category}.vcf.gz >${chrom}.${category}.tsv
        """
    } else if ( category == 'Damaging' ) {
        """
        #!/bin/bash
        # Filter gnomad
        bcftools view -e 'AF_grpmax > 0.005 || AF_nfe > 0.005' ${gnomad_file} | \
        bcftools +split-vep -a vep -s worst -c IMPACT,LoF | \
        bcftools view -i 'IMPACT="HIGH"' | \
        bcftools view -i 'LoF="HC"' | \
        bcftools view -e ID==@${benign_file} | \
        bcftools view --threads ${task.cpu} -Oz -o ${chrom}.${category}.vcf.gz

        tabix ${chrom}.${category}.vcf.gz

        bcftools query -f '%ID\n' ${chrom}.${category}.vcf.gz > ${chrom}.${category}.tsv
        """
    } else if ( category == 'High' ) {
        """
        #!/bin/bash
        # Filter gnomad
        bcftools view -e 'AF_grpmax > 0.005 || AF_nfe > 0.005' ${gnomad_file} | \
        bcftools +split-vep -a vep -s worst -c IMPACT,LoF | \
        bcftools view -i 'IMPACT="HIGH"' | \
        bcftools view -e ID==@${benign_file} | \
        bcftools view --threads ${task.cpu} -Oz -o ${chrom}.${category}.vcf.gz

        tabix ${chrom}.${category}.vcf.gz

        bcftools query -f '%ID\n' ${chrom}.${category}.vcf.gz > ${chrom}.${category}.tsv
        """
    } else if ( category == 'Stop' ) {
        """
        #!/bin/bash
        # Filter gnomad
        bcftools view -e 'AF_grpmax > 0.005 || AF_nfe > 0.005' ${gnomad_file} | \
        bcftools +split-vep -a vep -s worst -c Consequence,IMPACT,LoF | \
        bcftools view -i 'Consequence~"stop_gained"' | \
        bcftools view -e ID==@${benign_file} | \
        bcftools view --threads ${task.cpu} -Oz -o ${chrom}.${category}.vcf.gz

        tabix ${chrom}.${category}.vcf.gz

        bcftools query -f '%ID\n' ${chrom}.${category}.vcf.gz > ${chrom}.${category}.tsv
        """
    } else if ( category == 'PTV' ) {
        """
        #!/bin/bash
        # Filter gnomad
        bcftools view -e 'AF_grpmax > 0.005 || AF_nfe > 0.005' ${gnomad_file} | \
        bcftools +split-vep -a vep -s worst -c Consequence,IMPACT,LoF | \
        bcftools view -i 'Consequence~"stop_gained" || Consequence~"frameshift_variant" || Consequence~"splice_acceptor_variant"' | \
        bcftools view -e ID==@${benign_file} | \
        bcftools view --threads ${task.cpu} -Oz -o ${chrom}.${category}.vcf.gz

        tabix ${chrom}.${category}.vcf.gz

        bcftools query -f '%ID\n' ${chrom}.${category}.vcf.gz > ${chrom}.${category}.tsv
        """
    } else if ( category == 'Scored' ) {
        """
        #!/bin/bash
        # Filter gnomad
        bcftools view -e 'AF_grpmax > 0.005 || AF_nfe > 0.005' ${gnomad_file} | \
        bcftools +split-vep -a vep -s worst -c IMPACT,LoF | \
        bcftools view -i 'IMPACT="HIGH" || IMPACT="MODERATE"' | \
        bcftools view -i 'LoF="HC"' | \
        bcftools view -e ID==@${benign_file} | \
        bcftools view --threads ${task.cpu} -Oz -o ${chrom}.${category}.vcf.gz

        tabix ${chrom}.${category}.vcf.gz

        bcftools query -f '%ID\n' ${chrom}.${category}.vcf.gz > ${chrom}.${category}.tsv
        """
    } else if ( category == 'NotScored' ) {
        """
        #!/bin/bash
        # Filter gnomad
        bcftools view -e 'AF_grpmax > 0.005 || AF_nfe > 0.005' ${gnomad_file} | \
        bcftools +split-vep -a vep -s worst -c IMPACT,LoF | \
        bcftools view -i 'IMPACT="HIGH" || IMPACT="MODERATE"' | \
        bcftools view -e ID==@${benign_file} | \
        bcftools view --threads ${task.cpu} -Oz -o ${chrom}.${category}.vcf.gz

        tabix ${chrom}.${category}.vcf.gz

        bcftools query -f '%ID\n' ${chrom}.${category}.vcf.gz > ${chrom}.${category}.tsv
        """
    } else if ( category == 'Splicing' ) {
        """
        #!/bin/bash
        # Filter gnomad
        bcftools view -e 'AF_grpmax > 0.005 || AF_nfe > 0.005' ${gnomad_file} | \
        bcftools view -i 'spliceai_ds_max > 0.8'  | \
        bcftools view -e ID==@${benign_file} | \
        bcftools view --threads ${task.cpu} -Oz -o ${chrom}.${category}.vcf.gz

        tabix ${chrom}.${category}.vcf.gz

        bcftools query -f '%ID\n' ${chrom}.${category}.vcf.gz > ${chrom}.${category}.tsv
        """
    }
}
