[![Runs successfully](https://github.com/houlstonlab/tabulate-gnomad-variants/actions/workflows/runs-successfully.yml/badge.svg)](https://github.com/houlstonlab/tabulate-gnomad-variants/actions/workflows/runs-successfully.yml)

### Introduction

This workflow tests the filters and counts variants meeting certain criteria in GnomAD. Qualifying
variants are then tabulated by gene and information about thier frequency in the population is 
exported. Variants with different profiles can be extracted and checked against known variants in
ClinVar.

### Usage

The typical command looks like the following. `--gnomad_files` and `--clinvar_files` are required inputs. 
Different versions of the workflow can be called using `-r` and output directed to `--output_dir`

```bash
nextflow run houlstonlab/tabulate-gnomad-variants \
    -r main \
    --output_dir results/ \
    --gnomad_files input/cohorts_info.csv \
    --clinvar_files input/clinvar.20200520.vcf.gz{,.tbi} \
```

### Inputs & Parameters

- `gnomad_files`  : a csv file with three columns `chrom`, `file`, and `index`
- `clinvar_files` : a VCF file with ClinVar variants

Other paramters include:
- `genome`: genome version (default is'hg38') 
- `style` : chromosome names style 'UCSC' or 'NCBI'
- `categories`: selection categories. One or more of `'Pathogenic,Damaging,Splicing,High,PTV,Stop'`
- `AC`      : minimum allele count. Default `> 0`
- `AF_MAX`  : maximum group allele frequence. Default `> 0.005`
- `AF`      : maximum allele frequence. Default `> 0.005`
- `PAF_MAX` : maximum pathogenic variants group allele frequence. Default `> 0.01`
- `PAF`     : maximum pathogenic variants allele frequence. Default `> 0.01`
- `DS`      : minimum delta score for spliceAI
  
### Output

- `clinvar/`    : clinvar annotations
- `coordiantes/`: coding gene coordinates
- `subsets/`    : split vcf file by chromosome
- `filtered/`   : filtered vcfs
- `combined/`   : a combined vcf
- `variants/`   : extracted variant frequence
- `aggregate/`  : aggregated variant frequency by gene
- `reports/`    : summary reports
