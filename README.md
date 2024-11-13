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
    --gnomad_files input/gnomad.v4.vcf.gz{,.tbi} \
    --clinvar_files input/clinvar.20200520.vcf.gz{,.tbi} \
```

### Inputs & Parameters

- `gnomad_files`  : a VCF file with GnomAD variants
- `clinvar_files` : a VCF file with ClinVar variants
- `genome`: genome version (default is'hg38') 
- `style` : chromosome names style 'UCSC' or 'NCBI'
    
### Output

- `filtered/`  : filtered variants
- `frequency/` : the frequency of qualifying variants
- `aggregated/`: aggregated variant info by gene
- The final ouput in `summary/` is a tsv file with `gene`, `nvar`, `ac`, `an`, `af`, and `nhom` columns
