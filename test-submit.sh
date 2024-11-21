#!/bin/bash

#SBATCH -o test/test.out
#SBATCH -e test/test.err
#SBATCH -J test
#SBATCH -p master-worker
#SBATCH -t 120:00:00

# Setup test directory
mkdir -p test/ test/input
cd test/

# Download test data
URL="https://figshare.com/ndownloader/files"

wget -c $URL/50357535 -O input/gnomad.v4.vcf.gz
wget -c $URL/50357526 -O input/gnomad.v4.vcf.gz.tbi
wget -c $URL/50357532 -O input/clinvar.20200520.vcf.gz
wget -c $URL/50357529 -O input/clinvar.20200520.vcf.gz.tbi

# Load required modules
module load Nextflow

# Run nextflow
# nextflow run houlstonlab/tabulate-gnomad-variants -r main \
nextflow run ../main.nf \
    --output_dir ./results/ \
    -profile local,gha \
    -resume

# usage: nextflow run [ local_dir/main.nf | git_url ]  
# These are the required arguments:
#     -r            {main,dev} to run specific branch
#     -profile      {local,cluster} to run using differens resources
#     -params-file  params.json to pass parameters to the pipeline
#     -resume       To resume the pipeline from the last checkpoint
