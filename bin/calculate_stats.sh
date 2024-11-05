#!/bin/bash
chrom=$1
category=$2
file=$3

# Annotations: AF_grpmax AF_nfe spliceai_ds_max
# AF_grpmax=$( bcftools query -f "%AF_grpmax\n" ${file} | awk '{print $1 + 0}' | grep -v '^0$' | sort -uV | tail -1 )
# AF_nfe=$( bcftools query -f "%AF_grpmax\n" ${file} | awk '{print $1 + 0}' | grep -v '^0$' | sort -uV | tail -1 )
# spliceai_ds_max=$( bcftools query -f "%AF_grpmax\n" ${file} | awk '{print $1 + 0}' | grep -v '^0$' | sort -uV | tail -1 )

# Function to query bcftools and get the maximum value
get_val() {
    local query=$1
    local file=$2
    local which=$3

    bcftools query -f "$query" $file | \
        awk '{print $1 + 0}' | \
        grep -v '^0$' | \
        sort -ug > tmp.txt
    
    if [ "$which" == "max" ]; then
        tail -1 tmp.txt
    elif [ "$which" == "min" ]; then
        head -1 tmp.txt
    else
        echo "Invalid option. Please choose either 'max' or 'min'."
        exit 1
    fi
    rm tmp.txt
}

# Apply get_max function to the provided lines
AF_grpmax=$(get_val "%AF_grpmax\n" ${file} "max")
AF_nfe=$(get_val "%AF_nfe\n" ${file} "max")
spliceai_ds_max=$(get_val "%spliceai_ds_max\n" ${file} "min") 

# Print headers
echo -e "Chrom\tCategory\tAF_grpmax\tAF_nfe\tspliceai_ds_max"

# Print values
echo -e "$chrom\t$category\t$AF_grpmax\t$AF_nfe\t$spliceai_ds_max"
