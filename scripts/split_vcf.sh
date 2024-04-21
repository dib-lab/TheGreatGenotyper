#!/bin/bash

# Define the input VCF file
input_vcf=$1

# Define the maximum number of variants per file
max_variants=$2

# Define the output file prefix
output_prefix=$3

# Count the number of lines (variants) in the VCF file
num_variants=$(bcftools view -H $input_vcf | wc -l)

# Calculate the number of parts to split the VCF into
num_parts=$((num_variants / max_variants + 1))

# Extract header
bcftools view -h $input_vcf > header.$$.vcf

# Split the VCF file
for ((i=0; i<$num_parts; i++)); do
    start=$((i * max_variants + 1))
    end=$((start + max_variants - 1))
    output_vcf="${output_prefix}_split_${i}.vcf"
    bcftools view -H $input_vcf | sed -n ${start},${end}p | cat header.$$.vcf - > $output_vcf
done

# Remove the header file
rm header.$$.vcf
