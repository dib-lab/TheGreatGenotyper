#!/bin/bash

# Input VCF file without samples
input_vcf=$1

# Output VCF file with added dummy sample
output_vcf=$2


bgzip -c ${input_vcf} > tmp.vcf.gz
tabix -p vcf tmp.vcf.gz

bcftools norm -m -any tmp.vcf.gz |  bcftools view -i 'ALT!="."' > tmp2.vcf


# Name of the dummy sample
sample_name="DummySample"

grep  '^##' "${input_vcf}" > "${output_vcf}"
# Add the sample header line to the output VCF file
echo -e "##SAMPLE=<ID=${sample_name},Genomes=1>" >> "${output_vcf}"
echo -e "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">" >> "${output_vcf}"

# Add the dummy sample genotype header line to the output VCF file
echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t${sample_name}" >> "${output_vcf}"

bgzip tmp2.vcf
tabix -p vcf tmp2.vcf.gz
# Process each variant in the input VCF file
bcftools view tmp2.vcf.gz   | grep -v '^#' | awk -v sample_name="${sample_name}" 'BEGIN{OFS="\t"}{print $1, $2, $3, $4, $5, $6, $7, $8, "GT", "0/1"}' >> "${output_vcf}"


rm tmp2.vcf.gz*
echo "Dummy sample added successfully to the VCF file."
