#!/bin/bash

# Input VCF file without samples
input_vcf="clinvar.vcf"

# Output VCF file with added dummy sample
output_vcf="clinvar.normalized.noalt.dummySamples.vcf"


bgzip -c ${input_vcf} > tmp.vcf.gz
tabix -p vcf tmp.vcf.gz

bcftools norm -m +any tmp.vcf.gz |  bcftools view -i 'ALT!="."' > tmp2.vcf



bgzip tmp2.vcf
tabix -p vcf tmp2.vcf.gz

bcftools view tmp2.vcf.gz 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X  > tmp.3.vcf

bcftools annotate --force --rename-chrs chr_name_conv.txt tmp.3.vcf >  tmp.4.vcf
 
python addSample_clivar.py tmp.4.vcf $output_vcf


rm tmp.4.vcf tmp.3.vcf tmp2.vcf.gz tmp2.vcf.gz.tbi tmp.vcf.gz* 
echo "Dummy sample added successfully to the VCF file."
