output=$1
input="${@:2}"
bcftools merge $input |sed -e 's/GT:[^\t]*/GT/' |sed -e 's/\([0-9]\/[0-9]\)[^\t]*/\1/g'|sed -e 's/.:.:[^\t]*/0\/0/g' |  bcftools +fill-tags - -Ov  -- -t all | bgzip > $output

tabix -p vcf $output
