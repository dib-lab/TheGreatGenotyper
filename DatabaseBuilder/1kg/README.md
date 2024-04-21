1. fix the resources in the Snakefile.build_2.smk
2. edit config.yaml
3. copy sample_table 
```
cp sample_tables_build/sample_table_all.csv sample_table.csv
```

4. create subsample_table.csv so that each sample has three files sample.fasta.gz, sample.kmer_count.gz, and sample.histo.

helpful commands
```
find {downloaded folder 1} |grep -P "(histo)$|(fasta.gz)$|(kmer_counts.gz)$" > files
find {downloaded folder 2} |grep -P "(histo)$|(fasta.gz)$|(kmer_counts.gz)$" >> files

grep -oP "/[^/]*$" files | cut -f1 -d'.' |tr -d '/' > tmp
echo -e "sample_name,file" > subsample_table.csv
paste -d, tmp files >> subsample_table.csv
```
take care that SGDP has smooth_1 and smooth_10000 contigs
and HGDP has onle smooth_100000