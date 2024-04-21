#!/bin/bash
outputFolder=$1

# Define the list of populations
populations="Africa.1 America.1 EUR.FIN East_Asia.CHB East_Asia.mix South_Asia.PJL Africa.2 America.3 EUR.GBR East_Asia.CHS South_Asia.1 South_Asia.STU Africa.ACB_PUR America.4 EUR.IBS East_Asia.China South_Asia.BEB West_Eurasia Africa.ASW_MSL America_Oceanic EUR.TSI East_Asia.JPT South_Asia.GIH middle_east Africa_middle_east EUR.CEU East_Asia.CDX East_Asia.KHV South_Asia.ITU"

# Loop through each population
for population in $populations; do
  mkdir -p $outputFolder/$population 
  curl -o $outputFolder/$population/graph.desc.tsv https://farm.cse.ucdavis.edu/~tahmed/GG_index/$population/graph.desc.tsv
  curl -o $outputFolder/$population/graph.dbg https://farm.cse.ucdavis.edu/~tahmed/GG_index/$population/graph.dbg 
  curl -o $outputFolder/$population/annotation.relaxed.row_diff_int_brwt.annodbg https://farm.cse.ucdavis.edu/~tahmed/GG_index/$population/annotation.relaxed.row_diff_int_brwt.annodbg 
done
