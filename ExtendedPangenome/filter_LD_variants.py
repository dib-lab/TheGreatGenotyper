import sys

for line in sys.stdin:
    l=line.strip().split("\t")
    allele_1,allele_2=l[22],l[26]
    if(len(allele_1) <15 and len(allele_2) <15):
        continue
    allele_1_id,allele_2_id=l[0],l[25]
    if "rs" not in allele_1_id and "rs" not in allele_2_id:
        continue
    INFO=l[5].split(";")
    INFO=[x.split("=") for x in INFO]
    MAF=list(filter(lambda x:x[0]=="MAF", INFO))[0]
    if float(MAF[1]) <0.01 :
        continue
    gene_position=l[16]
    gene_position_parts = gene_position.split("-")
    if "intron" in gene_position and gene_position_parts[0] ==gene_position_parts[1]:
        continue
    print(line.strip())
