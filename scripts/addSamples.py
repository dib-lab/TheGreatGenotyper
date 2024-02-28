import sys
from itertools import groupby
from pysam import VariantFile,VariantHeader
import math

vcfFile =VariantFile(sys.argv[1])
outVCF= sys.argv[2]
maxAlleles=1

for rec in vcfFile.fetch():
    maxAlleles=max(maxAlleles,len(rec.alts))





numSamples= int(math.ceil(float(maxAlleles+1)/2.0))
print(f"Max Alleles = {maxAlleles} and number of samples = {numSamples}")
vcfFile2 =VariantFile(sys.argv[1])


new_header = VariantHeader()
for record in vcfFile2.header.records:
    if record.type != 'SAMPLE':
        new_header.add_record(record)
        



for s in range(0,numSamples):
    new_header.samples.add(f"Sample_{s}")

    
out=VariantFile(outVCF,'w',header=new_header)



for old_record in vcfFile2.fetch():
#    record= record.translate(out.header)
    record= out.new_record(contig=old_record.contig, start=old_record.start, stop=old_record.stop, alleles=old_record.alleles, id=old_record.id, qual=old_record.qual, filter=old_record.filter, info=old_record.info)
    currentAllele=0
    numAlleles=len(old_record.alleles)
    for s in range(0,numSamples):
        sample=f"Sample_{s}"
        a=currentAllele
        currentAllele= (currentAllele +1) % numAlleles 
        b=currentAllele
        currentAllele= (currentAllele +1) % numAlleles
        record.samples[sample]["GT"] = (a, b)
#        record.samples[sample]["GT"].phased = True

    out.write(record)
        
