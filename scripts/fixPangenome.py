import sys
from itertools import groupby
from pysam import VariantFile



vcfFile =VariantFile(sys.argv[1])
out=VariantFile(sys.argv[2],'w',header=vcfFile.header)


for rec in vcfFile.fetch():
    for sample in rec.samples.keys():    
        if None in list(rec.samples[sample]["GT"]):
            rec.samples[sample].phased=False
            
    out.write(rec)
