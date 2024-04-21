import sys
from itertools import groupby
from pysam import VariantFile



vcfFile =VariantFile(sys.argv[1])
minQual=int(sys.argv[2])
out=VariantFile(sys.argv[3],'w',header=vcfFile.header)


for rec in vcfFile.fetch():
    for sample in rec.samples.keys():    
        if rec.samples[sample]["GQ"] < minQual:
            rec.samples[sample]["GT"]=(None,None)
    out.write(rec)
