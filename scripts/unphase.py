import sys
from itertools import groupby
from pysam import VariantFile



vcfFile =VariantFile(sys.argv[1])
out=VariantFile(sys.argv[2],'w',header=vcfFile.header)


for rec in vcfFile.fetch():
    for sample in rec.samples.keys():
        rec.samples[sample].phased=False
        gt=list(rec.samples[sample]["GT"])
        if None not in gt:
            gt.sort()
        else:
            print(rec)
        rec.samples[sample]["GT"]=(gt[0],gt[1])
	
    out.write(rec)
    
