import sys
from itertools import groupby
from pysam import VariantFile



vcfFile =VariantFile(sys.argv[1])
#minQual=int(sys.argv[2])
out=VariantFile(sys.argv[2],'w',header=vcfFile.header)

HWE_threshold=10 ** -50
EXC_threshold= 0.1
for rec in vcfFile.fetch():
    stats={}
    qs={}
    qualThreshold={}
    # if rec.info['ExcHet'][0] < EXC_threshold:
    #     for sample in rec.samples.keys():
    #         rec.samples[sample]["GT"]=(None,None)
    #     out.write(rec)
    #     continue

    for sample in rec.samples.keys():
        gt=rec.samples[sample]["GT"]
        if gt not in stats:
            stats[gt]={}
            qs[gt]=[]
            qualThreshold[gt]=0
        gq=(rec.samples[sample]["GQ"])
        if gq == None:
            gq=0
        qs[gt].append(gq)
        # if gq not in stats[gt]:
        #     stats[gt][gq]=0
        # stats[gt][gq]+=1

    for k,v in stats.items():
        freqs=list(v.items())
        freqs=sorted(freqs,key=lambda x: -x[1] )
        qs[k].sort()
        qualThreshold[k] = qs[k][int(len(qs[k])/2)]
        
    for sample in rec.samples.keys():
        gt=rec.samples[sample]["GT"]
        gq=(rec.samples[sample]["GQ"])    
        if gq ==None or gq < qualThreshold[gt]:
            rec.samples[sample]["GT"]=(None,None)

    out.write(rec)
  
