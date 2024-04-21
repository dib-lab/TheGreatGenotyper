import sys
from itertools import groupby
from pysam import VariantFile

def fasta_iter(fasta_name):
    fh = open(fasta_name)

    # ditch the boolean (x[0]) and just keep the header or sequence since
    # we know they alternate.
    faiter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))

    for header in faiter:
        # drop the ">"
        headerStr = header.__next__()[1:].strip()
        # join all sequence lines to one.
        seq = "".join(s.strip() for s in faiter.__next__())
        yield (headerStr, seq)



vcfFile =VariantFile(sys.argv[1])
refIter=fasta_iter(sys.argv[2])
out=VariantFile(sys.argv[3],'w',header=vcfFile.header)

reference={} 

for header,seq in refIter:
    reference[header.split(" ")[0]]=seq


for rec in vcfFile.fetch():
    if "SVTYPE" in rec.info and rec.info["SVTYPE"] =="DEL" and rec.alts == ('<DEL>',):
        rec.alts=rec.ref[0]
        del_len=abs(rec.info["SVLEN"][0])
#        print(rec.start,rec.stop,)
        rec.ref = reference[rec.chrom][rec.start:rec.start+del_len]
    if "SVTYPE" in rec.info and rec.info["SVTYPE"] =="INS" and rec.alts == ('<INS>',):
        left=rec.info["LEFT_SVINSSEQ"][0]
        right=rec.info["RIGHT_SVINSSEQ"][0]
        print(left,right)
        rec.alts=(left + "N" +right,)        
    out.write(rec)
    
