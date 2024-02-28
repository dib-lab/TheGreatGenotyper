import os
import glob
import random




input_folders=["/group/ctbrowngrp/mshokrof/1kg_clinvar/byChromosome/" , "/group/ctbrowngrp/mshokrof/1kg_GWAS/byChromosome/","/group/ctbrowngrp/mshokrof/pangenome_GG_GRCh38/merged/","/group/ctbrowngrp/mshokrof/1kg_GG_GRCh38/merged_plink_filtered/"]
tempFolder = "/scratch/mshokrof/"
input_reference = "/group/ctbrowngrp/mshokrof/GRCh38_full_analysis_set_plus_decoy_hla.fa"
beagle = "/home/mshokrof/genotyping-experiments/SV_genotyping/beagle.22Jul22.46e.jar"
beagleMap =  "/home/mshokrof/genotyping-experiments/SV_genotyping/plink.autsomal.map"
OUTPUT_DIR = "/group/ctbrowngrp/mshokrof/Extended_pangenome/"




## Get the list of chromosomes based on the input files
CHROMOSOMES = ["chr"+str(x) for x in range(1,23)] +["chrX"]
CHROMOSOMES = ["chr"+str(x) for x in range(1,23)] 
#CHROMOSOMES = [ "chr21"]


rule all:
    input:
        expand(os.path.join(OUTPUT_DIR, "annotate_sv/{chrom}.vcor"), chrom=CHROMOSOMES)




ruleorder: filltags > run_beagle 


# localrules: indexVCF


# rule indexVCF:
#     input:
#         vcf="{prefix}.vcf.gz",
#     output:
#         index="{prefix}.vcf.gz.tbi"
#     shell:
#         """
# 	    tabix -p vcf {input.vcf}
#         """


rule fix_ids_for_jasmine:
    input:
        OUTPUT_DIR+ "split/{chrom}_{id}.sv.vcf"
    output:
        OUTPUT_DIR+ "split/{chrom}_idfixed_{id}.sv.vcf"
    log:
        OUTPUT_DIR+ "jasmine/{chrom}_idfixed_{id}.log"
    threads: 1
    retries: 0
    resources:
        mem_mb=4 *1024,
        cores=32,
        nodes = 1,
        meduim=1,
        runtime = 60 * 16,
        partition =  "med2",
        tmpdir= lambda wildcards: "%sjasmine%s/"%(tempFolder,f"{wildcards.chrom}")
    conda:
        "jasmine_env.yaml"
    shell:
        """
        perl fix_id.pl < {input} > {output}
        """



rule run_jasmine:
    input:
        expand(OUTPUT_DIR+ "split/{{chrom}}_idfixed_{id}.sv.vcf", id = range(len(input_folders)) )
    output:
        temp(OUTPUT_DIR+ "jasmine/{chrom}.bcf")
    log:
        OUTPUT_DIR+ "jasmine/{chrom}.log"
    threads: 16
    retries: 0
    resources:
        mem_mb=80 *1024,
        cores=32,
        nodes = 1,
        meduim=1,
        runtime = 60 * 16,
        partition =  "med2",
        tmpdir= lambda wildcards: "%sjasmine%s/"%(tempFolder,f"{wildcards.chrom}")
    conda:
        "jasmine_env.yaml"
    shell:
        """
        mkdir -p {resources.tmpdir}
        rm -rf {resources.tmpdir}*
        echo {input} |  tr -s ' ' $'\n'  > {resources.tmpdir}jasmine_list.txt

        BINDIR=$(which jasmine| sed -e 's/\/jasmine$//')
        java -Xmx80G -cp $BINDIR/jasmine_iris.jar:$BINDIR/jasmine.jar  Main iris_args=samtools_path=samtools,racon_path=racon,minimap_path=minimap2 file_list={resources.tmpdir}jasmine_list.txt out_file={resources.tmpdir}out.vcf threads={threads}  2> {log}
        
        bcftools view {resources.tmpdir}out.vcf -O b -o {output}

        rm -rf {resources.tmpdir}
        """


        

# bcftools concat {resources.tmpdir}out.vcf {resources.tmpdir}merged.small.vcf | \
#             bcftools sort -T {resources.tmpdir}  --max-mem 20G | \
#             bcftools norm -m +any | \
#             bgzip -c  > {output}  

rule split_to_sv_small:
    input:
        lambda wildcards: input_folders[int(f"{wildcards.id}")]+ "{chrom}.vcf.gz"
    output:
        small= temp(OUTPUT_DIR+ "split/{chrom}_{id}.small.bcf"),
        sv = temp(OUTPUT_DIR+ "split/{chrom}_{id}.sv.vcf")
    log:
        OUTPUT_DIR+ "split/{chrom}_{id}.log"
    threads: 16
    retries: 0
    resources:
        mem_mb=10 *1024,
        cores=32,
        nodes = 1,
        meduim=1,
        runtime = 60 * 16,
        partition =  "med2",
        tmpdir= lambda wildcards: "%ssplit%s/"%(tempFolder,f"{wildcards.chrom}")
    conda:
        "jasmine_env.yaml"
    shell:
        """
        mkdir -p {resources.tmpdir}
        rm -rf {resources.tmpdir}*

        bcftools query -l  {input} | sort > {resources.tmpdir}samples.txt

        mkfifo  {resources.tmpdir}sv.fifo

        cat {resources.tmpdir}sv.fifo | \
            bcftools view -i 'INFO/VARLEN >= 30' > {output.sv} &

        bcftools view -S {resources.tmpdir}samples.txt {input} | \
        bcftools norm -m -any  | \
            python ~/TheGreatGenotyper/scripts/assign-variant-len.py   | \
            tee {resources.tmpdir}sv.fifo | \
            bcftools view -i 'INFO/VARLEN < 30' -O b -o {output.small} 
        rm -rf {resources.tmpdir}
        """


rule merge_all:
    input:
        small = expand(OUTPUT_DIR+ "split/{{chrom}}_{id}.small.bcf", id = range(0,len(input_folders)) ),
        sv = OUTPUT_DIR+ "jasmine/{chrom}.bcf"
    output:
        OUTPUT_DIR+ "merged/{chrom}.vcf.gz"
    log:
        OUTPUT_DIR+ "merged/{chrom}.log"
    threads: 16
    retries: 0
    resources:
        mem_mb=40 *1024,
        cores=32,
        nodes = 1,
        meduim=1,
        runtime = 60 * 16,
        partition =  "med2",
        tmpdir= lambda wildcards: "%smerge_all%s/"%(tempFolder,f"{wildcards.chrom}")
    conda:
        "jasmine_env.yaml"
    shell:
        """
        mkdir -p {resources.tmpdir}
        bcftools concat -Ou {input.small} {input.sv} | \
            bcftools sort -T  {resources.tmpdir} -m 30G  -Ou | \
            bcftools norm -m +any -Oz -o  {output}
            tabix -p vcf {output}
        rm -rf {resources.tmpdir}
        """


 


rule run_beagle:
    input:
        os.path.join(OUTPUT_DIR, "merged/{chrom}.vcf.gz")
    output:
        vcf=OUTPUT_DIR+ "beagle/{chrom}.vcf.gz"
    params:
        prefix=OUTPUT_DIR+ "beagle/{chrom}"
    log:
        OUTPUT_DIR+ "beagle/{chrom}.log"
    threads: 64
    retries: 0
    resources:
        mem_mb=110 *1024,
        cores=32,
        nodes = 1,
        meduim=1,
        runtime = 60 * 24,
        partition =  "med2",
        tmpdir= lambda wildcards: "%sbeagle%s/"%(tempFolder,f"{wildcards.chrom}")
    shell:
        """
        mkdir -p {resources.tmpdir}
        bcftools norm -m -any -Oz -o  {resources.tmpdir}ref.vcf.gz {input} 
        tabix -p vcf {resources.tmpdir}ref.vcf.gz
        java -Xmx100G -jar {beagle} gp=true ap=true gt={resources.tmpdir}ref.vcf.gz out={params.prefix} nthreads={threads}  map={beagleMap} &> {log}
        rm -rf {resources.tmpdir}
        """

rule fix_beagle_header:
    input:
        vcf=OUTPUT_DIR+ "beagle/{chrom}.vcf.gz"
    output:
        vcf=OUTPUT_DIR+ "beagle/{chrom}.reheader.vcf.gz",
        tbi= OUTPUT_DIR+ "beagle/{chrom}.reheader.vcf.gz.tbi"
    log:
        OUTPUT_DIR+ "beagle/{chrom}.reheader.log"
    threads: 1
    retries: 0
    resources:
        mem_mb=1 *1024,
        cores=32,
        nodes = 1,
        meduim=1,
        runtime = 60 * 1,
        partition =  "med2",
        tmpdir= lambda wildcards: "%sbeagle_reheader%s/"%(tempFolder,f"{wildcards.chrom}")
    shell:
        """
        mkdir -p {resources.tmpdir}
        bcftools view  --header-only {input.vcf} |perl fixHeader.pl > {resources.tmpdir}header
        bcftools reheader -h  {resources.tmpdir}header {input.vcf}|bgzip -dc   |python unique_ID.py|bgzip >  {output.vcf}
	    #tabix -p vcf {resources.tmpdir}ref.beagle.vcf.gz 
	    #bcftools norm -m +any -Oz -o {output.vcf} {resources.tmpdir}ref.beagle.vcf.gz
        tabix -p vcf {output.vcf}
        rm -rf {resources.tmpdir}
        """





rule filltags:
    input:
        vcf="{OUTPUT_DIR}/beagle/{chrom}.reheader.vcf.gz",
	    index="{OUTPUT_DIR}/beagle/{chrom}.reheader.vcf.gz.tbi"
    output:
        vcf="{OUTPUT_DIR}/beagle/{chrom}.tagged.vcf.gz",
        tbi="{OUTPUT_DIR}/beagle/{chrom}.tagged.vcf.gz.tbi",
    threads: 4
    resources:
        mem_mb=10 *1024,
        cores=4,
        nodes = 1,
        meduim=1,
        runtime = 60 * 8,
        partition =  "med2",
    shell:
        """
        bcftools +fill-tags {input.vcf} -Ov -- -t all | bgzip > {output.vcf}
    	tabix -p vcf {output.vcf}
        """


rule annotSV:
    input:
        vcf="{OUTPUT_DIR}/beagle/{chrom}.tagged.vcf.gz",
	    index="{OUTPUT_DIR}/beagle/{chrom}.tagged.vcf.gz.tbi"
    output:
        annotated_tsv= "{OUTPUT_DIR}/annotate_sv/{chrom}.annotated.tsv",
        unannotated_tsv= "{OUTPUT_DIR}/annotate_sv/{chrom}.unannotated.tsv"
    threads: 4
    resources:
        mem_mb=50 *1024,
        cores=4,
        nodes = 1,
        meduim=1,
        runtime = 60 * 10,
        partition =  "med2",
    shell:
        """
        mkdir -p Annot_SV_{wildcards.chrom}
        bcftools view -s HGDP00845 {input.vcf}| bgzip > Annot_SV_{wildcards.chrom}/chr21.vcf.gz
        tabix -p vcf Annot_SV_{wildcards.chrom}/chr21.vcf.gz
        ~/AnnotSV/bin/AnnotSV  -SVminSize 20 -vcf 0 -SVinputFile Annot_SV_{wildcards.chrom}/chr21.vcf.gz -outputDir Annot_SV_{wildcards.chrom} -benignAF 0.001
        cp Annot_SV_{wildcards.chrom}/chr21.annotated.tsv {output.annotated_tsv}
        cp Annot_SV_{wildcards.chrom}/chr21.unannotated.tsv {output.unannotated_tsv}
        """

rule filterAffectingGenes:
    input:
        annotated_tsv= "{OUTPUT_DIR}/annotate_sv/{chrom}.annotated.tsv"
    output:
        annotated_tsv= "{OUTPUT_DIR}/annotate_sv/{chrom}.annotated.affecting_genes.tsv"
    threads: 1
    resources:
        mem_mb=1 *1024,
        cores=1,
        nodes = 1,
        meduim=1,
        runtime = 60 * 1,
        partition =  "med2",
    shell:
        """
        cut -f8,2,3,4,5,13,18-29,34,108 {input} |awk '{{if(NF>=12){{print $0}}}}' >{output}
        """

rule calulateLD:
    input:
        vcf="{OUTPUT_DIR}/beagle/{chrom}.tagged.vcf.gz",
	    index="{OUTPUT_DIR}/beagle/{chrom}.tagged.vcf.gz.tbi",
        annotated_tsv= "{OUTPUT_DIR}/annotate_sv/{chrom}.annotated.affecting_genes.tsv"
    output:
        vcor= "{OUTPUT_DIR}/annotate_sv/{chrom}.vcor"
    params:
        out_prefix= "{OUTPUT_DIR}/annotate_sv/{chrom}"
    threads: 8
    resources:
        mem_mb=60 *1024,
        cores=8,
        nodes = 1,
        meduim=1,
        runtime = 60 * 10,
        partition =  "med2",
        tmpdir= lambda wildcards: "%scalculateLD%s/"%(tempFolder,f"{wildcards.chrom}")
    shell:
        """
        mkdir -p {resources.tmpdir}
        tail -n +2  {input.annotated_tsv} |cut -f5 |sort |uniq  > {wildcards.chrom}.sv.ids

        ./plink2 --indep-pairwise  100kb 0.8 \
                 --indep-preferred {wildcards.chrom}.sv.ids \
                 --vcf  {input.vcf} \
                 --threads {threads} \
                 --r-phased \
                 --ld-snp-list {wildcards.chrom}.sv.ids \
                 --out  {resources.tmpdir}tmp 
        
       mv {resources.tmpdir}tmp.vcor {params.out_prefix}.vcor
       mv {resources.tmpdir}tmp.log {params.out_prefix}.log
       mv {resources.tmpdir}tmp.prune.in {params.out_prefix}.prune.in
       mv {resources.tmpdir}tmp.prune.out {params.out_prefix}.prune.out
       

        rm -rf  {resources.tmpdir}
        """


rule GWAS_spy:
    input:
        gt_vcf= OUTPUT_DIR+ "beagle/{chrom}.tagged.vcf.gz"
    output:
        gt_vcf= OUTPUT_DIR+"GWAS_spy/{chrom}.vcf.gz",
        tbi= OUTPUT_DIR+"GWAS_spy/{chrom}.vcf.gz.tbi"   
    resources:
        mem_mb=50*1024,
        cores=32,
        nodes = 1,
        meduim=1,
        runtime = 60 * 16,
        tmpdir= lambda wildcards: "%sgenotyper_gwaspy_%s/"%(tempFolder,f"{wildcards.chrom}"),
        partition =  "med2",
    threads: 16
    log: OUTPUT_DIR+"GWAS_spy/{chrom}.log"	
    shell:
        '''
	    mkdir -p {resources.tmpdir}
        preimp_qc --dirname {OUTPUT_DIR}beagle/ --basename {wildcards.chrom} --annotations annotation.all.tsv --input-type vcf --report False --reference GRCh38  --out-dir {resources.tmpdir} --export-type vcf  &> {log} || true
	    cp {resources.tmpdir}GWASpy/Preimp_QC/{wildcards.chrom}_qced.vcf.bgz {output.gt_vcf}
        tabix -p vcf {output.gt_vcf}
	    rm  -rf {resources.tmpdir}
	'''

