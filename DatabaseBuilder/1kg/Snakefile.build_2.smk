pepfile: "project_config.yaml"
configfile: "config.yaml" 

tmpFolder=config["tempFolder"]
outputFolder=config["outputFolder"]
kSize=config["kSize"]
urlsFile=config["downloadUrls"]
batchSize=config["batchSize"]

kmcMinCount=3
if "kmcMinCount" in config:
   kmcMinCount= config["kmcMinCount"]



globus_ebi_id= ""
globus_client_id= ""

if "globus_client_id" in config:
   globus_ebi_id= config["globus_ebi_id"]
   globus_client_id= config["globus_client_id"]
   

import random

def getHighPartition(sample):
    partitions=["high2","bmh"]
    p=hash(sample)%len(partitions)
    return partitions[p]

def getMeduimPartition(sample):
    partitions=["med2"]
    p=hash(sample)%len(partitions)
    return partitions[p]

def getLowPartition(sample,attempt):
    if attempt>1:
       return getMeduimPartition(sample)
    partitions=["bml"]
    p=hash(sample)%len(partitions)
    return partitions[p]




allSamples=[s.sample_name for s in pep.samples]
clusters={}
#list(set([s.cluster for s in pep.samples]))

allColumns= []
allColumnsCounts= []

allRowDiff= []
allRowDiffCounts= []


for s in pep.samples:
    if s.cluster not in clusters:
        clusters[s.cluster] = {
            "samples": [],
            "unitigs":[] ,
            "histograms": [], 
            "batches": [],
            "columns": [],
            "columnsCount": [],
            "rowdiff": [],
            "rowdiffCount": []
            }
    
    clusters[s.cluster]["samples"].append(s.sample_name)
    unitig= list(filter(lambda x: x[-9:] == ".fasta.gz", s.file))[0]
    clusters[s.cluster]["unitigs"].append(unitig)
    
    histo= list(filter(lambda x: x[-6:] == ".histo", s.file))[0]
    clusters[s.cluster]["histograms"].append(unitig)
    column=      f"{outputFolder}{s.cluster}/columns/{s.sample_name}.fasta.gz.column.annodbg"
    columnCount= f"{outputFolder}{s.cluster}/columns/{s.sample_name}.fasta.gz.column.annodbg.counts"
    rowdiff=     f"{outputFolder}{s.cluster}/rowDiff/{s.sample_name}.fasta.gz.column.annodbg"
    rowdiffCount=     f"{outputFolder}{s.cluster}/rowDiff/{s.sample_name}.fasta.gz.column.annodbg.counts"
        
    clusters[s.cluster]["columns"].append(column)
    clusters[s.cluster]["columnsCount"].append(columnCount)
    allColumns.append(column)
    allColumns.append(columnCount)



    clusters[s.cluster]["rowdiff"].append(rowdiff)
    clusters[s.cluster]["rowdiffCount"].append(rowdiffCount)

    allRowDiff.append(rowdiff)
    allRowDiffCounts.append(rowdiffCount)



batchFlags=[]

for c in clusters.keys():
    numSamples = len(clusters[c]["samples"])
    for i in range(0,numSamples,batchSize):
        start=i
        end=i+batchSize
        end=min(end,numSamples)
        clusters[c]["batches"]= clusters[c]["unitigs"][start:end]
        flag = f"{outputFolder}{c}/columns/{i}/done"
        batchFlags.append(flag)


print(clusters)

localrules: all, createGraphDescriptor
indexFiles=["annotation.relaxed.row_diff_int_brwt.annodbg", "graph.desc.tsv", "graph.dbg"]
out = expand("{outputFolder}{cluster}/{file}",outputFolder=outputFolder, cluster=clusters,file= indexFiles)

rule all:
    input:
        out



rule createGraphDescriptor:
    input:
         histograms = lambda wildcards: clusters[f"{wildcards.cluster}"]["histograms"],
	     unitigs= lambda wildcards: clusters[f"{wildcards.cluster}"]["unitigs"],
         graph=outputFolder+"{cluster}/graph.dbg"
    params:
         labels =  lambda wildcards: clusters[f"{wildcards.cluster}"]["samples"]
    output:
         description=outputFolder+"{cluster}/graph.desc.tsv"
    conda:
         "env.yaml"
    threads: 1
    resources:
        mem_mb=3000,cores=1,mem_gb=2
    log: outputFolder+"{cluster}/graph.desc.tsv.log"
    shell:
       	"""
        echo {params.labels} | tr -s ' ' $'\n' > tmp1
        parallel --gnu -j1 -k "grep  'parameters' {{}} | cut -f3" ::: {input.histograms} > tmp2
	echo {input.unitigs}| tr -s ' ' $'\n' > tmp3
        paste -d'\t' tmp1 tmp2 tmp3 > {output.description}
	"""


rule metagraph_buildSmoothMainGraph_1:
    input: lambda wildcards: clusters[f"{wildcards.cluster}"]["unitigs"]
    output:
         graph=temp(outputFolder+"{cluster}/graph_canonical.dbg")
    params:
         outputPrefix=outputFolder+"{cluster}/graph_canonical.dbg"
    threads: 16
    resources:
        mem_mb= 64 * 1024,
        cores=32,
        mem_gb= 60,
	nodes = 1,
        runtime = lambda wildcards, attempt: 60 * 12 * attempt,
        partition = lambda wildcards , attempt: getMeduimPartition(f"{wildcards.cluster}"),
        tmp= lambda wildcards: "%smaingraph_smooth%s/"%(tmpFolder,f"{wildcards.cluster}")
    log: outputFolder+"{cluster}/buildMainGraphlog"
    shell:
       	"""
	mamba create -p ./metagraph.$$ -c bioconda -c conda-forge metagraph
        source ~/miniconda3/etc/profile.d/conda.sh && conda activate ./metagraph.$$
	#conda activate ./metagraph.$$ 

        mkdir  -p {tmpFolder}$$/
        metagraph build  \
        -k {kSize} \
        --mode canonical  \
        --mem-cap-gb {resources.mem_gb} \
        --disk-swap {tmpFolder}$$/  \
        -p {threads} \
        -o {params.outputPrefix} \
        {input} &> {log}
	    rm -rf {tmpFolder}$$/
        """

rule metagraph_buildSmoothMainGraph_2:
    input:
        graph=outputFolder+"{cluster}/graph_canonical.dbg"
    output:
         graph_fasta=temp(outputFolder+"{cluster}/graph_primary.fasta.gz")
    params:
         outputPrefix=outputFolder+"{cluster}/graph_primary"
#    conda:
#         "../../environment.yml"
    threads: 8
    resources:
        mem_mb= 24 * 1024,
        cores=32,
        mem_gb=22,
	nodes = 1,
        runtime = lambda wildcards, attempt: 60 * 2 * attempt,
        partition = lambda wildcards , attempt: getMeduimPartition(f"{wildcards.cluster}"),
        tmp= lambda wildcards: "%smaingraph_smooth%s/"%(tmpFolder,f"{wildcards.cluster}")
    log: outputFolder+"{cluster}/buildMainGraphlog"
    shell:
       	"""
	mamba create -p ./metagraph.$$ -c bioconda -c conda-forge metagraph
        source ~/miniconda3/etc/profile.d/conda.sh && conda activate ./metagraph.$$
	#conda activate ./metagraph.$$
	


        metagraph transform  \
            --to-fasta --primary-kmers \
	        -p {threads} \
            -o  {params.outputPrefix} \
            {input.graph} &>> {log}
	
        """

rule metagraph_buildSmoothMainGraph_3:
    input:
        graph_fasta=outputFolder+"{cluster}/graph_primary.fasta.gz"
    output:
        graph=outputFolder+"{cluster}/graph.dbg"
    params:
        outputPrefix=outputFolder+"{cluster}/graph.dbg"
#    conda:
#        "../../environment.yml"
    threads: 8
    resources:
        mem_mb= 48 * 1024,
        cores=32,
        mem_gb=45,
	    nodes = 1,
        runtime = lambda wildcards, attempt: 60 * 4 * attempt,
        partition = lambda wildcards , attempt: getMeduimPartition(f"{wildcards.cluster}"),
        tmp= lambda wildcards: "%smaingraph_smooth%s/"%(tmpFolder,f"{wildcards.cluster}")
    log: outputFolder+"{cluster}/buildMainGraphlog"
    shell:
       	"""
        mamba create -p ./metagraph.$$ -c bioconda -c conda-forge metagraph
        source ~/miniconda3/etc/profile.d/conda.sh && conda activate ./metagraph.$$
        #conda activate ./metagraph.$$

        mkdir  -p {tmpFolder}$$/
        metagraph build  \
            -k {kSize} \
            --mode primary \
            --mem-cap-gb {resources.mem_gb} \
            --disk-swap {tmpFolder}$$/ \
            -p {threads} \
            -o {params.outputPrefix} \
            {input.graph_fasta}  &>> {log}
	    
	    rm -rf {tmpFolder}$$/
        """




rule createAnnotationColumns:
    input:
         unitigs=lambda wildcards: clusters[f"{wildcards.cluster}"]["batches"][int(f"{wildcards.batchID}")],
         graph=outputFolder+"{cluster}/graph.dbg"
    output:
         batchFlag = outputFolder+"{cluster}/columns/{batchID}/done",
    params:
         outputPrefix=outputFolder+"{cluster}/columns/{batchID}/"
    wildcard_constraints:
        cluster= "[^/]*"
#    conda:
#         "../../environment.yml"
    threads: 8
    resources:
        mem_mb= 64 * 1024,
        cores=32,
        mem_gb=60,
	nodes = 1,
        runtime = lambda wildcards, attempt: 60 * 4 * attempt,
        partition = lambda wildcards , attempt: getMeduimPartition(f"{wildcards.batchID}"),
        tmp= lambda wildcards: "%screateAnnotation_smooth%s_batch%s/"%(tmpFolder,f"{wildcards.cluster}",f"{wildcards.batchID}")
    log: outputFolder+"{cluster}/columns/{batchID}.done.log"
    shell:
       	"""
	# mamba create   -p ./metagraph.$$ 
    #     source ~/miniconda3/etc/profile.d/conda.sh && conda activate ./metagraph.$$
	# mamba install --file build_tools -c bioconda  -c conda-forge
	# currDir=$(pwd)
    #     cd /mnt/gs21/scratch/mansourt/TheGreatGenotyper/
    #     cmake -Bbuild.$$
    #     cmake --build build.$$ -j{threads}
    #     cd ${{currDir}}

	/mnt/gs21/scratch/mansourt/TheGreatGenotyper/build/metagraph annotate  \
                -i {input.graph} \
                --anno-filename \
                --separately \
                --count-kmers --count-width 16 \
                -o {params.outputPrefix} \
		--threads-each 4 \
                -p 2 \
                {input.unitigs} &> {log}

	touch {output.batchFlag}
	"""


print(batchFlags)
rule prepareAnnnotationColumns:
    input:
         flags=batchFlags
    output:
         columns= allColumns,
         columnsCounts= allColumnsCounts
    conda:
         "env.yaml"
    threads: 1
    resources:
        mem_mb=2048,cores=1,mem_gb=2
    shell:
        """
        ls {input.flags} |sed -e 's/done//' > input.$$
        ls {input.flags} |sed -e 's/[0-9]*\/done//' > output.$$
	    paste input.$$ output.$$ | parallel --col-sep $"\t" --gnu -j1 "mv {{1}}/* {{2}}"
	
	    """





rule optimizeAnnotationColumns_1:
    input:
         columns=      lambda wildcards: clusters[f"{wildcards.cluster}"]["columns"],
         columnsCounts=lambda wildcards: clusters[f"{wildcards.cluster}"]["columnsCount"],
         graph=outputFolder+"{cluster}/graph.dbg"
    output:
         temp(outputFolder+"{cluster}/rowDiff/graph.row_count")
    wildcard_constraints:
        cluster= "[^/]*"
    params:
         rowDiffPrefix = outputFolder + "{cluster}/rowDiff/",
	 outputPrefix  = outputFolder + "{cluster}/annotation",
	 outputFolder  = outputFolder + "{cluster}/"
#    conda:
#         "../../environment.yml"
    threads: 8
    resources:
        mem_mb=64 * 1024,
        cores=32,
        mem_gb=60,
	nodes = 1,
        runtime = lambda wildcards, attempt: 60 * 12 * attempt,
        partition = lambda wildcards , attempt: getMeduimPartition(f"{wildcards.cluster}"),
        tmp= lambda wildcards: "%soptomize_smooth%s/"%(tmpFolder,f"{wildcards.cluster}")
    log: outputFolder+"{cluster}/optimize_stage0.log"
    shell:
       	"""
    #     mamba create -p ./metagraph.$$ -c bioconda -c conda-forge metagraph
    #     source ~/miniconda3/etc/profile.d/conda.sh && conda activate ./metagraph.$$

	# mkdir  -p {tmpFolder}$$/
	# mkdir -p {params.rowDiffPrefix}

	/mnt/gs21/scratch/mansourt/TheGreatGenotyper/build/metagraph transform_anno  \
            --anno-type row_diff --count-kmers \
            --row-diff-stage 0 \
            --mem-cap-gb {resources.mem_gb} \
            --disk-swap {tmpFolder}$$ \
            -i {input.graph} \
            -o {output} \
            -p {threads} \
	    {input.columns} &> {log}
 	echo "stage 0 finished" >> {log}
	"""


rule optimizeAnnotationColumns_2:
    input:
         columns=      lambda wildcards: clusters[f"{wildcards.cluster}"]["columns"],
         columnsCounts=lambda wildcards: clusters[f"{wildcards.cluster}"]["columnsCount"],
         graph=outputFolder+"{cluster}/graph.dbg",
         row_count=outputFolder+"{cluster}/rowDiff/graph.row_count"
    output:
         temp(outputFolder+"{cluster}/rowDiff/graph.row_reduction")
    params:
         rowDiffPrefix = outputFolder + "{cluster}/rowDiff/",
	 outputPrefix  = outputFolder + "{cluster}/annotation",
	 outputFolder  = outputFolder + "{cluster}/"
    wildcard_constraints:
        cluster= "[^/]*"
    conda:
         "../../environment.yml"
    threads: 32
    resources:
        mem_mb=110000,
        cores=32,
        mem_gb=100,
	nodes = 1,
        runtime = lambda wildcards, attempt: 60 * 48 * attempt,
        partition = lambda wildcards , attempt: getMeduimPartition(f"{wildcards.cluster}"),
        tmp= lambda wildcards: "%soptomize_smooth%s/"%(tmpFolder,f"{wildcards.cluster}")
    log: outputFolder+"{cluster}/optimize_stage1.log"
    shell:
       	"""
	mkdir  -p {tmpFolder}$$/

    # mkdir  -p {tmpFolder}$$/
    #     currDir=$(pwd)
    #     cd /mnt/gs21/scratch/mansourt/TheGreatGenotyper/
    #     cmake -Bbuild.$$
    #     cmake --build build.$$ -j{threads}
    #     cd ${{currDir}}


	/mnt/gs21/scratch/mansourt/TheGreatGenotyper/build/metagraph transform_anno  \
            --anno-type row_diff --count-kmers \
            --row-diff-stage 1 \
            --mem-cap-gb {resources.mem_gb} \
            --disk-swap {tmpFolder}$$ \
            -i {input.graph} \
            -o {output} \
            -p {threads} \
	    {input.columns} &>> {log}
 	echo "stage 1 finished" >> {log}
	"""

rule optimizeAnnotationColumns_3:
    input:
         columns=      lambda wildcards: clusters[f"{wildcards.cluster}"]["columns"],
         columnsCounts=lambda wildcards: clusters[f"{wildcards.cluster}"]["columnsCount"],
         graph=outputFolder+"{cluster}/graph.dbg",
         row_reduction= outputFolder+"{cluster}/rowDiff/graph.row_reduction"
    output:
         outputFolder + "{cluster}/rowDiff_tmp/done" 
    params:
         rowDiffPrefix = outputFolder + "{cluster}/rowDiff_tmp/"
    wildcard_constraints:
        cluster= "[^/]*"
    conda:
         "../../environment.yml"
    threads: 32
    resources:
        mem_mb=110000,
        cores=32,
        mem_gb=100,
	nodes = 1,
        runtime = lambda wildcards, attempt: 60 * 48 * attempt,
        partition = lambda wildcards , attempt: getMeduimPartition(f"{wildcards.cluster}"),
        tmp= lambda wildcards: "%soptomize_smooth%s/"%(tmpFolder,f"{wildcards.cluster}")
    log: outputFolder+"{cluster}/optimize_stage2.log"
    shell:
       	"""
	# mkdir  -p {tmpFolder}$$/

    # mkdir  -p {tmpFolder}$$/
    #     currDir=$(pwd)
    #     cd /mnt/gs21/scratch/mansourt/TheGreatGenotyper/
    #     cmake -Bbuild.$$
    #     cmake --build build.$$ -j{threads}
    #     cd ${{currDir}}

   
    /mnt/gs21/scratch/mansourt/TheGreatGenotyper/build/metagraph transform_anno  \
            --anno-type row_diff --count-kmers \
            --row-diff-stage 2 \
            --mem-cap-gb {resources.mem_gb} \
            --disk-swap {tmpFolder}$$ \
            -i {input.graph} \
            -o {params.rowDiffPrefix}out \
            -p {threads} \
	    {input.columns} &>> {log}
	echo "stage 2 finished" >> {log}
    touch {output}
	"""



rule prepareRowDiffs:
    input:
         flags=expand("{outputFolder}{cluster}/rowDiff_tmp/done",outputFolder=outputFolder,cluster=clusters.keys())
    output:
         columns=allRowDiff,
         columnsCounts=allRowDiffCounts
    conda:
         "env.yaml"
    threads: 1
    resources:
        mem_mb=2048,cores=1,mem_gb=2
    shell:
       	"""
        ls {input.flags} |sed -e 's/rowDiff_tmp\/done/rowDiff_tmp\//' > input.$$
        ls {input.flags} |sed -e 's/rowDiff_tmp\/done/rowDiff\//' > output.$$
	    paste input.$$ output.$$ | parallel --col-sep $"\t" --gnu -j1 "mv {{1}}/* {{2}}"
	
    	"""




rule optimizeAnnotationColumns_4:
    input:
         graph=outputFolder+"{cluster}/graph.dbg",
         columns=      lambda wildcards: clusters[f"{wildcards.cluster}"]["rowdiff"],
         columnsCounts=lambda wildcards: clusters[f"{wildcards.cluster}"]["rowdiffCount"]
    output:
         temp(outputFolder+"{cluster}/annotation.row_diff_int_brwt.annodbg")
    params:
         rowDiffPrefix = outputFolder + "{cluster}/rowDiff/",
	 outputPrefix  = outputFolder + "{cluster}/annotation",
	 outputFolder  = outputFolder + "{cluster}/"
    wildcard_constraints:
        cluster= "[^/]*"
    conda:
         "../../environment.yml"
    threads: 32
    resources:
        mem_mb=110000,
        cores=32,
        mem_gb=100,
	nodes = 1,
        runtime = lambda wildcards, attempt: 60 * 48 * attempt,
        partition = lambda wildcards , attempt: getMeduimPartition(f"{wildcards.cluster}"),
        tmp= lambda wildcards: "%soptomize_smooth%s/"%(tmpFolder,f"{wildcards.cluster}")
    log: outputFolder+"{cluster}/optimize_stage3.log"
    shell:
       	"""
	# mkdir  -p {tmpFolder}$$/

    # mkdir  -p {tmpFolder}$$/
    #     currDir=$(pwd)
    #     cd /mnt/gs21/scratch/mansourt/TheGreatGenotyper/
    #     cmake -Bbuild.$$
    #     cmake --build build.$$ -j{threads}
    #     cd ${{currDir}}


        echo {input.columns} |tr -s ' ' $'\n'  \
         | /mnt/gs21/scratch/mansourt/TheGreatGenotyper/build/metagraph transform_anno  \
            --anno-type row_diff_int_brwt \
            --greedy --fast --subsample 1000000  \
            -i {input.graph} \
            -o {params.outputPrefix} \
            -p {threads} --parallel-nodes 10  &>> {log}
	echo "transform to brwt  finished" >> {log}
	"""

rule optimizeAnnotationColumns_relax:
    input:
        graph = outputFolder+"{cluster}/graph.dbg",
        annotation = outputFolder+"{cluster}/annotation.row_diff_int_brwt.annodbg"
    output:
        outputFolder+"{cluster}/annotation.relaxed.row_diff_int_brwt.annodbg"
    params:
        outputPrefix  = outputFolder + "{cluster}/annotation"
    wildcard_constraints:
        cluster= "[^/]*"
    conda:
        "../../environment.yml"
    threads: 32
    resources:
        mem_mb=110000,
        cores=32,
        mem_gb=100,
	    nodes = 1,
        runtime = lambda wildcards, attempt: 60 * 48 * attempt,
        partition = lambda wildcards , attempt: getMeduimPartition(f"{wildcards.cluster}"),
        tmp= lambda wildcards: "%soptomize_smooth%s/"%(tmpFolder,f"{wildcards.cluster}")
    log: outputFolder+"{cluster}/optimize_relax.log"
    shell:
       	"""

        # mkdir  -p {tmpFolder}$$/
        # currDir=$(pwd)
        # cd /mnt/gs21/scratch/mansourt/TheGreatGenotyper/
        # cmake -Bbuild.$$
        # cmake --build build.$$ -j{threads}
        # cd ${{currDir}}


        /mnt/gs21/scratch/mansourt/TheGreatGenotyper/build/metagraph relax_brwt  \
            -p {threads} \
            --relax-arity 32 \
            -o {params.outputPrefix}.relaxed \
           {input.annotation}  &>> {log}
	    echo "relax finished" >> {log}
	    
	    rm {input.graph}.*
	"""

