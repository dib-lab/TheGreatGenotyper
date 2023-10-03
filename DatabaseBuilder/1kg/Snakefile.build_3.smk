pepfile: "project_config.yaml"
configfile: "config.yaml" 

tmpFolder=config["tempFolder"]
outputFolder=config["outputFolder"]
kSize=config["kSize"]
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
    partitions=["bml","low2"]
    p=hash(sample)%len(partitions)
    return partitions[p]




allSamples=[s.sample_name for s in pep.samples]
clusters={}
samples_files={}
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
            "smoothed_unitigs": [],
            "batches": [],
            "columns": [],
            "columnsCount": [],
            "rowdiff": [],
            "rowdiffCount": []
            }
    samples_files[s.sample_name] = {}

    clusters[s.cluster]["samples"].append(s.sample_name)
    unitig= list(filter(lambda x: x[-9:] == ".fasta.gz", s.file))[0]
    samples_files[s.sample_name]['unitig'] = unitig 
    clusters[s.cluster]["unitigs"].append(unitig)

    smoothed_unitigs = f"{outputFolder}unitigs/smooth_10000000/{s.sample_name}.fasta.gz"

    clusters[s.cluster]["smoothed_unitigs"].append(smoothed_unitigs)

    
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
    clusters[c]["batches"] = {}
    for i in range(0,numSamples,batchSize):
        start=i
        end=i+batchSize
        end=min(end,numSamples)
        batchID = int(i/batchSize)
        clusters[c]["batches"][batchID]= clusters[c]["smoothed_unitigs"][start:end]
        flag = f"{outputFolder}{c}/columns/{batchID}/done"
        batchFlags.append(flag)




localrules: all, createGraphDescriptor, prepareAnnnotationColumns, prepareRowDiffs
indexFiles=["annotation.relaxed.row_diff_int_brwt.annodbg", "graph.desc.tsv", "graph.dbg"]
out = expand("{outputFolder}{cluster}/{file}",outputFolder=outputFolder, cluster=clusters,file= indexFiles)

rule all:
    input:
        out



rule createGraphDescriptor:
    input:
         coverages = "coverage.tsv",
	     unitigs= lambda wildcards: clusters[f"{wildcards.cluster}"]["smoothed_unitigs"],
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
        parallel --gnu -j1 -k "grep  '{{}}' {input.coverages}| cut -f2" ::: {params.labels} > tmp2
	    echo {input.unitigs}| tr -s ' ' $'\n' > tmp3
        paste -d'\t' tmp1 tmp2 tmp3 > {output.description}
	"""

# rule fasta_to_graph:
#     input: unitig=lambda wildcards: samples_files[f"{wildcards.sample}"]['unitig']
#     output:
#          graph= temp(outputFolder+"new_graphs/{sample}.dbg"),
#          weights= temp(outputFolder+"new_graphs/{sample}.dbg.weights"),	 
#     params:
#          outputFolder= outputFolder+"new_graphs/"
#     conda:
#          "env.yaml"
#     threads: 16
#     priority: 1
#     resources:
#         mem_mb=32 *1024,
#         cores=16,
#         mem_gb=30,
# 	nodes = 1,
#         runtime = lambda wildcards, attempt: 60 * 8 * attempt,
# #        partition = lambda wildcards , attempt: getLowPartition(f"{wildcards.sample}", attempt),
#         partition = "med2",
# 	tmp= lambda wildcards: "%sgraph_%s/"%(tmpFolder,f"{wildcards.sample}")
#     log: outputFolder+"new_graphs/{sample}.log"
#     shell:
#        	"""	    
# 	    mkdir  -p {resources.tmp}
# 	    function cleanup() {{ rm -rf {resources.tmp}; }}
# 	    trap cleanup EXIT
# 	    trap cleanup USR1

	   

# 	    metagraph build  \
#             -k {kSize} \
#             --mode canonical \
#             --count-kmers --count-width 16 \
#             --mem-cap-gb {resources.mem_gb} \
#             --disk-swap  {resources.tmp}\
#             -p {threads} \
#             -o  {resources.tmp}graph \
#             {input.unitig} &> {log}
	    
# 	    mv {resources.tmp}graph.dbg {params.outputFolder}{wildcards.sample}.dbg
# 	    mv {resources.tmp}graph.dbg.weights {params.outputFolder}{wildcards.sample}.dbg.weights
# 	    ls -lsah {output.graph} >> {log}
#         rm -rf {resources.tmp}
#         """



 
# rule clean_and_smooth:
#     input:
#          graph= ancient(outputFolder+"new_graphs/{sample}.dbg"),
#          weights= ancient(outputFolder+"new_graphs/{sample}.dbg.weights")
#     output:
#          unitigs= outputFolder+"unitigs/smooth_10000000/{sample}.fasta.gz",
#          counts= outputFolder+"unitigs/smooth_10000000/{sample}.kmer_counts.gz",
#     params:
#          outputPrefix=outputFolder+"unitigs/smooth_10000000/{sample}"
#     conda:
#          "env.yaml"
#     threads: 16
#     priority: 1
#     resources:
#         mem_mb=20000,
#         cores=16,
#         mem_gb=18,
# 	nodes = 1,
#         runtime = lambda wildcards, attempt: 60 * 4 * attempt,
#         partition = lambda wildcards , attempt: getLowPartition(f"{wildcards.sample}", attempt),
#         tmp= lambda wildcards: "%sclean_smooth%s_%s/"%(tmpFolder,10000000,f"{wildcards.sample}")
#     log: outputFolder+"unitigs/smooth_10000000/{sample}.log"
#     shell:
#        	"""
#             mkdir  -p {resources.tmp}
# 	    function cleanup() {{ rm -rf {resources.tmp}; }}
# 	    trap cleanup EXIT
# 	    trap cleanup USR1

# 	    cp {input.graph} {input.weights} {resources.tmp}

#         metagraph clean  \
#         --to-fasta --primary-kmers \
#         --smoothing-window 1 \
#         --prune-tips 0 \
#         --prune-unitigs 1 \
#         --smoothing-window 10000000 \
#         -p {threads} \
#         -o {resources.tmp}{wildcards.sample} \
#         {resources.tmp}{wildcards.sample}.dbg  &> {log} 

# 	    mv {resources.tmp}{wildcards.sample}.fasta.gz {output.unitigs}
# 	    mv {resources.tmp}{wildcards.sample}.kmer_counts.gz {output.counts}
#  	    rm -rf {resources.tmp}
#         """

rule smooth_unitigs_fasta:
    input: unitigs=lambda wildcards: samples_files[f"{wildcards.sample}"]['unitig']
    output:
         unitigs= outputFolder+"unitigs/smooth_10000000/{sample}.fasta.gz",
         counts= outputFolder+"unitigs/smooth_10000000/{sample}.kmer_counts.gz",	 
    params:
         outputPrefix= outputFolder+"unitigs/smooth_10000000/{sample}"
#    conda:
#         "env.yaml"
    threads: 1
    priority: 1
    resources:
        mem_mb=8 *1024,
        cores=1,
        mem_gb=8,
	    nodes = 1,
        runtime = lambda wildcards, attempt: 60 * 1 * attempt,
#        partition = lambda wildcards , attempt: getLowPartition(f"{wildcards.sample}", attempt),
        partition = "med2",
	    tmp= lambda wildcards: "%sgraph_%s/"%(tmpFolder,f"{wildcards.sample}")
    log: outputFolder+"unitigs/smooth_10000000/{sample}.log"
    shell:
       	"""
	  /mnt/gs21/scratch/mansourt/TheGreatGenotyper/build_noName/metagraph smoothCounts \
          -k {kSize}  \
          -o {params.outputPrefix} \
          {input.unitigs} &> {log}
        """




rule metagraph_buildSmoothMainGraph_1:
    input: lambda wildcards: clusters[f"{wildcards.cluster}"]["smoothed_unitigs"]
    output:
         graph=temp(outputFolder+"{cluster}/graph_canonical.dbg")
    params:
         outputPrefix=outputFolder+"{cluster}/graph_canonical.dbg"
    threads: 16
    resources:
        mem_mb= 64 * 1024,
        cores=16,
        mem_gb= 60,
	    nodes = 1,
        runtime = lambda wildcards, attempt: 60 * 6 * attempt,
        partition = lambda wildcards , attempt: getMeduimPartition(f"{wildcards.cluster}"),
        tmp= lambda wildcards: "%smaingraph_smooth%s/"%(tmpFolder,f"{wildcards.cluster}")
    log: outputFolder+"{cluster}/buildMainGraphlog"
    shell:
       	"""
	#    mamba create -p ./metagraph.$$ -c bioconda -c conda-forge metagraph
        #source ~/miniconda3/etc/profile.d/conda.sh && conda activate ./metagraph.$$ 

        mkdir  -p {tmpFolder}$$/
        /mnt/gs21/scratch/mansourt/TheGreatGenotyper/build_noName/metagraph build  \
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
        #mamba create -p ./metagraph.$$ -c bioconda -c conda-forge metagraph
        #source ~/miniconda3/etc/profile.d/conda.sh && conda activate ./metagraph.$$ 

        /mnt/gs21/scratch/mansourt/TheGreatGenotyper/build_noName/metagraph transform  \
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

	#    mamba create -p ./metagraph.$$ -c bioconda -c conda-forge metagraph
        #source ~/miniconda3/etc/profile.d/conda.sh && conda activate ./metagraph.$$ 

        mkdir  -p {tmpFolder}$$/
        /mnt/gs21/scratch/mansourt/TheGreatGenotyper/build_noName/metagraph build  \
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
    threads: 24
    resources:
        mem_mb= 64 * 1024,
        cores=32,
        mem_gb=60,
	    nodes = 1,
        runtime = lambda wildcards, attempt: 60 * 6 * attempt,
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

	/mnt/gs21/scratch/mansourt/TheGreatGenotyper/build_noName/metagraph annotate  \
                -i {input.graph} \
                --anno-filename \
                --separately \
                --count-kmers --count-width 16 \
                -o {params.outputPrefix} \
		        --threads-each 6 \
                -p 4 \
                {input.unitigs} &> {log}

	touch {output.batchFlag}
	"""



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
    threads: 16
    resources:
        mem_mb= 32 * 1024,
        cores=32,
        mem_gb= 30,
	    nodes = 1,
        runtime = lambda wildcards, attempt: 60 * 24 * attempt,
        partition = lambda wildcards , attempt: getMeduimPartition(f"{wildcards.cluster}"),
        tmp= lambda wildcards: "%soptomize_smooth%s/"%(tmpFolder,f"{wildcards.cluster}")
    log: outputFolder+"{cluster}/optimize_stage0.log"
    shell:
       	"""
    #     mamba create -p ./metagraph.$$ -c bioconda -c conda-forge metagraph
    #     source ~/miniconda3/etc/profile.d/conda.sh && conda activate ./metagraph.$$

	mkdir  -p {tmpFolder}$$/
	# mkdir -p {params.rowDiffPrefix}

	/mnt/gs21/scratch/mansourt/TheGreatGenotyper/build_noName/metagraph transform_anno  \
            --anno-type row_diff --count-kmers \
            --row-diff-stage 0 \
            --mem-cap-gb {resources.mem_gb} \
            --disk-swap {tmpFolder}$$ \
            -i {input.graph} \
            -o {output} \
            -p {threads} \
	    {input.columns} &> {log}

    rm -rf {tmpFolder}$$/
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
    threads: 16
    resources:
        mem_mb= 64 * 1024,
        cores= 16,
        mem_gb=62,
	    nodes = 1,
        runtime = lambda wildcards, attempt: 60 * 24 * attempt,
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


	/mnt/gs21/scratch/mansourt/TheGreatGenotyper/build_noName/metagraph transform_anno  \
            --anno-type row_diff --count-kmers \
            --row-diff-stage 1 \
            --mem-cap-gb {resources.mem_gb} \
            --disk-swap {tmpFolder}$$ \
            -i {input.graph} \
            -o {output} \
            -p {threads} \
	    {input.columns} &>> {log}

    rm -rf {tmpFolder}$$/
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
    threads: 16
    resources:
        mem_mb=64 * 1024,
        cores=16,
        mem_gb=62,
	    nodes = 1,
        runtime = lambda wildcards, attempt: 60 * 24 * attempt,
        partition = lambda wildcards , attempt: getMeduimPartition(f"{wildcards.cluster}"),
        tmp= lambda wildcards: "%soptomize_smooth%s/"%(tmpFolder,f"{wildcards.cluster}")
    log: outputFolder+"{cluster}/optimize_stage2.log"
    shell:
       	"""
	 mkdir  -p {tmpFolder}$$/

    # mkdir  -p {tmpFolder}$$/
    #     currDir=$(pwd)
    #     cd /mnt/gs21/scratch/mansourt/TheGreatGenotyper/
    #     cmake -Bbuild.$$
    #     cmake --build build.$$ -j{threads}
    #     cd ${{currDir}}

   
    /mnt/gs21/scratch/mansourt/TheGreatGenotyper/build_noName/metagraph transform_anno  \
            --anno-type row_diff --count-kmers \
            --row-diff-stage 2 \
            --mem-cap-gb {resources.mem_gb} \
            --disk-swap {tmpFolder}$$ \
            -i {input.graph} \
            -o {params.rowDiffPrefix}out \
            -p {threads} \
	    {input.columns} &>> {log}
	echo "stage 2 finished" >> {log}

    rm -rf {tmpFolder}$$/
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
#    conda:
#        "env.yaml"
    threads: 16
    resources:
        mem_mb=256 * 1024,
        cores=32,
        mem_gb=62,
	nodes = 1,
        runtime = lambda wildcards, attempt: 60 * 8 * attempt,
        partition = lambda wildcards , attempt: getMeduimPartition(f"{wildcards.cluster}"),
        tmp= lambda wildcards: "%soptomize_smooth%s/"%(tmpFolder,f"{wildcards.cluster}")
    log: outputFolder+"{cluster}/optimize_stage3.log"
    shell:
       	"""
	
    # mkdir  -p {tmpFolder}$$/
    #     currDir=$(pwd)
    #     cd /mnt/gs21/scratch/mansourt/TheGreatGenotyper/
    #     cmake -Bbuild.$$
    #     cmake --build build.$$ -j{threads}
    #     cd ${{currDir}}


        echo {input.columns} |tr -s ' ' $'\n'  \
         | /mnt/gs21/scratch/mansourt/TheGreatGenotyper/build_noName/metagraph transform_anno  \
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
#    conda:
#        "env.yaml"
    threads: 16
    resources:
        mem_mb=64 * 1024,
        cores=16,
        mem_gb=62,
	    nodes = 1,
        runtime = lambda wildcards, attempt: 60 * 4 * attempt,
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


        /mnt/gs21/scratch/mansourt/TheGreatGenotyper/build_noName/metagraph relax_brwt  \
            -p {threads} \
            --relax-arity 32 \
            -o {params.outputPrefix}.relaxed \
           {input.annotation}  &>> {log}
	    echo "relax finished" >> {log}
	    
	    rm {input.graph}.*
	"""

