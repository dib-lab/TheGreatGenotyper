# Database Builder
* This script constructs a Metagraph index and is a Snakemake adaptation of the original [script](https://github.com/ratschlab/metagraph/tree/master/projects/kingsford)([explained here](https://metagraph.ethz.ch/static/docs/quick_start.html#index-k-mer-counts)) created by [Mikhail Karasikov](https://github.com/karasikov)


* To run this workflow, make sure you have Snakemake version 7.0.0 or higher and Conda for installing necessary tools: Metagraph and KMC.

* After building TheGreatGenotyper, Go to the DatabaseBuilder folder in the build directory

  ```
  cd build/DatabaseBuilder
  ``` 


* We are using [PEP](https://pep.databio.org/en/2.1.0/specification/) specification to describe the input for the workflow, and you'll need to edit three files before running the Snakemake workflow:

  1. config.yaml: This file contains various configuration options for the workflow.
  2. sample_table.csv: This CSV file should contain information about your samples.
  3. subsample_table.csv: This CSV file should contain information about your subsamples.

### config.yaml
```
outputFolder:  result/
tempFolder: /scratch/mshokrof/
kSize: 31
kmcMinCount: 3
batchSize: 16
downloadUrls: 
```
1. outputFolder: The location where the generated index will be saved, specifically at result/smooth_1 within this folder.
2. tempFolder: A directory used for writing temporary files during the script's execution. Ideally, this should be a high-speed disk for efficient processing.
3. kSize: This refers to the kmer size that will be used while creating the kmer index.
4. kmcMinCount: The minimum count for kmers. Any kmer with a count below this value will be discarded during the process.
5. batchSize: This is the batch size used by Metagraph when processing samples. Larger batches will process faster but will also require more memory. Note that the batch size doesn't affect the final output.
6. downloadUrls: This feature allows the workflow to download samples from external servers. Detailed documentation on this feature will be provided later.

### sample_table.csv
```
sample_name
sample1
sample2
sample3
```

The "sample_table.csv" file contains a list of sample names that you want to use to build the index. Edit this file and include your sample namesand don't remove the header, which should be "sample_name". 

### subsample_table.csv
```
sample_name,file
sample1,/path/sample1_1.fastq.gz
sample1,/path/sample1_2.fastq.gz
sample2,/path/sample2_1.fastq.gz
sample2,/path/sample2_2.fastq.gz
sample3,/path/sample3.fastq.gz
```

The "subsample_table.csv" file contains the file paths for the FASTQ files associated with your samples. If a sample has multiple FASTQ files, you should add a new line for each file in this .csv file.


## running the workflow

Running the workflow is very simple after finishing the configurations
```
snakemake -j 64  --use-conda
```

## Alignment-free QC and partitioing 

It is always advised to check for the sequencing data quality, extent of genome coverage and sequencing depth. We uses [Snipe](https://github.com/snipe-bio/snipe) to sketch sequencing datasets and perform an alignment-free quality control before indexing. Moreover, for large population cohorts, we partition the data into groups, each has up to 350 samples. Then, each parition of samples is indexed independently. Later, users can genotype one or more indexes at a time. The Great Genotyper can eventually merge the genotype output files from multiple indexes then apply population-level QC and imputation. The larger the size of each index, the more memory is needed for both construction and genotyping with a trade in the time needed to genotype the whole population. We can further reduce the memeory needed for indexing by partitioning the samples based on their sequence similarity. Snipe has a [plugin](https://github.com/snipe-bio/snipe-popclust) that calculates pairwise similarities between all samples based on their sketches and creates a dendrogram that can be used to partition the data into homogenous groups. 


   

