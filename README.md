
<h1 align="center"> @DIB-LAB/The Great Genotyper </h1>
 <img src="https://cdn.pixabay.com/photo/2017/06/04/16/33/pyramids-2371501_1280.jpg" alt="Logo"/>
<details>
<summary>📖 Table of Contents</summary>
<br />

[![-----------------------------------------------------](https://raw.githubusercontent.com/andreasbm/readme/master/assets/lines/colored.png)](#table-of-contents)

## ➤ Table of Contents

- [➤ Table of Contents](#-table-of-contents)
- [➤ Introduction](#-introduction)
- [➤ Quick Installation](#-quick-installation-pip)
- [➤ Build from source](#-build-from-source)
  - [Clone](#clone)
  - [Install dependencies](#install-dependencies)
  - [Build](#build)
    - [CMake options](#cmake-options)
    - [**Build The kProcessor Library**](#build-the-kprocessor-library)
    - [**Build Everything**](#build-everything)
- [➤ Manually build the Python bindings](#-manually-build-the-python-bindings)
  - [Generate bindings](#generate-bindings)
- [➤ Contributors](#-contributors)
- [➤ License](#-license)

</details>

[![-----------------------------------------------------](https://raw.githubusercontent.com/andreasbm/readme/master/assets/lines/colored.png)](#introduction)

## ➤ Introduction

**The Great Genotyper** is  a mapping-free population genotyping. GG uses metagraph to index the population of samples using its kmer content; indexing runs once per population. For each SV to be genotyped, GG uses pangenie model to genotype SV based on the kmer counts for its marker kmers. Finally, GG will improve the results further using the imputation model based on the genotyping results. In conclusion, GG can accurately estimate AF using the same resources as variant calling pipelines.

[![-----------------------------------------------------](https://raw.githubusercontent.com/andreasbm/readme/master/assets/lines/colored.png)](#quick_installation)

## ➤ Quick Installation

Download the [precompiled binary](https://github.com/dib-lab/TheGreatGenotyper/releases/download/untagged-6f4ca3f2f787ecf0e1e0/TheGreatGenotyper)

[![-----------------------------------------------------](https://raw.githubusercontent.com/andreasbm/readme/master/assets/lines/colored.png)](#build_source)

## ➤ Build from source

### Clone

```bash
git clone https://github.com/dib-lab/TheGreatGenotyper.git
cd TheGreatGenotyper/
```


### Install dependencies

```bash
conda env create -f environment.yml
conda activate gg
conda env config vars set CPATH=${CONDA_PREFIX}/include:${CPATH}
conda activate
conda activate gg
```

### Build

```bash=
# Run CMake configure
cmake -Bbuild

# Run make with parallel execution.
cmake --build build -j4 # -j4 = execute 4 recipes simultaneously.
```


</details>

[![-----------------------------------------------------](https://raw.githubusercontent.com/andreasbm/readme/master/assets/lines/colored.png)](#manual_build_python)

## ➤ Run 
uncompress test_data
```
gzip -d test_data/GRCh38_chr21.fa.gz test_data/test.vcf.gz
```
Download index(50G) for Simons Genome Diversity population([SGDP](https://www.nature.com/articles/nature18964))

```
mkdir -p index/SGDP/
cd index/SGDP/
wget https://farm.cse.ucdavis.edu/~mshokrof/SGDP/graph.dbg
wget https://farm.cse.ucdavis.edu/~mshokrof/SGDP/graph.desc.tsv
wget https://farm.cse.ucdavis.edu/~mshokrof/SGDP/annotation.relaxed.row_diff_int_brwt.annodbg
wget https://farm.cse.ucdavis.edu/~mshokrof/SGDP/samples.csv
cd -
```

Run The great genotyper
```
mkdir test_output
./build/pangenie/src/TheGreatGenotyper -j 32 -t 32 -a index/SGDP/annotation.relaxed.row_diff_int_brwt.annodbg -f index/SGDP/graph.desc.tsv -g -i index/SGDP/graph.dbg -r  test_data/GRCh38_chr21.fa  -v test_data/test.vcf -o test_output/test  &> test_output/log
```

Postprocessing
```
parallel --gnu 'bgzip {} && tabix -p vcf {}.gz' ::: test_output/test_*_genotyping.vcf
scripts/populationMerge.sh test_output/merged.vcf.bgz test_output/test_*_genotyping.vcf.gz
```

Query the frequent variants
```
bcftools view  -q 0.9 test_output/merged.vcf.bgz |grep -vP "^#" |head
```

Query the rare variants
```
bcftools view  -Q 0.1 test_output/merged.vcf.bgz |grep -vP "^#" |head
```





[![-----------------------------------------------------](https://raw.githubusercontent.com/andreasbm/readme/master/assets/lines/colored.png)](#contributors)

## ➤ Contributors

Moustafa Shokrof
C.Titus Brown
Tamer Mansour

[![-----------------------------------------------------](https://raw.githubusercontent.com/andreasbm/readme/master/assets/lines/colored.png)](#license)

## ➤ License

Licensed under [BSD-3-Clause](https://opensource.org/licenses/BSD-3-Clause).
