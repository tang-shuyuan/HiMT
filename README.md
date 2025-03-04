# HiMT 
![GitHub](https://github.com/user-attachments/assets/9c086c04-09b7-4958-9abe-804b0d5fd27e)


An Integrative Toolkit for Assembling Plant Mitochondrial Genomes Using HiFi Reads 

HiMT is a rapid plastid genome assembly tool for HiFi data, which can be used for the assembly of plant mitochondrial, chloroplast, and animal mitochondrial genomes. It supports usage on Linux, Windows, and macOS, and can be used through both command line and GUI.

Find source codes and documentation at (https://github.com/tang-shuyuan/HiMT)  
Find detailed documentation at (https://www.yuque.com/yuqueyonghuwrgkbo/tonqgq?#)  
For any question about SynGAP, please contact 2421496996@qq.com
## Installation under Linux
conda (recommended)
```
conda create -n himt
conda activate himt
conda install shuyuan_tang::himt -c bioconda
```
## Installation under Windows

click the [link](https://github.com/tang-shuyuan/HiMT/releases/download/untagged-591749e6e2f1267a3b1a/himt_windows.tar.gz)


## Installation under MacOS
conda (recommended)
```
conda create -n himt
conda activate himt
conda install shuyuan_tang::himt -c bioconda
```
## Dependence
```
    - python >=3.10
    - blast >=2.14.0
    - flye >=2.9.4
    - miniprot >=0.13
    - plotly >=5.24.0
    - numpy >=1.26.0

```
# Quickly start
#### Introduction
HiMT includes four subprograms
```
usage: himt [function] [argument]

a toolkit for assembling mitochondrial genome
version 1.0.7
for any questions, please put them forward on https://github.com/tang-shuyuan/himt or https://bioanno.com.

options:
  -h, --help     show this help message and exit
  -v, --version  show program's version number and exit

function:

      assemble      Assemble mitochondrial genome with HiFi data
      assess        Assess the assembly quality of the mitochondrial genome
      filter        Filter low-depth nuclear genome sequencing reads
      compare       Visualize the collinearity between two genomes
```

### Assemble
Run the command below to assemble mitochondrial and chloroplast genomes simultaneously.
```
usage: himt assemble [argument]
please use 'himt assemble -h or --help' to show help information

Assembling mitochondrail genome with HiFi sequcncing data

Required arguments:
  -i INPUT_FILE, --input_file INPUT_FILE
                        input a fasta or fastq file,gz compressed files are supported.
  -o OUTPUT_DIR, --output_dir OUTPUT_DIR
                        output directory.

Optional arguments:
  -h, --help            Show this help message and exit
  -n HEAD_NUMBER, --head_number HEAD_NUMBER
                        default=4,The number of kmer species randomly selected.
  -t THREAD, --thread THREAD
                        default=2 The number of thread used during code execution.
  -b {3,4}, --base_number {3,4}
                        default=3,only can be 3 and 4,The number of bases at the beginning of kmer.
  -fd FILTER_DEPTH, --filter_depth FILTER_DEPTH
                        read depths below this value will be filtered. You must input the -p parameter to enable the
                        use of the -fd parameter
  -fp FILTER_PERCENTAGE, --filter_percentage FILTER_PERCENTAGE
                        default=0.3,The depth of the mitochondrial genome obtained by blast, the proportion adjusted
                        downwards on this value.
  -p PROPORTION, --proportion PROPORTION
                        default=1,The percentage of the selected dataset from the entire file, choose a value from
                        0-1.
  -c ACCURACY, --accuracy ACCURACY
                        default=0.8,If one read has a high-frequency kmer ratio exceeding this value, it will be
                        considered as a high-frequency read,choose a value from 0-1.
  -s {plant,animal}, --species {plant,animal}
                        default=plant,Species category,only can be plant or animal.
  --no_flye_meta        By default, we use flye Meta pattern to assemble the mitochondrial genome. If you don't want
                        to use meta pattern, add this parameter.
  -x NORMALIZE_DEPTH, --normalize_depth NORMALIZE_DEPTH
                        Normalize the mitochondrial genome depth to a value.If the input value exceeds the
                        mitochondrial genome depth, retain the maximum mitochondrial genome depth,the default
                        mitochondrial genome depth ranges between 15 and 50. input a value less than 0 (such as:-1) to
                        retain the maximum mitogenome depth
```

#### Examples
Demo data download

click the [link](https://github.com/tang-shuyuan/HiMT/releases/download/untagged-591749e6e2f1267a3b1a/demo.fa) 
or
```
wget https://github.com/tang-shuyuan/HiMT/releases/download/untagged-591749e6e2f1267a3b1a/demo.fa
```


1. run himt to assemble plant organelle genomes(HiMT support fa/fq/fa.gz/fq/gz)

```
himt assemble -i pineapple.LY.fa -o output -t 10
```

2. If you occasionally notice incomplete mitochondrial genome assemblies with HiMT, you can try the following commands:

```
himt assemble -i input.file -o output -fp 0.2 -x 50
```

3. To assemble animal mitogenome

```
himt assemble -i input_file -o output_dir -s animal
```

#### Main output files
| file | description |
| --- | --- |
| extract.fa |The filtered high-depth data can be directly used for mitochondrial genome assembly. |
| himt_mitochondrial.gfa | Mitochondrial genome result file |
| himt_mitochondrial.html | Evaluation of mitochondrial genome assembly results |
| himt_chloroplast.gfa | Chloroplast genome result file |
| himt_chloroplast.html | Evaluation of chloroplast genome assembly results |
| chloroplast_hap1.fa | Chloroplast haplotype 1 |
| chloroplast_hap2.fa | Chloroplast haplotype 2 |


### Filter 
only filter data but not assemble organelle genome

```
usage: himt filter [argument]
please use 'himt filter -h or --help' to show help information

Filtering the low-depth nuclear genome reads

Required arguments:
  -i INPUT_FILE, --input_file INPUT_FILE
                        input a fasta or fastq file,gz compressed files are supported.
  -o OUTPUT_DIR, --output_dir OUTPUT_DIR
                        output directory.

Optional arguments:
  -h, --help            Show this help message and exit
  -n HEAD_NUMBER, --head_number HEAD_NUMBER
                        default=4,The number of kmer species randomly selected.
  -t THREAD, --thread THREAD
                        default=2 The number of thread used during code execution.
  -b {3,4}, --base_number {3,4}
                        default=3,only can be 3 and 4,The number of bases at the beginning of kmer.
  -fd FILTER_DEPTH, --filter_depth FILTER_DEPTH
                        read depths below this value will be filtered. You must input the -p parameter to enable the
                        use of the -fd parameter
  -fp FILTER_PERCENTAGE, --filter_percentage FILTER_PERCENTAGE
                        default=0.3,The depth of the mitochondrial genome obtained by blast, the proportion adjusted
                        downwards on this value.
  -p PROPORTION, --proportion PROPORTION
                        The percentage of the selected dataset from the entire file, choose a value from 0-1.
  -c ACCURACY, --accuracy ACCURACY
                        default=0.8,If one read has a high-frequency kmer ratio exceeding this value, it will be
                        considered as a high-frequency read,choose a value from 0-1.
  -s {plant,animal}, --species {plant,animal}
                        default=plant,Species category,only can be plant or animal.
  -x NORMALIZE_DEPTH, --normalize_depth NORMALIZE_DEPTH
                        Normalize the mitochondrial genome depth to a value.If the input value exceeds the
                        mitochondrial genome depth, retain the maximum mitochondrial genome depth,the default
                        mitochondrial genome depth ranges between 15 and 50. input a value less than 0 (such as:-1) to
                        retain the maximum mitogenome depth
```

#### Examples
```
himt filter -i input_file -o output -t 10
```

### Assess
If you have an assembled plant mitochondrial genome and wish to evaluate the assembly quality.

```
usage: himt assess [argument]
please use 'himt assess -h or --help' to show help information

Assessing the assembly quality of the mitochondrial genome

options:
  -h, --help            show this help message and exit
  -c {mitochondrial,chloroplast}, --category {mitochondrial,chloroplast}
                        default=mitochondrial,Choose the category of organelles you want to assess.

Required arguments:
  -i INPUT_FILE, --input_file INPUT_FILE
                        input a fasta or gfa file.
  -o OUTPUT_DIR, --output_dir OUTPUT_DIR
                        output directory.
```

#### Examples
Support fasta format and gfa format

```
himt assess -i input.fa/gfa -o assess_output
```

### Compare 


```
usage: himt compare [argument]
please use 'himt compare -h or --help' to show help information

compare the collinearity between two genomes based on alignment results from minimap2

Required arguments:
  -r REFERENCE, --reference REFERENCE
                        input reference genome
  -q QUERY, --query QUERY
                        input query genome
  -o OUTPUT_DIR, --output_dir OUTPUT_DIR
                        output directory

Optional arguments:
  -c {mitochondrial,chloroplast,other}, --category {mitochondrial,chloroplast,other}
                        default=mitochondrial,choose the category of genome you want to compare
  -h, --help            Show this help message and exit
```

#### Examples

Support fasta format and gfa format
```
himt compare -q genome1.fa/gfa -r genome2.fa/gfa -o output
```

