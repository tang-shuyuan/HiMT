# HiMT![himt_logo](https://github.com/user-attachments/assets/240e7aa3-5ea4-45f4-af66-f29ab79ce732)

HiMT is a rapid plastid genome assembly tool for HiFi data, which can be used for the assembly of plant mitochondrial, chloroplast, and animal mitochondrial genomes. It supports usage on Linux, Windows, and macOS, and can be used through both command line and GUI. 
# installation under Linux
## conda(recommended)

```
conda install shuyuan_tang::himt -c bioconda
```
# Installation under Windows

click the link [HiMT download](https://github.com/tang-shuyuan/HiMT/releases/download/untagged-591749e6e2f1267a3b1a/himt_windows.tar.gz)

# quick start

## Assemble
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

### Example
Download test data
click the link [demo download](https://github.com/tang-shuyuan/HiMT/releases/download/untagged-591749e6e2f1267a3b1a/demo.fa) 

run himt to assemble organelle genomes
```
himt assemble -i demo.fa -o output -t 10
```

To assemble animal mitogenome
```
himt assemble -i input.fa -o output -t 10 -s animal
```

## Filter 
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
                        read depths below this value will be filtered. You must input the -p parameter to enable the use of the -fd
                        parameter
  -fp FILTER_PERCENTAGE, --filter_percentage FILTER_PERCENTAGE
                        default=0.3,The depth of the mitochondrial genome obtained by blast, the proportion adjusted downwards on
                        this value.
  -p PROPORTION, --proportion PROPORTION
                        The percentage of the selected dataset from the entire file, choose a value from 0-1.
  -c ACCURACY, --accuracy ACCURACY
                        default=0.8,If one read has a high-frequency kmer ratio exceeding this value, it will be considered as a
                        high-frequency read,choose a value from 0-1.
  -s {plant,animal}, --species {plant,animal}
                        default=plant,Species category,only can be plant or animal.
  -x NORMALIZE_DEPTH, --normalize_depth NORMALIZE_DEPTH
                        Normalize the mitochondrial genome depth to a value.If the input value exceeds the mitochondrial genome
                        depth, retain the maximum mitochondrial genome depth,the default mitochondrial genome depth ranges between 15
                        and 50. input a value less than 0 (such as:-1) to retain the maximum mitogenome depth
```
### examples
```
himt assemble -i demo.fa -o output -t 10
```

## Assess
If you have an assembled plant mitochondrial genome, you wish to evaluate the quality of the assembly
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
### examples
to assess mitogenome
 ```
himt assess -i input.fa/gfa -o output -c mitochondrial
```
to assess ptgenome
```
himt assess -i chloroplast
```

## Compare
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
examples
Supports Fasta and GFA files
```
himt compare -q genome1。fa/gfa -r genome2.fa/gfa -o output
```

