# HiMT
HiMT is a plant mitochondrial genome assembly toolkit written in Python3. It can extract high-depth reads from HiFi sequencing data and use Flye to assemble the mitochondria genome and chloroplasts genome of plants.
# installation
## conda(recommended)
```
conda install shuyuan_tang::himt -c bioconda
```
# Usage
```
usage: himt assemble [argument]

Assembling mitochondrail genome with HiFi sequcncing data

Required arguments:
  -i INPUT_FILE, --input_file INPUT_FILE
                        input a fasta or fastq file,gz compressed files are supported.
  -o OUTPUT_DIR, --output_dir OUTPUT_DIR
                        output directory.

options:
  -h, --help            show this help message and exit
  -n HEAD_NUMBER, --head_number HEAD_NUMBER
                        default=4,The number of kmer species randomly selected.
  -t THREAD, --thread THREAD
                        default=2 The number of thread used during code execution.
  -b {3,4}, --base_number {3,4}
                        default=3,only can be 3 and 4,The number of bases at the beginning of kmer.
  -fd FILTER_DEPTH, --filter_depth FILTER_DEPTH
                        read depths below this value will be filtered.
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
                        mitochondrial genome depth ranges between 15 and 50, input a value less than 0 (such as:-1) to
                        retain the maximum mitogenome depth.
```

## Assembly of plant mitogenomes

## general usage
```
himt assemble -i pineapple.LY.fa -o output -t 10
```
## Assembly of animal mitogenomes(Not mature)
```
himt assemble -i input.fa -o output -t 10 -s plant
```
