# EeayMT
EasyMT is a plant mitochondrial genome assembly toolkit written in Python3. It can extract high-depth reads from HiFi sequencing data and use Flye to assemble the mitochondria genome and chloroplasts genome of plants.
# installation
## conda(recommended)
```
conda install -c bioconda easymt
```
## docker image
```
docker pull dacongmian/easymt:1.0.0
docker run -it dacongmian/easymt:1.0.0
conda init
source .bashrc
conda activate easymt
```
# Usage
```
options:
  -h, --help            show this help message and exit
  -v, --version         show program's version number and exit

Required arguments:
  -i INPUT_FILE, --input_file INPUT_FILE
                        input a fasta or fastq file,gz compressed files are supported
  -o OUTPUT_DIR, --output_dir OUTPUT_DIR
                        output directory

Options arguments:
  -n HEAD_NUMBER, --head_number HEAD_NUMBER
                        default=4,The number of kmer species randomly selected
  -t THREAD, --thread THREAD
                        default=2 The number of thread used during code execution.
  -b {3,4}, --base_number {3,4}
                        default=3,only can be 3 and 4,The number of bases at the beginning of kmer
  -fd FLITER_DEPTH, --fliter_depth FLITER_DEPTH
                        read depths below this value will be filtered
  -fp FLITER_PERCENTAGE, --fliter_percentage FLITER_PERCENTAGE
                        default=0.3,The depth of the mitochondrial genome obtained by blast, the proportion adjusted downwards on this value.
  -p PROPORTION, --proportion PROPORTION
                        default=1,The percentage of the selected dataset from the entire file
  -s {plant,animal}, --species {plant,animal}
                        default=plant,Species category,only can be plant or animal
```
## example
## Assembly of plant mitogenomes
## fast modle
Randomly select a part of all reads and then extract high depth reads
```
easymt -i hifi_data.fasta -o output -t 12 -p 0.2
```
## general model
```
conda create -n env_name easymt
conda activate env_name
easymt -i hifi_data.fasta -o output -t 12
```
